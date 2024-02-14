#include<pch.h>
#include"CalculationLayer.h"
#include"imgui.h"

#include"Scene/Camera.hpp"
#include"3D/VertexArray.h"
#include"3D/ShaderProgram.h"
#include"Maths/Geometry.hpp"
#include"Maths/Isosurfaces.h"

#include"Quantum/Molecule.h"
#include"Quantum/Basis/Gaussian.h"
#include"Quantum/HF/HartreeFock.h"
#include"Quantum/DFT/DFT.h"

#include"Application/Application.h"

#include"GUI/ImSpinner.h"

#include"Utilities/PathUtils.h"
#include"CalculationLayer.h"




namespace SCESS
{
	void CalculationLayer::UploadTriSurface(VertexArray& vArray, const TriSurf& tri)
	{
#ifdef _DEBUG
		if (!tri.HasFinished())
		{
			throw std::exception("Trisurface generation thread has not completed.");
		}
#endif

		if (vArray.GetNoOfVBuffers() != 0)
		{
			const sRef<VertexBuffer>& vBuffer = vArray.GetVBuffer(0);
			vBuffer->StoreData((const float*)(tri.GetTriangles().data()),
				tri.GetTriangles().size() * 3 * 3);
			vArray.UpdateElementLength(tri.GetTriangles().size() * 3);
		}
	}

	void CalculationLayer::OnAttach()
	{
		mDiffuseProg = ShaderProgram::CreateFromFile(PathUtils::GetAssetPath(R"(Shaders\ColourDiffuse.glsl)").string());
		mDiffuseProg->Bind();
		mDiffuseProg->UploadVecF("uLightPos", { 0.0f, 10.0f, 0.0f });
		mDiffuseProg->Unbind();

		mPosTrisurfVAO = VertexArray::Create();
		mPosTrisurfVAO->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));
		mNegTrisurfVAO = VertexArray::Create();
		mNegTrisurfVAO->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));

		mSphereVAO = VertexArray::Create();
		sRef<Mesh> sphereMesh = Geometry::CreateQuadsphere(1.0f, 30);
		mSphereVAO->AddVertexBuffer(VertexBuffer::Create(sphereMesh->GetVertexArray(),
			sphereMesh->GetVertexFloatCount()));
		mSphereVAO->AddIndexBuffer(IndexBuffer::Create(sphereMesh->GetIndiciesArray(),
			sphereMesh->GetIndexIntCount()));
	}

	void CalculationLayer::OnDetach()
	{
	}

	void CalculationLayer::OnUpdate(const sRef<SCESS::Framebuffer>& screenBuffer, const Camera& cam,
		const glm::mat4& projectionMatrix, const glm::vec2& nearFarPlanes)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDisable(GL_CULL_FACE);

		glm::mat4 pvMatrix = projectionMatrix * cam.GetViewMatrix();

		mDiffuseProg->Bind();
		mDiffuseProg->UploadVecF("uLightPos", cam.GetActualPosition());

		if (mRenderMolecule)
		{
			mDiffuseProg->UploadF("uOpacity", mMoleculeOpacity);
			mDiffuseProg->UploadMatF("uMPVMatrix", pvMatrix);

			if (mMolecule != nullptr)
			{
				mSphereVAO->Bind();
				for (const Atom& a : mMolecule->GetAtoms())
				{
					const AtomStyle& style = AtomStyle::GetStyle(a.GetAtomSymbol());
					mDiffuseProg->UploadMatF("uMPVMatrix", pvMatrix * glm::translate(glm::vec3{ a.GetPos().x(), a.GetPos().y() ,a.GetPos().z() }) *
						glm::scale(glm::vec3{ style.mRadius }));
					mDiffuseProg->UploadVecF("uColour", style.mColour);


					mSphereVAO->SubmitArrays();
				}
				mSphereVAO->Unbind();
			}
		}

		glPolygonMode(GL_FRONT_AND_BACK, mRenderStyle);
		mDiffuseProg->UploadF("uOpacity", mIsoOpacity);
		mDiffuseProg->UploadMatF("uMPVMatrix", pvMatrix);

		mDiffuseProg->UploadVecF("uColour", { 1.0f, 0.0f, 0.0f });
		mPosTrisurfVAO->Bind();
		mPosTrisurfVAO->SubmitArrays();
		mPosTrisurfVAO->Unbind();

		mDiffuseProg->UploadVecF("uColour", { 0.0f, 0.0f, 1.0f });
		mNegTrisurfVAO->Bind();
		mNegTrisurfVAO->SubmitArrays();
		mNegTrisurfVAO->Unbind();

		mDiffuseProg->Unbind();

		glEnable(GL_CULL_FACE);
	}

	void CalculationLayer::RenderMOIsosurface()
	{
		const Eigen::MatrixXd* coeffsM = nullptr;

		mPosTrisurf.Stop();
		mNegTrisurf.Stop();

		if (mMethod == CalcMethod::RHF)
		{
			coeffsM = &mRHFResult->Coeffs;
		}
		else
		{
			coeffsM = mSelectedSpin == 1 ? &mUHFResult->CoeffsAlpha : &mUHFResult->CoeffsBeta;
		}

		std::shared_ptr<ScalarField> field = std::make_shared<ScalarField>();
		field->Initilise(mIsosurfaceResolution,
			GenerateMolecularWavefunction<float>(coeffsM->col(mSelectedMO), *mMolecule),
			mIsosurfaceVolume);
		field->Evaluate(true);

		mPosTrisurf.Initilise(field, mIsosurfaceVolume, mIsosurfaceVal);
		mNegTrisurf.Initilise(field, mIsosurfaceVolume, -mIsosurfaceVal);

		mPosTrisurf.Evaluate(true);
		mNegTrisurf.Evaluate(true);
	}

	constexpr const char* PROG_SCREEN_NAME = "Calc_Progress";

	void CalculationLayer::SelectMolecule()
	{
		static SCESS::Molecule(*moleculeFuncs[])() =
		{
			{&GenerateWaterMolecule},
			{&GenerateAmmoniaMolecule},
			{&GenerateCOMolecule},
			{&GenerateEtheneMolecule},
			{&GenerateMethaneMolecule},
			{&GenerateLithiumHydrideMolecule},
			{&GenerateHydrogenMolecule}
		};

		mMolecule = std::make_unique<SCESS::Molecule>(moleculeFuncs[mSelectedMolecule]());
	}

	void CalculationLayer::ExecuteCalculation()
	{
		if (mSelectedMolecule == -1)
			return;

		mPosTrisurf.Reset();
		mNegTrisurf.Reset();

		mOrderedLabels = false;
		mSelectedMO = -1;

		mMolecule->MakeOrbitals();
		mConsoleOutput.str(std::string());

		mCalculationThread.Start([this]() {
			switch (mMethod)
			{
			case CalcMethod::RHF:
				mRHFResult = std::make_shared<RHFResult>(
					std::move(SCESS::HF::RestrictedHF(*mMolecule, 1E-5, 100, &mConsoleOutput)));
				break;
			case CalcMethod::UHF:
				mUHFResult = std::make_shared<UHFResult>(
					std::move(SCESS::HF::UnrestrictedHF(*mMolecule, 1E-5, 100, &mConsoleOutput)));
				break;
			case CalcMethod::DFT:
				mDFTEnergy = SCESS::DFT::RunDFT(*mMolecule, &mConsoleOutput);
				break;
			}
			mCalculationThread.MarkAsFinished();
			});
	}

	void CalculationLayer::DrawProgressScreen()
	{
		ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Once, ImVec2(0.5f, 0.5f));
		ImGui::SetNextWindowSize({ 600,300 }, ImGuiCond_Once);
		if (ImGui::BeginPopupModal(PROG_SCREEN_NAME))
		{
			ImSpinner::SpinnerAtom("##atomSpinner", 20.0, 2.5);
			ImGui::SameLine();
			ImGui::Text("Running %s single-point energy. Please wait...", [this]() {
				switch (mMethod)
				{
				case CalcMethod::RHF:
					return "RHF";
				case CalcMethod::UHF:
					return "UHF";
				case CalcMethod::DFT:
					return "DFT";
				};
				}());


			const std::string& consoleStr = mConsoleOutput.str();
			ImGui::InputTextMultiline("##consoleOutput", (char*)consoleStr.c_str(), consoleStr.length(),
				{ 0,0 }, ImGuiInputTextFlags_ReadOnly);

			if (mCalculationThread.HasFinished())
				ImGui::CloseCurrentPopup();

			ImGui::EndPopup();
		}
	}

	void CalculationLayer::SortMOLabels(std::vector<const char*>& moLabels)
	{
		auto sortCoeffs = [this](Eigen::VectorXd& Energies, Eigen::MatrixXd& Coeffs,
			const Eigen::VectorXi& Occupation) {
				std::unordered_map<double, uint32_t> sortMap;
				for (size_t i = 0; i < Energies.rows(); i++)
				{
					sortMap.insert({ Energies[i], i });
				}
				std::sort(Energies.begin(), Energies.end());

				Eigen::MatrixXd sortedCoeffs = Eigen::MatrixXd::Zero(Coeffs.rows(), Coeffs.rows());

				mMOLabels.reserve(Energies.rows());
				for (size_t i = 0; i < Energies.rows(); i++)
				{
					uint32_t index = sortMap.at(Energies[i]);
					sortedCoeffs.col(i) = Coeffs.col(index);

					std::stringstream ss;
					ss.precision(3);
					ss << "MO " << i + 1 << "; E=" << Energies[i] << " hartrees. ";

					if (i + 1 < Occupation.rows())
						if (Occupation(i + 1) == 0 && Occupation(i) > 0)
							ss << "HOMO";
					if (i != 0)
						if (Occupation(i - 1) > 0 && Occupation(i) == 0)
							ss << "LUMO";

					mMOLabels.push_back(ss.str());
				}
				Coeffs = sortedCoeffs;
			};

		mMOLabels.clear();
		if (mMethod == CalcMethod::RHF)
			sortCoeffs(mRHFResult->Energies, mRHFResult->Coeffs, mRHFResult->Occupation);
		else
		{
			sortCoeffs(mUHFResult->EnergiesAlpha, mUHFResult->CoeffsAlpha, mUHFResult->OccupationAlpha);
			sortCoeffs(mUHFResult->EnergiesBeta, mUHFResult->CoeffsBeta, mUHFResult->OccupationBeta);
		}

		moLabels.clear();
		moLabels.reserve(mMOLabels.size());
		for (std::string& s : mMOLabels)
		{
			moLabels.push_back(s.c_str());
		}
	}

	void CalculationLayer::OnImGuiRender()
	{
		static std::vector<const char*> moLabels;
		static const char* moleculeNames[] =
		{
			"Water",
			"Ammonia",
			"Carbon Monoxide",
			"Ethene",
			"Methane",
			"Lithium Hydride",
			"Hydrogen"
		};

		if (ImGui::Begin("Select Molecule"))
		{
			ImGui::SeparatorText("Molecules");
			ImGui::BulletText("Selected Molecule: %s", mMoleculeName.c_str());
			int selectedMolecule = -1;
			ImGui::ListBox("##selectMoleculeLB", &selectedMolecule, moleculeNames, IM_ARRAYSIZE(moleculeNames));
			ImGui::Checkbox("Draw Molecule", &mRenderMolecule);
			if (mRenderMolecule)
				ImGui::SliderFloat("##molOpacity", &mMoleculeOpacity, 0.0f, 1.0f, "%.2f");


			if (selectedMolecule != -1)
			{
				mSelectedMolecule = selectedMolecule;
				mMoleculeName = moleculeNames[mSelectedMolecule];
				SelectMolecule();
			}

			//if (ImGui::Button("Test"))
			//{
			//	mIsosurfaceVolume = 25.0f;
			//	std::shared_ptr<ScalarField> field = std::make_shared<ScalarField>();
			//	field->Initilise(50,
			//		[](float x, float y, float z)->float {
			//			return (3*z*z-x*x-y*y) * exp(-0.5 * sqrt(x * x + y * y + z * z));
			//		},
			//		mIsosurfaceVolume);
			//	field->Evaluate(false);

			//	mPosTrisurf.Initilise(field, mIsosurfaceVolume, mIsosurfaceVal);
			//	mNegTrisurf.Initilise(field, mIsosurfaceVolume, -mIsosurfaceVal);

			//	mPosTrisurf.Evaluate(false);
			//	mNegTrisurf.Evaluate(false);

			//	UploadTriSurface(*mPosTrisurfVAO, mPosTrisurf);
			//	UploadTriSurface(*mNegTrisurfVAO, mNegTrisurf);
			//	mPosTrisurf.Reset();
			//	mNegTrisurf.Reset();
			//}

		}
		ImGui::End();

		DrawProgressScreen();
		if (mCalculationThread.IsRunning())
			ImGui::OpenPopup(PROG_SCREEN_NAME);


		if (ImGui::Begin("Quantum Calculation"))
		{
			if (mSelectedMolecule != -1)
			{
				ImGui::SeparatorText("Hartree Fock");
				ImGui::SetCursorPosX((ImGui::GetWindowWidth() - 100.0) * 0.5);
				if (ImGui::Button("Run RHF", { 100,0 }))
				{
					mMethod = CalcMethod::RHF;
					ExecuteCalculation();
				}

				ImGui::SetCursorPosX((ImGui::GetWindowWidth() - 100.0) * 0.5);
				if (ImGui::Button("Run UHF", { 100,0 }))
				{
					mMethod = CalcMethod::UHF;
					ExecuteCalculation();
				}

				if (mRHFResult != nullptr)
				{
					ImGui::SeparatorText("Moller-Plesset");

					ImGui::SetCursorPosX((ImGui::GetWindowWidth() - 200.0) * 0.5);
					if (ImGui::Button("Calculate MP2 Energy", { 200,0 }))
					{
						if (mMethod == CalcMethod::RHF)
						{
							mMP2Energy = HF::CalculateMP2Energy(*mRHFResult);
						}
					}

					if (mMP2Energy != 0)
					{
						ImGui::Text("Moller-Plesset Order 2 Correction= %.8f hartrees", mMP2Energy);
						ImGui::TextColored({ 0.0, 1.0, 0.0, 1.0 }, "Total Energy (HF+MP2) = %.8f hartrees", mMP2Energy + mRHFResult->Energy);
					}
				}

				ImGui::SeparatorText("Density Functional Theory");
				ImGui::SetCursorPosX((ImGui::GetWindowWidth() - 100.0) * 0.5);
				if (ImGui::Button("Run DFT", { 100,0 }))
				{
					mMethod = CalcMethod::DFT;
					ExecuteCalculation();
				}
			}

			if (mCalculationThread.HasFinished())
			{
				if (ImGui::CollapsingHeader("Console", ImGuiTreeNodeFlags_DefaultOpen))
				{
					const std::string& consoleStr = mConsoleOutput.str();
					ImGui::InputTextMultiline("##consoleOutput", (char*)consoleStr.c_str(), consoleStr.length(),
						{ 0,0 }, ImGuiInputTextFlags_ReadOnly);
				}

				if (ImGui::CollapsingHeader("Result", ImGuiTreeNodeFlags_DefaultOpen))
				{
					if (mMethod == CalcMethod::RHF || mMethod == CalcMethod::UHF)
					{
						if (!mOrderedLabels)
						{
							SortMOLabels(moLabels);
							mOrderedLabels = true;
						}

						ImGui::BulletText("Basis Set: STO-3G");
						ImGui::BulletText("Total Hartree Fock Energy: %.4f hartrees",
							mMethod == CalcMethod::RHF ? mRHFResult->Energy : mUHFResult->Energy);

						ImGui::Indent();
						if (ImGui::CollapsingHeader("Orbitals", ImGuiTreeNodeFlags_DefaultOpen))
						{
							if (mSelectedMO != -1)
							{
								if (mMethod == CalcMethod::RHF)
								{
									ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "Previewing MO %i", mSelectedMO + 1);
									if (mSelectedMO + 1 < mRHFResult->Occupation.rows())
										if (mRHFResult->Occupation(mSelectedMO + 1) == 0
											&& mRHFResult->Occupation(mSelectedMO) > 0)
											ImGui::TextColored({ 1.0f, 0.0f, 0.0f, 1.0f }, "HOMO");
									if (mSelectedMO - 1 > 0)
										if (mRHFResult->Occupation(mSelectedMO - 1) > 0
											&& mRHFResult->Occupation(mSelectedMO) == 0)
											ImGui::TextColored({ 0.0f, 0.5f, 0.5f, 1.0f }, "LUMO");
									ImGui::BulletText("Energy: %.3f hartrees", mRHFResult->Energies(mSelectedMO));
									ImGui::BulletText("Occupation: %i e-", mRHFResult->Occupation(mSelectedMO));
								}
								else
								{
									const auto& energies =
										mSelectedSpin == 1 ? mUHFResult->EnergiesAlpha : mUHFResult->EnergiesBeta;
									const auto& occupation = mSelectedSpin == 1 ?
										mUHFResult->OccupationAlpha
										: mUHFResult->OccupationBeta;
									ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "Previewing MO %i", mSelectedMO + 1);
									if (mSelectedMO + 1 < occupation.rows())
										if (occupation(mSelectedMO + 1) == 0 && occupation(mSelectedMO) > 0)
											ImGui::TextColored({ 1.0f, 0.0f, 0.0f, 1.0f }, "HOMO");
									if (mSelectedMO - 1 > 0)
										if (occupation(mSelectedMO - 1) > 0 && occupation(mSelectedMO) == 0)
											ImGui::TextColored({ 0.0f, 0.5f, 0.5f, 1.0f }, "LUMO");
									ImGui::BulletText("Spin: %s", mSelectedSpin == 1 ? "Alpha (up)" : "Beta (down)");
									ImGui::BulletText("Energy: %.3f hartrees", energies(mSelectedMO));
									ImGui::BulletText("Occupation: %i e-", occupation(mSelectedMO));

								}
								ImGui::Separator();
							}

							int selectedOption = -1;
							int selectedOption2 = -1;
							if (!(mPosTrisurf.IsRunning() || mNegTrisurf.IsRunning()))
							{
								ImGui::Text("Select an MO to generate an isosurface");
								if (mMethod == CalcMethod::RHF)
									ImGui::ListBox("##listbox1", &selectedOption, &moLabels[0], moLabels.size());
								else
								{
									float width = ImGui::GetColumnWidth() * 0.5f;

									ImGui::Text("Alpha (spin up) Orbitals");
									ImGui::SameLine(width, 10.0f);
									ImGui::Text("Beta (spin down) Orbitals");
									ImGui::PushItemWidth(width);
									ImGui::ListBox("##listbox1", &selectedOption, &moLabels[0], mUHFResult->EnergiesAlpha.rows());
									ImGui::SameLine();
									ImGui::ListBox("##listbox2", &selectedOption2, &moLabels[mUHFResult->EnergiesAlpha.rows()], mUHFResult->EnergiesBeta.rows());
									ImGui::PopItemWidth();
								}
							}

							if (selectedOption != -1 || selectedOption2 != -1)
							{
								if (mMethod == CalcMethod::RHF)
									mSelectedMO = selectedOption;
								else
								{
									if (selectedOption != -1)
									{
										mSelectedMO = selectedOption;
										mSelectedSpin = 1;
									}
									else
									{
										mSelectedMO = selectedOption2;
										mSelectedSpin = -1;
									}

								}

								mIsosurfaceUploaded = false;
								mIsoWasNull = false;
								RenderMOIsosurface();
							}

							if ((mPosTrisurf.HasFinished() && mNegTrisurf.HasFinished()) && !mIsosurfaceUploaded)
							{
								//if (mPosTrisurf.triangles.size() != 0)
								UploadTriSurface(*mPosTrisurfVAO, mPosTrisurf);
								//if (mNegTrisurf.triangles.size() != 0)
								UploadTriSurface(*mNegTrisurfVAO, mNegTrisurf);

								mIsoWasNull = mPosTrisurf.GetTriangles().size() == 0 && mNegTrisurf.GetTriangles().size() == 0;

								mPosTrisurf.Reset();
								mNegTrisurf.Reset();
								mIsosurfaceUploaded = true;
							}

							if (mPosTrisurf.IsRunning() || mNegTrisurf.IsRunning())
							{
								ImSpinner::SpinnerTwinBall("spinner1", 16.0, 16.0f, 3.0, 6.0, { 1.0f, 0.0f, 0.0f, 1.0f }, ImSpinner::white, 3.0f, 1);
								ImGui::SameLine();
								ImGui::Text("Generating isosurface...");
								if (ImGui::Button("Cancel", { 100,0 }))
								{
									mPosTrisurf.Reset();
									mNegTrisurf.Reset();
								}
							}
							else
							{
								if (mIsosurfaceUploaded)
								{
									if (!mIsoWasNull)
										ImGui::TextColored({ 0.0f,1.0f,0.0f,1.0f }, "Isosurface generated.");
									else
										ImGui::TextColored({ 1.0f,0.0f,0.0f,1.0f }, "Wavefunction has no valid isosurface in the given parameters.");
								}

								if (ImGui::CollapsingHeader("Isosurface Options", ImGuiTreeNodeFlags_DefaultOpen))
								{
									ImGui::Text("Isosurface Value:");
									ImGui::DragFloat("##isosurfVal", &mIsosurfaceVal, 0.001, 1E-10, 1.0f);
									ImGui::Text("Isosurface Resolution:");
									ImGui::DragInt("##isoResolutionVal", &mIsosurfaceResolution, 1.0f, 10, 1000);
									ImGui::Text("Isosurface Volume:");
									ImGui::DragFloat("bohr##isoVolumeVal", &mIsosurfaceVolume, 0.5f, 0.01f, 100.0f);
									ImGui::Text("Isosurface Opacity:");
									ImGui::DragFloat("##isoOpacity", &mIsoOpacity, 0.0f, 0.01f, 1.0f);

								}
							}
						}
					}
					else if (mMethod == CalcMethod::DFT)
					{
						ImGui::TextColored({ 0.0, 1.0, 0.0, 1.0 }, "Density Functional calculation finished!");
						ImGui::BulletText("Total DFT Energy : %.4f hartrees", mDFTEnergy);
					}
				}
			}
		}
		ImGui::End();
	}

	void CalculationLayer::OnEvent(Event& e)
	{
	}
}
