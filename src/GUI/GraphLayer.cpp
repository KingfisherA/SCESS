#include<pch.h>
#include "GraphLayer.h"
#include"imgui.h"

#include "Scene/Camera.hpp"
#include "3D/VertexArray.h"
#include"3D/ShaderProgram.h"
#include"Maths/Geometry.hpp"
#include "Maths/MarchingCubes.h"

#include "Quantum/Molecule.h"
#include"Quantum/Basis/Gaussian.h"
#include "Quantum/HF/HartreeFock.h"

#include"Application/Application.h"

#include "GUI/ImSpinner.h"

#include"Utilities/PathUtils.h"


namespace SCESS
{
	//static float OrbitalPZ(float x, float y, float z, const glm::vec3& pos) {
	//	x = x - pos.x;
	//	y = y - pos.y;
	//	z = z - pos.z;
	//	float r = std::sqrt(x * x + y * y + z * z);
	//	return z * std::exp(-r / 2) / std::sqrt(32 * std::numbers::pi);
	//}

	//static float OrbitalDZ2(float x, float y, float z)
	//{
	//	float r = std::sqrt(x * x + y * y + z * z);
	//	return (1.0 / (81 * std::sqrt(6 * std::numbers::pi))) * (3 * z * z - r * r) * std::exp(-r / 3);
	//}

	GraphLayer::~GraphLayer()
	{
	}

	std::function<float(float, float, float)> GraphLayer::GenerateGaussianFunction(const GaussianCollection& gc)
	{
		auto returnFunc = [gc](float x, float y, float z) ->float
			{
				float result = 0.0f;
				for (const Gaussian& g : gc.GetGaussians())
				{
					glm::vec3 pos = { x - g.mPosition.x(), y - g.mPosition.y(), z - g.mPosition.z() };
					result += g.mNorm * g.mCoeff *
						std::pow(pos.x, g.mAngular[0])
						* std::pow(pos.y, g.mAngular[1])
						* std::pow(pos.z, g.mAngular[2])
						* std::exp(-std::abs(glm::dot(pos, pos)));
				}
				return result;
			};

		return returnFunc;
	}

	float GraphLayer::EvaluateGaussianFunction(const GaussianCollection& gc, float x, float y, float z)
	{
		float result = 0.0f;
		for (const Gaussian& g : gc.GetGaussians())
		{
			glm::vec3 pos(
				x - (float)g.mPosition.x(),
				y - (float)g.mPosition.y(),
				z - (float)g.mPosition.z());

			result += g.mNorm * g.mCoeff *
				std::pow(pos.x, g.mAngular[0])
				* std::pow(pos.y, g.mAngular[1])
				* std::pow(pos.z, g.mAngular[2])
				* std::exp(-float(g.mExponent) * glm::dot(pos, pos));
		}
		return result;
	}

	void GraphLayer::GenerateIsosurface(const std::function<float(float, float, float)>& func,
		TriSurf* outSurface, size_t resolution, bool* threadShouldTerminate, float isosurfaceVol, float isosurfaceVal)
	{
		std::unique_ptr<MarchingCubes> marchingCubes = std::make_unique<MarchingCubes>();

		std::vector<std::vector<std::vector<float>>> scalarField(
			resolution,
			std::vector<std::vector<float>>(resolution, std::vector<float>(resolution)));

		for (size_t i = 0; i < resolution; i++)
		{
			if (threadShouldTerminate != nullptr && *threadShouldTerminate)
				return;
			for (size_t j = 0; j < resolution; j++)
			{
				for (size_t k = 0; k < resolution; k++)
				{
					glm::vec3 pos(
						(2.0 * float(i) / float(resolution) - 1.0) * isosurfaceVol,
						(2.0 * float(j) / float(resolution) - 1.0) * isosurfaceVol,
						(2.0 * float(k) / float(resolution) - 1.0) * isosurfaceVol);
					scalarField[i][j][k] = func(pos.x, pos.y, pos.z);
				}
			}
		}

		marchingCubes->TriangulateField(scalarField, isosurfaceVal, &outSurface->triangles, threadShouldTerminate);

		//Rescale vertices
		for (size_t i = 0; i < outSurface->triangles.size(); i++)
		{
			for (uint8_t j = 0; j < 3; j++)
			{
				outSurface->triangles[i][j] = ((outSurface->triangles[i][j] * 2.0f / float(resolution)) - 1.0f) * isosurfaceVol;
			}
		}
	}

	void GraphLayer::UploadTriSurface(VertexArray& vArray, const TriSurf& tri)
	{
		if (vArray.GetNoOfVBuffers() != 0)
		{
			const sRef<VertexBuffer>& vBuffer = vArray.GetVBuffer(0);
			vBuffer->StoreData((const float*)(tri.triangles.data()),
				tri.triangles.size() * 3 * 3);
			vArray.UpdateElementLength(tri.triangles.size() * 3);
		}
	}

	void GraphLayer::OnAttach()
	{
		mDiffuseProg = ShaderProgram::CreateFromFile(PathUtils::GetAssetPath(R"(Shaders\ColourDiffuse.glsl)").string());
		mDiffuseProg->Bind();
		mDiffuseProg->UploadVecF("uLightPos", { 0.0f, 10.0f, 0.0f });
		mDiffuseProg->Unbind();

		mPosTrisurfVAO = VertexArray::Create();
		mPosTrisurfVAO->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));
		mNegTrisurfVAO = VertexArray::Create();
		mNegTrisurfVAO->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));
	}

	void GraphLayer::OnDetach()
	{
	}

	void GraphLayer::OnUpdate(const sRef<SCESS::Framebuffer>& screenBuffer, const Camera& cam,
		const glm::mat4& projectionMatrix, const glm::vec2& nearFarPlanes)
	{
		glPolygonMode(GL_FRONT_AND_BACK, mRenderStyle);
		glDisable(GL_CULL_FACE);

		mDiffuseProg->Bind();

		mDiffuseProg->UploadMatF("uMPVMatrix", projectionMatrix * cam.GetViewMatrix());
		mDiffuseProg->UploadVecF("uLightPos", cam.GetActualPosition());

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

	void GraphLayer::RenderMOIsosurface()
	{
		const Eigen::MatrixXd* coeffsM = nullptr;

		if (mMethod == HFMethod::Restricted)
		{
			coeffsM = &mRHFResult->Coeffs;
		}
		else
		{
			coeffsM = mSelectedSpin == 1 ? &mUHFResult->CoeffsAlpha : &mUHFResult->CoeffsBeta;
		}

		std::function<float(float, float, float)> molecularWavefunction = [coeffs = coeffsM->col(mSelectedMO),
			&molecule = static_cast<const Molecule&>(*mMolecule)](float x, float y, float z) -> float
			{
				float result = 0.0f;
				for (size_t i = 0; i < molecule.GetOrbitals().size(); i++)
				{
					result += (float)coeffs[i] * EvaluateGaussianFunction(molecule.GetOrbital(i), x, y, z);
				}
				return result;
			};

		GenerateIsosurface(molecularWavefunction, &mPosTrisurf, mIsosurfaceResolution,
			&mTrisurfThread.ShouldTerminate, mIsosurfaceVolume, mIsosurfaceVal);
		GenerateIsosurface(molecularWavefunction, &mNegTrisurf, mIsosurfaceResolution,
			&mTrisurfThread.ShouldTerminate, mIsosurfaceVolume, -mIsosurfaceVal);
		mTrisurfThread.Finished = true;
	}

	constexpr const char* PROG_SCREEN_NAME = "HF Progress";

	void GraphLayer::RunHF()
	{
		if (mTrisurfThread.Ptr != nullptr)
		{
			mTrisurfThread.ShouldTerminate = true;
			mTrisurfThread.Ptr = nullptr;
		}

		mHartreeFockThread.Ptr = nullptr;

		mMolecule = std::make_shared<Molecule>(SCESS::GenerateWaterMolecule());
		mHartreeFockThread.Finished = false;
		mHFConsoleOutput.str(std::string());
		mHartreeFockThread.Ptr = std::make_unique<std::jthread>([this]() {
			if (mMethod == HFMethod::Restricted)
				mRHFResult = std::make_shared<RHFResult>(SCESS::HF::RestrictedHF(*mMolecule, 1E-5, 100, &mHFConsoleOutput));
			else
				mUHFResult = std::make_shared<UHFResult>(SCESS::HF::UnrestrictedHF(*mMolecule, 1E-5, 100, &mHFConsoleOutput));

			mHartreeFockThread.Finished = true;
			});
	}

	void GraphLayer::DrawProgressScreen()
	{
		ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Once, ImVec2(0.5f, 0.5f));
		ImGui::SetNextWindowSize({ 600,300 }, ImGuiCond_Once);
		if (ImGui::BeginPopupModal(PROG_SCREEN_NAME))
		{
			ImSpinner::SpinnerAtom("##atomSpinner", 20.0, 2.5);
			ImGui::SameLine();
			ImGui::Text("Running RHF single-point energy. Please wait...");

			const std::string& consoleStr = mHFConsoleOutput.str();
			ImGui::InputTextMultiline("##hfConsoleOutput", (char*)consoleStr.c_str(), consoleStr.length(),
				{ 0,0 }, ImGuiInputTextFlags_ReadOnly);

			if (mHartreeFockThread.Finished)
				ImGui::CloseCurrentPopup();

			ImGui::EndPopup();
		}
	}

	void GraphLayer::SortMOLabels(std::vector<const char*>& moLabels)
	{
		auto sortCoeffs = [this](Eigen::VectorXd& Energies, Eigen::MatrixXd& Coeffs) {
			std::unordered_map<double, uint32_t> sortMap;
			for (size_t i = 0; i < Energies.rows(); i++)
			{
				sortMap.insert({ Energies[i], i });
			}
			std::sort(Energies.begin(), Energies.end());
			Energies.reverseInPlace();

			Eigen::MatrixXd sortedCoeffs = Eigen::MatrixXd::Zero(Coeffs.rows(), Coeffs.rows());

			mMOLabels.reserve(Energies.rows());
			for (size_t i = 0; i < Energies.rows(); i++)
			{
				uint32_t index = sortMap.at(Energies[i]);
				sortedCoeffs.col(i) = Coeffs.col(index);

				std::stringstream ss;
				ss.precision(3);
				ss << "MO " << Energies.rows() - i << "; E=" << Energies[i] << " hartrees";
				mMOLabels.push_back(ss.str());
			}
			Coeffs = sortedCoeffs;
			};

		mMOLabels.clear();
		if (mMethod == HFMethod::Restricted)
			sortCoeffs(mRHFResult->Energies, mRHFResult->Coeffs);
		else
		{
			sortCoeffs(mUHFResult->EnergiesAlpha, mUHFResult->CoeffsAlpha);
			sortCoeffs(mUHFResult->EnergiesBeta, mUHFResult->CoeffsBeta);
		}

		moLabels.clear();
		moLabels.reserve(mMOLabels.size());
		for (std::string& s : mMOLabels)
		{
			moLabels.push_back(s.c_str());
		}
	}

	void GraphLayer::OnImGuiRender()
	{
		static std::vector<const char*> moLabels;

		DrawProgressScreen();
		if (mHartreeFockThread.Ptr != nullptr && !mHartreeFockThread.Finished)
			ImGui::OpenPopup(PROG_SCREEN_NAME);

		if (!ImGui::Begin("Hartree Fock"))
		{
			ImGui::End();
			return;
		}

		ImGui::BulletText("Selected Molecule: Water");
		ImGui::Separator();
		ImGui::SetCursorPosX((ImGui::GetWindowWidth() - 100.0) * 0.5);
		if (ImGui::Button("Run RHF", { 100,0 }))
		{
			mMethod = HFMethod::Restricted;
			RunHF();
		}

		ImGui::SetCursorPosX((ImGui::GetWindowWidth() - 100.0) * 0.5);
		if (ImGui::Button("Run UHF", { 100,0 }))
		{
			mMethod = HFMethod::Unrestricted;
			RunHF();
		}

		if (mHartreeFockThread.Finished)
		{
			if (mHartreeFockThread.Ptr != nullptr)
			{
				SortMOLabels(moLabels);
				mHartreeFockThread.Ptr = nullptr;
			}

			if (ImGui::CollapsingHeader("Hartree Fock Console", ImGuiTreeNodeFlags_DefaultOpen))
			{
				const std::string& consoleStr = mHFConsoleOutput.str();
				ImGui::InputTextMultiline("##hfConsoleOutput", (char*)consoleStr.c_str(), consoleStr.length(),
					{ 0,0 }, ImGuiInputTextFlags_ReadOnly);
			}

			if (ImGui::CollapsingHeader("Result", ImGuiTreeNodeFlags_DefaultOpen))
			{
				ImGui::BulletText("Basis Set: STO-3G");
				ImGui::BulletText("Total Hartree Fock Energy: %.4f hartrees",
					mMethod == HFMethod::Restricted ? mRHFResult->Energy : mUHFResult->Energy);

				ImGui::Indent();
				if (ImGui::CollapsingHeader("Orbitals", ImGuiTreeNodeFlags_DefaultOpen))
				{
					if (mSelectedMO != -1)
					{
						if (mMethod == HFMethod::Restricted)
						{
							ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "Previewing MO %i", mRHFResult->Energies.rows() - mSelectedMO);
							ImGui::BulletText("Energy: %.3f hartrees", mRHFResult->Energies(mSelectedMO));
							ImGui::BulletText("Occupation: %i e-", mRHFResult->Occupation(mSelectedMO));
							if (mSelectedMO + 1 < mRHFResult->Occupation.rows())
								if (mRHFResult->Occupation(mSelectedMO + 1) == 0 
									&& mRHFResult->Occupation(mSelectedMO) > 0)
									ImGui::BulletText("HOMO");
							if (mSelectedMO - 1 > 0)
								if (mRHFResult->Occupation(mSelectedMO - 1) > 0 
									&& mRHFResult->Occupation(mSelectedMO) == 0)
									ImGui::BulletText("LUMO");
						}
						else
						{
							const auto& energies =
								mSelectedSpin == 1 ? mUHFResult->EnergiesAlpha : mUHFResult->EnergiesBeta;
							const auto& occupation = mSelectedSpin == 1 ?
								mUHFResult->OccupationAlpha
								: mUHFResult->OccupationBeta;
							ImGui::TextColored({ 1.0f, 1.0f, 0.0f, 1.0f }, "Previewing MO %i", energies.rows() - mSelectedMO);
							ImGui::BulletText("Spin: %s", mSelectedSpin == 1 ? "Alpha (up)" : "Beta (down)");
							ImGui::BulletText("Energy: %.3f hartrees", energies(mSelectedMO));
							ImGui::BulletText("Occupation: %i e-", occupation(mSelectedMO));
							if (mSelectedMO + 1 < occupation.rows())
								if (occupation(mSelectedMO + 1) == 0 && occupation(mSelectedMO)>0)
									ImGui::BulletText("HOMO");
							if (mSelectedMO - 1 > 0)
								if (occupation(mSelectedMO - 1) > 0 && occupation(mSelectedMO) == 0)
									ImGui::BulletText("LUMO");
						}
						ImGui::Separator();
					}

					int selectedOption = -1;
					int selectedOption2 = -1;
					if (mTrisurfThread.Ptr == nullptr)
					{
						ImGui::Text("Select an MO to generate an isosurface");
						if (mMethod == HFMethod::Restricted)
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
						mTrisurfThread.ShouldTerminate = false;
						mTrisurfThread.Finished = false;
						mIsosurfaceUploaded = false;

						if (mMethod == HFMethod::Restricted)
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
						mTrisurfThread.Ptr = std::make_unique<std::jthread>(
							std::bind(&GraphLayer::RenderMOIsosurface, this));
					}

					if (mTrisurfThread.Finished && !mIsosurfaceUploaded)
					{
						if (mPosTrisurf.triangles.size() != 0)
							UploadTriSurface(*mPosTrisurfVAO, mPosTrisurf);
						if (mPosTrisurf.triangles.size() != 0)
							UploadTriSurface(*mNegTrisurfVAO, mNegTrisurf);

						mIsoWasNull = mPosTrisurf.triangles.size() == 0 && mNegTrisurf.triangles.size() == 0;

						mPosTrisurf.triangles.clear();
						mNegTrisurf.triangles.clear();
						mIsosurfaceUploaded = true;

						mTrisurfThread.Ptr = nullptr;
					}

					if (mTrisurfThread.Ptr != nullptr)
					{
						ImSpinner::SpinnerTwinBall("spinner1", 16.0, 16.0f, 3.0, 6.0, { 1.0f, 0.0f, 0.0f, 1.0f }, ImSpinner::white, 3.0f, 1);
						ImGui::SameLine();
						ImGui::Text("Generating isosurface...");
						if (ImGui::Button("Cancel", { 100,0 }))
						{
							mTrisurfThread.ShouldTerminate = true;
							mTrisurfThread.Ptr = nullptr;
							mTrisurfThread.Finished = false;
						}
					}
					else
					{
						if (mIsosurfaceUploaded && mTrisurfThread.Finished)
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
						}
					}

				}
			}
		}

		ImGui::End();
	}

	void GraphLayer::OnEvent(Event& e)
	{

	}
}
