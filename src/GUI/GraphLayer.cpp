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


namespace SCESS
{
	static float OrbitalPZ(float x, float y, float z, const glm::vec3& pos) {
		x = x - pos.x;
		y = y - pos.y;
		z = z - pos.z;
		float r = std::sqrt(x * x + y * y + z * z);
		return z * std::exp(-r / 2) / std::sqrt(32 * std::numbers::pi);
	}

	static float OrbitalDZ2(float x, float y, float z)
	{
		float r = std::sqrt(x * x + y * y + z * z);
		return (1.0 / (81 * std::sqrt(6 * std::numbers::pi))) * (3 * z * z - r * r) * std::exp(-r / 3);
	}


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
			glm::vec3 pos = { x - g.mPosition.x(), y - g.mPosition.y(), z - g.mPosition.z() };
			result += g.mNorm * g.mCoeff *
				std::pow(pos.x, g.mAngular[0])
				* std::pow(pos.y, g.mAngular[1])
				* std::pow(pos.z, g.mAngular[2])
				* std::exp(-g.mExponent * glm::dot(pos, pos));
		}
		return result;
	}

	void GraphLayer::GenerateIsosurface(const std::function<float(float, float, float)>& func,
		TriSurf* outSurface, size_t resolution, bool* threadShouldTerminate, float axis, float isosurfaceVal)
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
						(2.0 * float(i) / float(resolution) - 1.0) * axis,
						(2.0 * float(j) / float(resolution) - 1.0) * axis,
						(2.0 * float(k) / float(resolution) - 1.0) * axis);
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
				outSurface->triangles[i][j] = ((outSurface->triangles[i][j] * 2.0f / float(resolution)) - 1.0f) * axis;
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
		mDiffuseProg = ShaderProgram::CreateFromFile(R"(C:\dev\Chemistry\SCESS\Shaders\ColourDiffuse.glsl)");
		mDiffuseProg->Bind();
		mDiffuseProg->UploadVecF("uLightPos", { 0.0f, 10.0f, 0.0f });
		mDiffuseProg->Unbind();

		mPosTriSurf = VertexArray::Create();
		mPosTriSurf->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));
		mNegTrisurf = VertexArray::Create();
		mNegTrisurf->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));
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
		mPosTriSurf->Bind();
		mPosTriSurf->SubmitArrays();
		mPosTriSurf->Unbind();

		mDiffuseProg->UploadVecF("uColour", { 0.0f, 0.0f, 1.0f });
		mNegTrisurf->Bind();
		mNegTrisurf->SubmitArrays();
		mNegTrisurf->Unbind();


		mDiffuseProg->Unbind();

		glEnable(GL_CULL_FACE);
	}

	void GraphLayer::SelectMO(size_t index)
	{
		auto molecularWavefunction = [coeffs = mHFResult->Coeffs.col(index),
			&molecule = static_cast<const Molecule&>(*mMolecule)](float x, float y, float z) -> float
			{
				float result = 0.0f;
				for (size_t i = 0; i < molecule.GetOrbitals().size(); i++)
				{
					result += coeffs[i] * EvaluateGaussianFunction(molecule.GetOrbital(i), x, y, z);
				}
				return result;
			};

		GenerateIsosurface(molecularWavefunction, &mPosSurfaceTri, mIsosurfaceResolution,
			&mTrisurfShouldTerminate, 10.0f, mIsosurfaceVal);
		GenerateIsosurface(molecularWavefunction, &mNegSurfaceTri, mIsosurfaceResolution,
			&mTrisurfShouldTerminate, 10.0f, -mIsosurfaceVal);
		mTrisurfFinished = true;
	}

	constexpr const char* PROG_SCREEN_NAME = "HF Progress";

	void GraphLayer::RunHF()
	{
		if (mTrisurfThread != nullptr)
		{
			mTrisurfShouldTerminate = true;
			mTrisurfThread = nullptr;
		}

		if (mHFThread != nullptr)
			mHFThread->join();

		mMolecule = std::make_shared<Molecule>(SCESS::GenerateCOMolecule());
		mHFThreadFinished = false;
		mHFConsoleOutput.str(std::string());
		mHFThread = std::make_unique<std::jthread>([this]() {
			mHFResult = std::make_shared<HFResult>(SCESS::HF::RestrictedHF(*mMolecule, 1E-5, 100, &mHFConsoleOutput));
			mHFThreadFinished = true;
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

			if (mHFThreadFinished)
				ImGui::CloseCurrentPopup();

			ImGui::EndPopup();
		}
	}

	void GraphLayer::SortMOLabels(std::vector<const char*>& moLabels)
	{
		std::unordered_map<double, uint32_t> sortMap;
		for (size_t i = 0; i < mHFResult->Energies.rows(); i++)
		{
			sortMap.insert({ mHFResult->Energies[i], i });
		}
		std::sort(mHFResult->Energies.begin(), mHFResult->Energies.end());
		mHFResult->Energies.reverseInPlace();

		Eigen::MatrixXd sortedCoeffs = Eigen::MatrixXd::Zero(mHFResult->Coeffs.rows(), mHFResult->Coeffs.rows());
		mOrbLabels.clear();
		mOrbLabels.reserve(mHFResult->Energies.rows());
		for (size_t i = 0; i < mHFResult->Energies.rows(); i++)
		{
			uint32_t index = sortMap.at(mHFResult->Energies[i]);
			sortedCoeffs.col(i) = mHFResult->Coeffs.col(index);

			std::stringstream ss;
			ss << "MO " << mHFResult->Energies.rows() - i << "; E=" << mHFResult->Energies[i] << " hartrees";
			mOrbLabels.push_back(ss.str());
		}
		mHFResult->Coeffs = sortedCoeffs;

		moLabels.clear();
		moLabels.reserve(mOrbLabels.size());
		for (std::string& s : mOrbLabels)
		{
			moLabels.push_back(s.c_str());
		}
	}

	void GraphLayer::OnImGuiRender()
	{
		static std::vector<const char*> moLabels;

		DrawProgressScreen();
		if (mHFThread != nullptr && !mHFThreadFinished)
			ImGui::OpenPopup(PROG_SCREEN_NAME);

		if (!ImGui::Begin("Hartree Fock"))
		{
			ImGui::End();
			return;
		}

		ImGui::BulletText("Selected Molecule: Carbon Monoxide");
		ImGui::Separator();
		ImGui::SetCursorPosX((ImGui::GetWindowWidth() - 100.0) * 0.5);
		if (ImGui::Button("Run RHF", { 100,0 }))
			RunHF();

		if (mHFThreadFinished)
		{
			if (mHFThread != nullptr)
			{
				SortMOLabels(moLabels);
				mHFThread = nullptr;
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
				ImGui::BulletText("Total Hartree Fock Energy: %f hartrees", mHFResult->Energy);

				ImGui::Indent();
				if (ImGui::CollapsingHeader("Orbitals", ImGuiTreeNodeFlags_DefaultOpen))
				{
					int selectedOption = -1;
					if (mTrisurfThread == nullptr)
					{
						ImGui::Text("Select an MO to generate an isosurface");
						ImGui::ListBox("##listbox1", &selectedOption, moLabels.data(), moLabels.size());
					}

					if (selectedOption != -1 && mTrisurfThread == nullptr)
					{
						mTrisurfShouldTerminate = false;
						mTrisurfFinished = false;
						mTrisurfThread = std::make_unique<std::jthread>(
							std::bind(&GraphLayer::SelectMO, this, selectedOption));
					}

					if (mTrisurfFinished)
					{
						if (mPosSurfaceTri.triangles.size() != 0
							|| mNegSurfaceTri.triangles.size() != 0)
						{
							mTrisurfThread = nullptr;
							UploadTriSurface(*mPosTriSurf, mPosSurfaceTri);
							UploadTriSurface(*mNegTrisurf, mNegSurfaceTri);
							mPosSurfaceTri.triangles.clear();
							mNegSurfaceTri.triangles.clear();
						}
					}

					if (mTrisurfThread != nullptr)
					{
						ImSpinner::SpinnerTwinBall("spinner1", 16.0, 16.0f, 3.0, 6.0, { 1.0f, 0.0f, 0.0f, 1.0f }, ImSpinner::white, 3.0f, 1);
						ImGui::SameLine();
						ImGui::Text("Generating isosurface...");
						if (ImGui::Button("Cancel", { 100,0 }))
						{
							mTrisurfShouldTerminate = true;
							mTrisurfThread = nullptr;
						}
					}
					else
					{
						if (mTrisurfFinished)
							ImGui::TextColored({ 0.0f,1.0f,0.0f,1.0f }, "Isosurface generated.");

						if (ImGui::CollapsingHeader("Isosurface Options", ImGuiTreeNodeFlags_DefaultOpen))
						{
							ImGui::Text("Isosurface Value:");
							ImGui::DragFloat("##isosurfVal", &mIsosurfaceVal, 0.001, 1E-10, 1.0f);
							ImGui::Text("Isosurface Resolution:");
							ImGui::DragInt("##isoResolutionVal", &mIsosurfaceResolution, 1.0f, 10, 1000);
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
