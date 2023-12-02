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
						*std::pow(pos.y, g.mAngular[1])
						* std::pow(pos.z, g.mAngular[2])
						* std::exp(-std::abs(glm::dot(pos,pos)));
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
				* std::exp(-std::abs(glm::dot(pos, pos)));
		}
		return result;
	}

	TriSurf GraphLayer::GenerateIsosurface(const std::function<float(float, float, float)>& func,
		size_t resolution, float axis, float isosurfaceVal)
	{
		std::unique_ptr<MarchingCubes> marchingCubes = std::make_unique<MarchingCubes>();

		std::vector<std::vector<std::vector<float>>> scalarField(
			resolution,
			std::vector<std::vector<float>>(resolution, std::vector<float>(resolution)));

		for (size_t i = 0; i < resolution; i++)
		{
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

		TriSurf surface;
		marchingCubes->TriangulateField(scalarField, isosurfaceVal, &surface.triangles);

		//Rescale vertices
		for (size_t i = 0; i < surface.triangles.size(); i++)
		{
			for (uint8_t j = 0; j < 3; j++)
			{
				surface.triangles[i][j] = ((surface.triangles[i][j] * 2.0f / float(resolution)) - 1.0f) * axis;
			}
		}

		return surface;
	}

	sRef<VertexArray> GraphLayer::TriSurfaceToVertexArray(const TriSurf& tri)
	{
		sRef<VertexArray> vArray = VertexArray::Create();
		vArray->AddVertexBuffer(VertexBuffer::Create((const float*)(tri.triangles.data()),
			tri.triangles.size() * 3 * 3));

		return vArray;
	}

	void GraphLayer::OnAttach()
	{
		mDiffuseProg = ShaderProgram::CreateFromFile(R"(C:\dev\Chemistry\SCESS\Shaders\ColourDiffuse.glsl)");
		mDiffuseProg->Bind();
		mDiffuseProg->UploadVecF("uLightPos", {0.0f, 10.0f, 0.0f});
		mDiffuseProg->Unbind();
	}

	void GraphLayer::OnDetach()
	{
	}

	void GraphLayer::OnUpdate(const sRef<SCESS::Framebuffer>& screenBuffer, const Camera& cam,
		const glm::mat4& projectionMatrix, const glm::vec2& nearFarPlanes)
	{
		if (mPosTriSurf == nullptr || mNegTrisurf == nullptr)
		{
			return;
		}

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

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	void GraphLayer::SelectWaterWavefunction(size_t index)
	{
		auto waterWavefunction = [coeffs = mHFResult->Coeffs.col(index), water=mWaterMol](float x, float y, float z)->float {
			float result = 0.0f;
			for (size_t i = 0; i < water->GetOrbitals().size(); i++)
			{
				result += coeffs[i] * EvaluateGaussianFunction(water->GetOrbital(i), x, y, z);
			}
			return result;
			};

		mPosTriSurf = TriSurfaceToVertexArray(GenerateIsosurface(waterWavefunction, 40, 5.0f, 0.02f));
		mNegTrisurf = TriSurfaceToVertexArray(GenerateIsosurface(waterWavefunction, 40, 5.0f, -0.02f));
	}

	void GraphLayer::OnImGuiRender()
	{
		static std::vector<const char*> labels;

		if (!ImGui::Begin("Water Quantum Calc."))
		{
			ImGui::End();
			return;
		}

		if (ImGui::CollapsingHeader("Hartree Fock", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Button("Execute Hartree Fock"))
			{
				mWaterMol = std::make_shared<Molecule>(GenerateWaterMolecule());
				mHFResult = std::make_shared<HFResult>(SCESS::HF::ReducedHF(*mWaterMol, 1E-5, 500));
				
				std::unordered_map<double, uint32_t> sortMap;
				for (size_t i = 0; i < mHFResult->Energies.rows(); i++)
				{
					sortMap.insert({ mHFResult->Energies[i], i});
				}
				std::sort(mHFResult->Energies.begin(), mHFResult->Energies.end());
				mHFResult->Energies.reverseInPlace();

				Eigen::MatrixXd sortedCoeffs = Eigen::MatrixXd::Zero(mHFResult->Coeffs.rows(), mHFResult->Coeffs.rows());
				for (size_t i = 0; i < mHFResult->Energies.rows(); i++)
				{
					uint32_t index = sortMap.at(mHFResult->Energies[i]);
					sortedCoeffs.col(i) = mHFResult->Coeffs.col(index);

					std::stringstream ss;
					ss <<"MO "<< mHFResult->Energies.rows()-i << "; E=" << mHFResult->Energies[i] << " hartrees";
					mOrbLabels.push_back(ss.str());
				}
				mHFResult->Coeffs = sortedCoeffs;

				labels.clear();
				for (std::string& s : mOrbLabels)
				{
					labels.push_back(s.c_str());
				}
			}
			ImGui::Text("See console for Hartree Fock calculation progress...");
			
			if (mHFResult != nullptr)
			{
				if (ImGui::CollapsingHeader("Details", ImGuiTreeNodeFlags_DefaultOpen))
				{
					ImGui::Text("Total Hartree Fock Energy: %f hartrees", mHFResult->Energy);
				}
				if (ImGui::CollapsingHeader("Water Molecular Orbitals", ImGuiTreeNodeFlags_DefaultOpen))
				{
					ImGui::BulletText("Click on an MO to generate its isosurface!");
					int selectedOption = -1;
					ImGui::ListBox("", &selectedOption, labels.data(), labels.size());
					if (selectedOption != -1)
					{
						SelectWaterWavefunction(selectedOption);
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
