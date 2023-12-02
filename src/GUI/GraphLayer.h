#pragma once

#include"Layers/Layer.h"

#include"3D/VertexArray.h"
#include"3D/ShaderProgram.h"

#include"Quantum/Molecule.h"

namespace SCESS
{
	class Camera;
	struct Triangle;
	class GaussianCollection;
	struct HFResult;

	struct TriSurf
	{
		std::vector<Triangle> triangles;
	};

	class GraphLayer : public Layer
	{
	public:
		GraphLayer():Layer("Graph Layer") {}
		virtual ~GraphLayer();

		static float EvaluateGaussianFunction(const GaussianCollection& gc, float x, float y, float z);
		static std::function<float(float, float, float)> GenerateGaussianFunction(const GaussianCollection& gc);

		virtual void OnAttach() override;
		virtual void OnDetach() override;
		virtual void OnUpdate(const sRef<SCESS::Framebuffer>& screenBuffer, const Camera& cam,
			const glm::mat4& projectionMatrix, const glm::vec2& nearFarPlanes) override;
		virtual void OnImGuiRender() override;
		virtual void OnEvent(Event& e) override;

		static TriSurf GenerateIsosurface(const std::function<float(float, float, float)>& func,
			size_t resolution, float axis = 10.0f, float isosurfaceVal = 0.01f);
		sRef<VertexArray> TriSurfaceToVertexArray(const TriSurf& tri);

		void SelectWaterWavefunction(size_t index);

	private:
		sRef<VertexArray> mPosTriSurf, mNegTrisurf;
		sRef<ShaderProgram> mDiffuseProg;
		sRef<HFResult> mHFResult;
		sRef<Molecule> mWaterMol;

		std::unordered_map<uint32_t, uint32_t> mSelectedIndexToOrb;
		std::vector<std::string> mOrbLabels;
	};

}