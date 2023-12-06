#pragma once

#include"Layers/Layer.h"

#include"3D/VertexArray.h"
#include"3D/ShaderProgram.h"

#include"Quantum/Molecule.h"
#include"Maths/MarchingCubes.h"


#include<thread>
#include<mutex>

namespace SCESS
{
	class Camera;
	class GaussianCollection;
	struct HFResult;

	struct TriSurf
	{
		TriSurf(){}
		TriSurf(const std::vector<Triangle>& trianglesIn):triangles(trianglesIn) {}
		TriSurf(const TriSurf& tri):triangles(tri.triangles){}
		TriSurf(TriSurf&& tri) noexcept :triangles(std::move(tri.triangles)) {}

		void operator=(const TriSurf& tri)
		{
			triangles = tri.triangles;
		}

		void operator=(TriSurf&& tri) noexcept
		{
			triangles = std::move(tri.triangles);
		}

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


	private:

		void UploadTriSurface(VertexArray& vArray, const TriSurf& tri);
		void SelectMO(size_t index);

		sRef<VertexArray> mPosTriSurf, mNegTrisurf;
		std::mutex mTrisurfMutex;
		std::unique_ptr<std::jthread> mTrisurfThread;
		bool mTrisurfCalculated = false;
		TriSurf mPosSurfaceTri, mNegSurfaceTri;

		sRef<ShaderProgram> mDiffuseProg;
		sRef<HFResult> mHFResult;
		sRef<Molecule> mMolecule;

		std::unordered_map<uint32_t, uint32_t> mSelectedIndexToOrb;
		std::vector<std::string> mOrbLabels;


		float mIsosurfaceVal=0.02;
		int mIsosurfaceResolution = 20;
	};

}