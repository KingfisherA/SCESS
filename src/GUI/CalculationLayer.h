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
	struct RHFResult;
	struct UHFResult;

	struct TriSurf
	{
		TriSurf() {}
		TriSurf(const std::vector<Triangle>& trianglesIn) :triangles(trianglesIn) {}
		TriSurf(const TriSurf& tri) :triangles(tri.triangles) {}
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

	class CalculationLayer : public Layer
	{
	public:
		CalculationLayer() :Layer("Graph Layer") {}
		//virtual ~CalculationLayer();

		static float EvaluateGaussianFunction(const GaussianCollection& gc, float x, float y, float z);
		static void EvaluateScalarField(const std::function<float(float, float, float)>& func, 
			std::vector<std::vector<std::vector<float>>>& scalarField, float isosurfaceVol, bool* threadShouldTerminate);
		static std::function<float(float, float, float)> GenerateGaussianFunction(const GaussianCollection& gc);
		std::function<float(float, float, float)> GenerateMolecularWavefunction();


		virtual void OnAttach() override;
		virtual void OnDetach() override;
		virtual void OnUpdate(const sRef<SCESS::Framebuffer>& screenBuffer, const Camera& cam,
			const glm::mat4& projectionMatrix, const glm::vec2& nearFarPlanes) override;
		virtual void OnImGuiRender() override;
		virtual void OnEvent(Event& e) override;

		static void GenerateIsosurface(const std::function<float(float, float, float)>& func, TriSurf* outSurface,
			size_t resolution, bool* threadShouldTerminate = nullptr, float isosurfaceVol = 10.0f, float isosurfaceVal = 0.01f);

	private:
		void UploadTriSurface(VertexArray& vArray, const TriSurf& tri);
		void RenderMOIsosurface();
		void SelectMolecule();

		void RunHF();
		void DrawProgressScreen();
		void SortMOLabels(std::vector<const char*>& moLabels);

		int mSelectedMolecule = -1;
		std::string mMoleculeName;
		uRef<Molecule> mMolecule;

		sRef<VertexArray> mPosTrisurfVAO, mNegTrisurfVAO;
		struct TrisurfaceThread
		{
			std::unique_ptr<std::jthread> Ptr;
			bool ShouldTerminate, Finished = false;
		} mTrisurfThread;
		bool mIsosurfaceUploaded, mIsoWasNull = false;
		TriSurf mPosTrisurf, mNegTrisurf;

		sRef<ShaderProgram> mDiffuseProg;

		std::stringstream mHFConsoleOutput;
		struct HartreeFockThread
		{
			std::unique_ptr<std::jthread> Ptr;
			bool Finished = false;
		} mHartreeFockThread;

		enum class HFMethod
		{
			Restricted,
			Unrestricted
		} mMethod;
		sRef<RHFResult> mRHFResult;
		sRef<UHFResult> mUHFResult;
		double mMP2Energy = 0.0;

		std::unordered_map<uint32_t, uint32_t> mSelectedIndexToOrb;
		std::vector<std::string> mMOLabels;

		char mSelectedSpin = -1;
		int mSelectedMO = -1;

		float mIsosurfaceVal = 0.02;
		int mIsosurfaceResolution = 20;
		float mIsosurfaceVolume = 5.0f;
	};

}