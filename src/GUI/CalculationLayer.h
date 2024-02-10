#pragma once

#include"Layers/Layer.h"

#include"3D/VertexArray.h"
#include"3D/ShaderProgram.h"

#include"Quantum/Molecule.h"
#include "Maths/Isosurfaces.h"

#include"Utilities/Threading.h"

namespace SCESS
{
	class Camera;
	class GaussianCollection;
	struct RHFResult;
	struct UHFResult;
	class TriSurf;


	class CalculationLayer : public Layer
	{
	public:
		CalculationLayer() :Layer("Graph Layer") {}
		//virtual ~CalculationLayer();

		virtual void OnAttach() override;
		virtual void OnDetach() override;
		virtual void OnUpdate(const sRef<SCESS::Framebuffer>& screenBuffer, const Camera& cam,
			const glm::mat4& projectionMatrix, const glm::vec2& nearFarPlanes) override;
		virtual void OnImGuiRender() override;
		virtual void OnEvent(Event& e) override;

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
		bool mIsosurfaceUploaded = false, mIsoWasNull = false;
		TriSurf mNegTrisurf, mPosTrisurf;

		sRef<ShaderProgram> mDiffuseProg;

		std::stringstream mHFConsoleOutput;
		bool mOrderedLabels = false;
		Thread mHartreeFockThread;

		enum class HFMethod
		{
			Restricted,
			Unrestricted
		} mMethod = HFMethod::Restricted;
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