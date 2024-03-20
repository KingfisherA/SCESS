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

	enum class CalcMethod
	{
		RHF,
		UHF,
		DFT
	};

	enum class MoleculeStyle:int
	{
		SpaceFill =0,
		BallStick = 1 
	};

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
		void UploadMolecule();
		void SelectMoleculeFromFile();

		glm::mat3 GenerateBondRotMatrix(const glm::vec3& bondDirNorm, const glm::vec3& startingDirNorm);

		void ExecuteCalculation();
		void DrawProgressScreen();
		void OrderOrbitals();

		std::vector<std::string> mBasisSets;
		int mSelectedBasisSet = 0;

		int mSelectedMolecule = -1;
		std::string mMoleculeName;
		uRef<Molecule> mMolecule;
		sRef<VertexArray> mMoleculeVAO;
		sRef<VertexBuffer> mMoleculeBuffer;
		sRef<VertexArray> mMoleculeBondVAO;
		MoleculeStyle mMolRenderStyle = MoleculeStyle::BallStick;
		bool mRenderMolecule = true;
		float mMoleculeOpacity = 1.0f;

		sRef<VertexArray> mPosTrisurfVAO, mNegTrisurfVAO;
		bool mIsosurfaceUploaded = false, mIsoWasNull = false;
		TriSurf mNegTrisurf, mPosTrisurf;
		float mIsoOpacity = 1.0f;

		sRef<ShaderProgram> mDiffuseProg;
		sRef<ShaderProgram> mMoleculeDiffuseProg;
		sRef<ShaderProgram> mSolidColourProg;

		std::stringstream mConsoleOutput;
		bool mOrbitalsOrdered = false;
		Thread mCalculationThread;

		CalcMethod mMethod = CalcMethod::RHF;
		sRef<RHFResult> mRHFResult;
		sRef<UHFResult> mUHFResult;
		double mDFTEnergy = 0.0;
		double mMP2Energy = 0.0;
		std::string mCalcBasiSet;

		std::unordered_map<uint32_t, uint32_t> mSelectedIndexToOrb;
		std::vector<std::string> mMOLabels;

		char mSelectedSpin = -1;
		int mSelectedMO = -1;

		float mIsosurfaceVal = 0.02f;
		int mIsosurfaceResolution = 20;
		float mIsosurfaceVolume = 5.0f;
	};

}
