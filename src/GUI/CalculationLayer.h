#pragma once

#include"GraphicsEngine/Layers/Layer.h"
#include"GraphicsEngine/3D/VertexArray.h"
#include"GraphicsEngine/3D/ShaderProgram.h"
#include"GraphicsEngine/Utilities/Threading.h"

#include"Chemistry/Quantum/Molecule.h"
#include "Chemistry/Quantum/CI/ConfigurationInteraction.h"
#include "Chemistry/Quantum/DFT/DensityFunctional.h"
#include "Maths/Isosurfaces.h"


namespace GEng {
    class Camera3D;
    class ShaderProgram;
    class VertexArray;
    class VertexBuffer;
}

namespace SCESS {
    struct CCSD_TWavefunction;
    class ContractedGaussian;
    struct RHFWavefunction;
    struct UHFWavefunction;
    struct KSWavefunction;
    struct CCSDWavefunction;
    class TriSurf;

    enum class SCFMethod:int {
        RHF = 0,
        UHF = 1,
        DFT = 2
    };

    enum class CorrelationMethod:int {
        None = 0,
        MP2 = 1,
        CCSD = 2,
        CCSD_T = 3,
        FCI = 4
    };

    inline const char *CorrelationMethodToStr(CorrelationMethod method) {
        switch (method) {
            case CorrelationMethod::None: return "";
            case CorrelationMethod::MP2: return "MP2";
            case CorrelationMethod::CCSD: return "CCSD";
            case CorrelationMethod::CCSD_T: return "CCSD(T)";
            case CorrelationMethod::FCI: return "FCI";
            default: return "";
        }
    }

    enum class MoleculeStyle:int {
        SpaceFill = 0,
        BallStick = 1
    };

    class CalculationLayer : public GEng::Layer {
    public:
        CalculationLayer() : Layer("Graph Layer") {
        }

        //virtual ~CalculationLayer();

        virtual void OnInit() override;
        virtual void OnDestroy() override;
        virtual void OnUpdate(const sRef<GEng::Framebuffer> &screenBuffer, const GEng::Camera3D &cam,
                              const glm::mat4 &projectionMatrix, const glm::vec2 &nearFarPlanes) override;
        virtual void OnImGuiRender() override;
        virtual void OnEvent(GEng::Event &e) override;

    private:
        void UploadTriSurface(GEng::VertexArray &vArray, const TriSurf &tri);
        void RenderMOIsosurface();
        void SelectMolecule();
        void UploadMolecule();
        bool SelectMoleculeFromFile();

        glm::mat3 GenerateBondRotMatrix(const glm::vec3 &bondDirNorm, const glm::vec3 &startingDirNorm);

        void ExecuteCalculation();
        void DrawProgressScreen();
        void OrderOrbitals();
        void DrawResultsPanel();
        void DrawQuantumCalcPanel();
        void DrawMolSelectPanel();
        void DrawConsolePanel();

        std::vector<std::string> mBasisSets;
        int mSelectedBasisSet = 0;

        int mSelectedMolecule = -1;
        uRef<Molecule> mMolecule;
        uRef<Molecule> mLastMolecule;
        sRef<GEng::VertexArray> mMoleculeVAO;
        sRef<GEng::VertexBuffer> mMoleculeBuffer;
        sRef<GEng::VertexArray> mMoleculeBondVAO;
        MoleculeStyle mMolRenderStyle = MoleculeStyle::BallStick;
        bool mRenderMolecule = true;
        float mMoleculeOpacity = 1.0f;

        sRef<GEng::VertexArray> mPosTrisurfVAO, mNegTrisurfVAO;
        bool mIsosurfaceUploaded = false, mIsoWasNull = false;
        TriSurf mNegTrisurf, mPosTrisurf;
        float mIsoOpacity = 1.0f;

        sRef<GEng::ShaderProgram> mDiffuseProg;
        sRef<GEng::ShaderProgram> mMoleculeDiffuseProg;
        sRef<GEng::ShaderProgram> mSolidColourProg;

        std::stringstream mConsoleOutput;
        bool mOrbitalsOrdered = false;
        GEng::Thread mCalculationThread;

        SCFMethod mSCFMethod = SCFMethod::RHF;
        CorrelationMethod mCorrelationMethod = CorrelationMethod::None;

        uRef<RHFWavefunction> mRHFWfn;
        uRef<UHFWavefunction> mUHFWfn;
        //Kohn-Sham DFT wavefunction
        uRef<KSWavefunction> mKSWfn;
        bool mUseLastSCF = true, mInvalidateLastSCF = false;

        uRef<CCSDWavefunction> mCCSDWfn;
        uRef<CCSD_TWavefunction> mCCSD_TWfn;
        std::optional<double> mMP2Energy;

        uRef<CIWavefunction> mFCIWfn;

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
