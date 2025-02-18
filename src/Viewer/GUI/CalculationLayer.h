#pragma once
//
#include"GraphicsEngine/Layers/Layer.h"
#include"GraphicsEngine/3D/VertexArray.h"
#include"GraphicsEngine/3D/ShaderProgram.h"
#include"GraphicsEngine/Tools/Threading.h"
//
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/HF/HF.h"
#include "Quantum/CI/CI.h"
#include "Quantum/DFT/DFT.h"
#include "Quantum/CC/CC.h"
#include "Viewer/Isosurfaces/Isosurfaces.h"
//
#include <boost/circular_buffer.hpp>


namespace GEng {
    class Camera3D;
    class ShaderProgram;
    class VertexArray;
    class VertexBuffer;
}

namespace SCESS {
    class CGF;
}

namespace SCESS::Viewer {
    class TriSurf;
}


//
//
namespace SCESS::Viewer {
    //
    //
    enum class SCFMethod : int {
        None = -1,
        RHF = 0,
        UHF = 1,
        KS_DFT = 2
    };
    //
    enum class CorrelationMethod : int {
        None = -1,
        MP2 = 0,
        CCSD = 1,
        CCSD_T = 2,
        FCI = 3,
        FCIQMC = 4,
        CCQMC = 5
    };
    //
    enum class DFTGridRes : int {
        Low = 0,
        Medium,
        High,
        Ultrafine
    };
    //
    inline const char *SCFMethosToStr(SCFMethod method) {
        switch (method) {
            default:
            case SCFMethod::None: return "";
            case SCFMethod::RHF: return "RHF";
            case SCFMethod::UHF: return "UHF";
            case SCFMethod::KS_DFT: return "KS DFT";
        }
    }
    //
    inline const char *CorrelationMethodToStr(CorrelationMethod method) {
        switch (method) {
            default:
            case CorrelationMethod::None: return "";
            case CorrelationMethod::MP2: return "MP2";
            case CorrelationMethod::CCSD: return "CCSD";
            case CorrelationMethod::CCSD_T: return "CCSD(T)";
            case CorrelationMethod::FCI: return "FCI";
            case CorrelationMethod::CCQMC: return "CC";
            case CorrelationMethod::FCIQMC: return "FCI";
        }
    }
    //
    inline Quadrature::MolGrid::Resolution DFTGridResToRes(DFTGridRes res) {
        switch (res) {
            case DFTGridRes::Low:
                return Quadrature::MolGrid::LOW_RES;
            default:
            case DFTGridRes::Medium:
                return Quadrature::MolGrid::MEDIUM_RES;
            case DFTGridRes::High:
                return Quadrature::MolGrid::HIGH_RES;
            case DFTGridRes::Ultrafine:
                return Quadrature::MolGrid::ULTRAFINE_RES;
        }
    }
    //
    enum class MoleculeStyle:int {
        SpaceFill = 0,
        BallStick = 1
    };
    //
    //
    //
    class CalculationLayer : public GEng::Layer {
    public:
        CalculationLayer();

        virtual void OnInit() override;
        virtual void OnDestroy() override;
        virtual void OnUpdate(const sRef<GEng::Framebuffer> &screenBuffer, const GEng::Camera3D &cam,
                              const glm::mat4 &projectionMatrix, const glm::vec2 &nearFarPlanes) override;
        virtual void OnImGuiRender() override;
        virtual void OnEvent(GEng::Event &e) override;

    private:
        static void UploadTriSurface(GEng::VertexArray &vArray, const TriSurf &tri);
        void RenderMOIsosurface();
        void SelectMolecule();
        void UploadMolecule();
        bool SelectMoleculeFromFile();

        inline bool HaveSCFResult() {
            return mRHFWfn.has_value() || mUHFWfn.has_value() || mKSWfn.has_value();
        }

        inline bool HaveCorrelationResult() {
            return mCCSDWfn.has_value() || mCCSD_TWfn.has_value() || mFCIWfn.has_value()
                   || mMP2Energy.has_value() || mQMCEnergy != 0.0;
        }

        void ExecuteCalculation();
        void StopCalculation();
        //void DrawProgressScreen();
        void OrderOrbitals();
        void DrawResultsPanel();
        void DrawQuantumCalcPanel();
        void DrawMolSelectPanel();
        void DrawConsolePanel();
        void DrawQMCGraph();

        std::vector<std::string> mBasisSets;
        int mSelectedBasisSet = 0;

        int mSelectedMolecule = -1;
        uRef<Molecule> mMolecule;
        uRef<Molecule> mLastMolecule;
        BasisSet mMoleculeBasisSet;
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

        SCFMethod mSCFMethod = SCFMethod::RHF,
                mCalcSCFMethod = SCFMethod::None;
        CorrelationMethod mCorrelationMethod = CorrelationMethod::None,
                mCalcCorrelationMethod = CorrelationMethod::None;

        DFTGridRes mDFTGridRes = DFTGridRes::Medium;

        Integrals::MolInts mints;
        std::optional<HF::HF_Wfn> mRHFWfn;
        std::optional<HF::HF_Wfn> mUHFWfn;
        //Kohn-Sham DFT wavefunction
        std::optional<DFT::KS_Wfn> mKSWfn;
        bool mUseLastSCF = true, mInvalidateLastSCF = false;

        std::optional<CC::CCSD_Wfn> mCCSDWfn;
        std::optional<CC::CCSD_Wfn> mCCSD_TWfn;
        std::optional<double> mMP2Energy;

        std::optional<CI::CI_Wfn> mFCIWfn;

        std::string mCalcBasisSet;

        boost::circular_buffer<glm::vec2> mQMCGraphData;
        double mQMCEnergy = 0.0;

        std::unordered_map<uint32_t, uint32_t> mSelectedIndexToOrb;
        std::vector<std::string> mMOLabels;

        char mSelectedSpin = -1;
        int mSelectedMO = -1;

        float mIsosurfaceVal = 0.02f;
        int mIsosurfaceResolution = 20;
        float mIsosurfaceVolume = 5.0f;
    };
    //
    //
}
