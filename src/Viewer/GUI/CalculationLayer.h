#pragma once
//
#include "Dahlia/Layers/Layer.h"
#include "Dahlia/3D/VertexArray.h"
#include "Dahlia/3D/ShaderProgram.h"
#include "Dahlia/Tools/Threading.h"
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
#include <atomic>
#include <mutex>
#include <streambuf>
#include <boost/circular_buffer.hpp>
#include <glm/vec2.hpp>


namespace DH {
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
    //A string buffer which may be written to by the calculation thread while the UI thread reads it. It has no put
    //area so that every write goes through the locked overflow/xsputn.
    class SyncStringBuf : public std::streambuf {
    public:
        std::string Snapshot() const {
            std::lock_guard<std::mutex> lock(mMutex);
            return mData;
        }

        void Clear() {
            std::lock_guard<std::mutex> lock(mMutex);
            mData.clear();
        }

    protected:
        int_type overflow(int_type ch) override {
            if (!traits_type::eq_int_type(ch, traits_type::eof())) {
                std::lock_guard<std::mutex> lock(mMutex);
                mData.push_back(traits_type::to_char_type(ch));
            }
            return traits_type::not_eof(ch);
        }

        std::streamsize xsputn(const char *s, std::streamsize n) override {
            std::lock_guard<std::mutex> lock(mMutex);
            mData.append(s, static_cast<size_t>(n));
            return n;
        }

    private:
        mutable std::mutex mMutex;
        std::string mData;
    };
    //
    //DH::Thread::Stop() blocks until the worker joins. This lets the UI request a stop and keep rendering; the
    //worker still reports completion with MarkAsFinished() and is joined by the next Start() or the destructor.
    class CalculationThread : public DH::Thread {
    public:
        void RequestStop() { mThreadRequestStop.store(true); }

        [[nodiscard]] bool IsStopping() const { return IsRunning() && ShouldStop(); }
    };
    //
    //
    //
    class CalculationLayer : public DH::Layer {
    public:
        CalculationLayer();

        virtual void OnInit() override;
        virtual void OnDestroy() override;
        virtual void OnUpdate(const sRef<DH::Framebuffer> &screenBuffer, const DH::Camera3D &cam,
                              const glm::mat4 &projectionMatrix, const glm::vec2 &nearFarPlanes) override;
        virtual void OnImGuiRender() override;
        virtual void OnEvent(DH::Event &e) override;

    private:
        static void UploadTriSurface(DH::VertexArray &vArray, const TriSurf &tri);
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
        BasisSet mMoleculeBasisSet;
        sRef<DH::VertexArray> mMoleculeVAO;
        sRef<DH::VertexBuffer> mMoleculeBuffer;
        sRef<DH::VertexArray> mMoleculeBondVAO;
        MoleculeStyle mMolRenderStyle = MoleculeStyle::BallStick;
        bool mRenderMolecule = true;
        float mMoleculeOpacity = 1.0f;

        sRef<DH::VertexArray> mPosTrisurfVAO, mNegTrisurfVAO;
        bool mIsosurfaceUploaded = false, mIsoWasNull = false;
        TriSurf mNegTrisurf, mPosTrisurf;
        float mIsoOpacity = 1.0f;

        sRef<DH::ShaderProgram> mDiffuseProg;
        sRef<DH::ShaderProgram> mMoleculeDiffuseProg;
        sRef<DH::ShaderProgram> mSolidColourProg;

        //Calculation output; the stream writes into mConsoleBuffer which the UI thread reads with Snapshot()
        SyncStringBuf mConsoleBuffer;
        std::stringstream mConsoleOutput;
        std::string mConsoleSnapshot;
        bool mOrbitalsOrdered = false;
        CalculationThread mCalculationThread;

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
        bool mUseLastSCF = true;
        //Written by the calculation thread when it finishes
        std::atomic<bool> mInvalidateLastSCF = false;

        std::optional<CC::CCSD_Wfn> mCCSDWfn;
        std::optional<CC::CCSD_Wfn> mCCSD_TWfn;
        std::optional<double> mMP2Energy;

        std::optional<CI::CI_Wfn> mFCIWfn;

        std::string mCalcBasisSet;

        //Projected energies pushed by the calculation thread and plotted by the UI thread
        std::mutex mQMCGraphMutex;
        boost::circular_buffer<glm::vec2> mQMCGraphData;
        double mQMCEnergy = 0.0;

        std::unordered_map<uint32_t, uint32_t> mSelectedIndexToOrb;
        std::vector<std::string> mMOLabels;

        char mSelectedSpin = -1;
        int mSelectedMO = -1;

        float mIsosurfaceVal = 0.02f;
        int mIsosurfaceResolution = 20;
        //The field is allocated as resolution^3 floats
        static constexpr int MAX_ISOSURFACE_RESOLUTION = 250;
        float mIsosurfaceVolume = 5.0f;
    };
    //
    //
}
