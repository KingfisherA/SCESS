#pragma once
//
#include "Dahlia/Layers/Layer.h"
#include "Dahlia/Graphics/GraphicsTypes.h"
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
#include <vector>
#include <boost/circular_buffer.hpp>
#include <glm/vec2.hpp>


namespace DH {
    class Buffer;
    class Camera3D;
    class CommandList;
    class GPUMesh;
    class GraphicsPipeline;
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
        //targetLayout: attachment formats of the framebuffer the scene is drawn into
        explicit CalculationLayer(const DH::RenderTargetLayout &targetLayout);

        virtual void OnInit() override;
        virtual void OnDestroy() override;
        virtual void OnImGuiRender() override;
        virtual void OnEvent(DH::Event &e) override;

        //Records the scene's draws. Must be called inside a render pass targeting a framebuffer with the layout
        //given to the constructor.
        void Render(DH::CommandList &cmd, const DH::Camera3D &cam, const glm::mat4 &projectionMatrix);

    private:
        //A vertex buffer whose contents are replaced from the CPU. Uploads are queued and applied by Flush() at the
        //start of Render(), before any draws are recorded, so a frame's draws always match its buffer contents.
        struct DynamicVertexBuffer {
            sRef<DH::Buffer> buffer;
            std::vector<float> pending;
            bool dirty = false;
            uint32_t elementCount = 0;
            uint32_t floatsPerElement = 0;

            void Queue(const float *data, size_t floatCount);
            void Flush();
        };

        void CreatePipelines();
        static void UploadTriSurface(DynamicVertexBuffer &vBuffer, const TriSurf &tri);
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
        sRef<DH::GPUMesh> mAtomMesh;
        //Per atom instance data, see UploadMolecule()
        DynamicVertexBuffer mMoleculeBuffer;
        sRef<DH::GPUMesh> mBondMesh;
        MoleculeStyle mMolRenderStyle = MoleculeStyle::BallStick;
        bool mRenderMolecule = true;
        float mMoleculeOpacity = 1.0f;

        DynamicVertexBuffer mPosTrisurfBuffer, mNegTrisurfBuffer;
        bool mIsosurfaceUploaded = false, mIsoWasNull = false;
        TriSurf mNegTrisurf, mPosTrisurf;
        float mIsoOpacity = 1.0f;

        DH::RenderTargetLayout mTargetLayout;
        sRef<DH::GraphicsPipeline> mMoleculePipeline;
        sRef<DH::GraphicsPipeline> mBondPipeline;
        //Polygon mode is pipeline state, so the wireframe isosurface needs its own pipeline
        sRef<DH::GraphicsPipeline> mIsosurfacePipeline, mIsosurfaceWirePipeline;

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
