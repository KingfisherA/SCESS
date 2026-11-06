#include "CalculationLayer.h"
//
#include <numeric>
//
#include <ImGui/imgui.h>
#include <nfd.h>
#ifndef GLM_ENABLE_EXPERIMENTAL // Defined by the Dahlia target
#define GLM_ENABLE_EXPERIMENTAL
#endif
#include <glm/gtx/transform.hpp>
#include<implot.h>
//
#include "AtomDraw.h"
#include "Dahlia/Scene/Camera3D.h"
#include "Dahlia/Graphics/Buffer.h"
#include "Dahlia/Graphics/CommandList.h"
#include "Dahlia/Graphics/GPUMesh.h"
#include "Dahlia/Graphics/GraphicsPipeline.h"
#include "Dahlia/Graphics/ShaderModule.h"
#include "Dahlia/Maths/Geometry.h"
#include "Viewer/GUI/ImSpinner.h"
//
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/Molecule.h"
#include "Quantum/Basis/CGF.h"
#include "Quantum/HF/HF.h"
#include "Quantum/DFT/DFT.h"
#include "Quantum/CC/CC.h"
#include "Quantum/CI/CI.h"
#include "Quantum/F12/F12.h"
//
#include "Viewer/Isosurfaces/Isosurfaces.h"
//
#include "Viewer/Application/Application.h"
#include "Resources/AssetFolder.h"
#include "Resources/ShaderSources.h"
#include "Quantum/QMC/QMC.h"


//
//
namespace SCESS::Viewer {
    using namespace DH;

    //Push constant blocks; each matches the block declared in both stages of its shader (std140 layout)
    namespace {
        //Resources/Shaders/MoleculeDiffuse.*
        struct MoleculePushConstants {
            glm::mat4 pvMatrix;
            glm::vec4 lightPos;
            glm::vec4 params; //x = opacity
        };

        //Resources/Shaders/IsosurfaceDiffuse.*
        struct IsosurfacePushConstants {
            glm::mat4 pvMatrix;
            glm::vec4 lightPos;
            glm::vec4 colour; //a = opacity
        };

        //Resources/Shaders/SolidColour.*
        struct SolidColourPushConstants {
            glm::mat4 mvpMatrix;
            glm::vec4 colour; //a = opacity
        };

        //Floats per atom in the molecule instance buffer: (scale) (translation xyz) (colour rgb)
        constexpr uint32_t MOLECULE_INSTANCE_FLOATS = 7;
        //Vertex buffer binding for the molecule instance data. GPUMesh::Draw binds the mesh's own attributes to
        //bindings 0-2, so this must lie outside that range.
        constexpr uint32_t MOLECULE_INSTANCE_BINDING = 3;
    }

    void CalculationLayer::DynamicVertexBuffer::Queue(const float *data, size_t floatCount) {
        pending.assign(data, data + floatCount);
        dirty = true;
    }

    void CalculationLayer::DynamicVertexBuffer::Flush() {
        if (!dirty)
            return;
        dirty = false;

        const size_t bytes = pending.size() * sizeof(float);
        elementCount = (uint32_t) (pending.size() / floatsPerElement);
        if (bytes == 0)
            return;

        if (buffer == nullptr)
            buffer = Buffer::Create({bytes, BufferUsage::Vertex, BufferMemory::Dynamic}, pending.data());
        else if (buffer->GetSize() != bytes)
            buffer->Resize(bytes, pending.data());
        else
            buffer->Update(pending.data(), bytes);
    }

    void CalculationLayer::UploadTriSurface(DynamicVertexBuffer &vBuffer, const TriSurf &tri) {
#ifdef _DEBUG
        if (!tri.HasFinished()) {
            throw std::runtime_error("Trisurface generation thread has not completed.");
        }
#endif

        vBuffer.Queue((const float *) (tri.GetTriangles().data()), tri.GetTriangles().size() * 3 * 3);
    }

    CalculationLayer::CalculationLayer(const RenderTargetLayout &targetLayout): Layer("Graph Layer"),
        mTargetLayout(targetLayout), mQMCGraphData(1000) {
        mConsoleOutput.std::ios::rdbuf(&mConsoleBuffer);

        mMoleculeBuffer.floatsPerElement = MOLECULE_INSTANCE_FLOATS;
        //One vertex is three floats
        mPosTrisurfBuffer.floatsPerElement = 3;
        mNegTrisurfBuffer.floatsPerElement = 3;
    }

    void CalculationLayer::OnInit() {
        ImPlot::CreateContext();

        mAtomMesh = GPUMesh::Create(Geometry::CreateQuadsphere(1.0f, 30));
        mBondMesh = GPUMesh::Create(Geometry::CreateCylinder(0.2f, 0.7f, 20));

        CreatePipelines();

        //Orbital basis sets exclude the CABS family; every basis set may be used as a CABS, CABS family first
        for (const BasisSetInfo &info: BasisLibrary::List()) {
            if (BasisLibrary::IsCABS(info))
                mCABSSets.insert(mCABSSets.begin() + std::ranges::count_if(mCABSSets, BasisLibrary::IsCABS), info);
            else {
                mBasisSets.push_back(info);
                mCABSSets.push_back(info);
            }
        }
        const BasisSetInfo *sto3gInfo = BasisLibrary::Find("sto-3g");
        const auto sto3g = std::ranges::find_if(mBasisSets, [sto3gInfo](const BasisSetInfo &info) {
            return sto3gInfo && info.name == sto3gInfo->name;
        });
        mSelectedBasisSet = sto3g == mBasisSets.end() ? 0 : (int) std::distance(mBasisSets.begin(), sto3g);
        SelectDefaultCABS();
    }

    void CalculationLayer::SelectDefaultCABS() {
        if (!mCABSFollowsOBS || mBasisSets.empty())
            return;

        const std::string defaultCABS = BasisLibrary::DefaultCABS(mBasisSets[mSelectedBasisSet].name);
        const auto it = std::ranges::find_if(mCABSSets, [&defaultCABS](const BasisSetInfo &info) {
            return info.name == defaultCABS;
        });
        //Without a default keep the previous choice
        if (it != mCABSSets.end())
            mSelectedCABS = (int) std::distance(mCABSSets.begin(), it);
    }

    void CalculationLayer::CreatePipelines() {
        namespace Shaders = SCESS::Res::ShaderSources;

        //State shared by every pipeline: depth tested, no face culling (the isosurface winding is not consistent),
        //and alpha blending so the opacity sliders work. The old OpenGL renderer asked for depth testing off, but
        //Dahlia's text renderer re-enabled it every frame, so it always rendered with depth testing on.
        GraphicsPipelineDesc base;
        base.topology = PrimitiveTopology::Triangles;
        base.raster.cullMode = CullMode::None;
        base.depth = {true, true, CompareOp::Less};
        base.blend = BlendPreset::Alpha;
        base.targetLayout = mTargetLayout;

        //All three shaders only read positions from the mesh (location 0)
        const VertexBinding positionBinding{GPUMesh::POSITION_LOCATION, 3 * sizeof(float), false};
        const VertexAttribute positionAttribute{GPUMesh::POSITION_LOCATION, GPUMesh::POSITION_LOCATION,
                                                Format::RGB32F, 0};

        //Atoms: an instanced sphere, one instance per atom
        {
            GraphicsPipelineDesc desc = base;
            desc.vertexShader = ShaderModule::Create(ShaderStage::Vertex, Shaders::GetMoleculeDiffuseVert());
            desc.fragmentShader = ShaderModule::Create(ShaderStage::Fragment, Shaders::GetMoleculeDiffuseFrag());
            desc.vertexBindings = {
                positionBinding,
                {MOLECULE_INSTANCE_BINDING, MOLECULE_INSTANCE_FLOATS * sizeof(float), true}
            };
            desc.vertexAttributes = {
                positionAttribute,
                // Atom Scale
                {1, MOLECULE_INSTANCE_BINDING, Format::R32F, 0},
                // Atom Translation
                {2, MOLECULE_INSTANCE_BINDING, Format::RGB32F, 1 * sizeof(float)},
                // Atom Colour
                {3, MOLECULE_INSTANCE_BINDING, Format::RGB32F, 4 * sizeof(float)}
            };
            desc.pushConstantSize = sizeof(MoleculePushConstants);
            mMoleculePipeline = GraphicsPipeline::Create(desc);
        }

        //Bonds: a cylinder drawn once per bond
        {
            GraphicsPipelineDesc desc = base;
            desc.vertexShader = ShaderModule::Create(ShaderStage::Vertex, Shaders::GetSolidColourVert());
            desc.fragmentShader = ShaderModule::Create(ShaderStage::Fragment, Shaders::GetSolidColourFrag());
            desc.vertexBindings = {positionBinding};
            desc.vertexAttributes = {positionAttribute};
            desc.pushConstantSize = sizeof(SolidColourPushConstants);
            mBondPipeline = GraphicsPipeline::Create(desc);
        }

        //Isosurfaces: unindexed triangle lists, solid or wireframe
        {
            GraphicsPipelineDesc desc = base;
            desc.vertexShader = ShaderModule::Create(ShaderStage::Vertex, Shaders::GetIsosurfaceDiffuseVert());
            desc.fragmentShader = ShaderModule::Create(ShaderStage::Fragment, Shaders::GetIsosurfaceDiffuseFrag());
            desc.vertexBindings = {positionBinding};
            desc.vertexAttributes = {positionAttribute};
            desc.pushConstantSize = sizeof(IsosurfacePushConstants);
            mIsosurfacePipeline = GraphicsPipeline::Create(desc);

            desc.raster.polygonMode = PolygonMode::Lines;
            mIsosurfaceWirePipeline = GraphicsPipeline::Create(desc);
        }
    }

    void CalculationLayer::OnDestroy() {
        ImPlot::DestroyContext();

        mMoleculePipeline.reset();
        mBondPipeline.reset();
        mIsosurfacePipeline.reset();
        mIsosurfaceWirePipeline.reset();
        mAtomMesh.reset();
        mBondMesh.reset();
        mMoleculeBuffer.buffer.reset();
        mPosTrisurfBuffer.buffer.reset();
        mNegTrisurfBuffer.buffer.reset();
    }

    void CalculationLayer::Render(CommandList &cmd, const Camera3D &cam, const glm::mat4 &projectionMatrix) {
        mMoleculeBuffer.Flush();
        mPosTrisurfBuffer.Flush();
        mNegTrisurfBuffer.Flush();

        const glm::mat4 pvMatrix = projectionMatrix * cam.GetViewMatrix();
        const glm::vec4 lightPos(cam.GetActualPosition(), 0.0f);

        const bool drawMolecule = mMolecule != nullptr && mRenderMolecule;

        if (drawMolecule && mMoleculeBuffer.elementCount != 0) {
            cmd.BindPipeline(*mMoleculePipeline);
            cmd.PushConstants(MoleculePushConstants{pvMatrix, lightPos, {mMoleculeOpacity, 0.0f, 0.0f, 0.0f}});
            cmd.BindVertexBuffer(MOLECULE_INSTANCE_BINDING, *mMoleculeBuffer.buffer);
            mAtomMesh->Draw(cmd, mMoleculeBuffer.elementCount);
        }

        if (mMolRenderStyle == MoleculeStyle::BallStick && drawMolecule) {
            cmd.BindPipeline(*mBondPipeline);
            const glm::vec4 bondColour(glm::vec3(0.1f), mMoleculeOpacity);

            for (const BondPair &bond: mMolecule->GetBonds()) {
                const Atom &a1 = mMolecule->GetAtom(bond.atom1);
                const Atom &a2 = mMolecule->GetAtom(bond.atom2);
                glm::vec3 pos1 = glm::vec3(a1.GetPos().x(), a1.GetPos().y(), a1.GetPos().z());
                glm::vec3 pos2 = glm::vec3(a2.GetPos().x(), a2.GetPos().y(), a2.GetPos().z());
                glm::vec3 dirVec = pos1 - pos2;
                glm::vec3 midpoint = dirVec * 0.5f + pos2;

                glm::mat4 transformMat =
                        glm::translate(midpoint) *
                        glm::mat4(CalculateBondRotMatrix(glm::normalize(dirVec), {0.0f, 0.0f, 1.0f}))
                        * glm::scale(glm::vec3(1.0f, 1.0f, glm::length(dirVec)));

                cmd.PushConstants(SolidColourPushConstants{pvMatrix * transformMat, bondColour});
                mBondMesh->Draw(cmd);
            }
        }

        const auto drawIsosurface = [&](const DynamicVertexBuffer &surface, const glm::vec3 &colour) {
            if (surface.elementCount == 0)
                return;
            cmd.PushConstants(IsosurfacePushConstants{pvMatrix, lightPos, glm::vec4(colour, mIsoOpacity)});
            cmd.BindVertexBuffer(GPUMesh::POSITION_LOCATION, *surface.buffer);
            cmd.Draw(surface.elementCount);
        };

        cmd.BindPipeline(mPolygonMode == PolygonMode::Lines ? *mIsosurfaceWirePipeline : *mIsosurfacePipeline);
        drawIsosurface(mPosTrisurfBuffer, {1.0f, 0.0f, 0.0f});
        drawIsosurface(mNegTrisurfBuffer, {0.0f, 0.0f, 1.0f});
    }

    void CalculationLayer::RenderMOIsosurface() {
        const Eigen::MatrixXd *coeffsM = nullptr;

        mPosTrisurf.Stop();
        mNegTrisurf.Stop();

        switch (mCalcSCFMethod) {
            default:
                return;
            case SCFMethod::RHF:
                coeffsM = &mRHFWfn->C_Alpha;
                break;
            case SCFMethod::UHF:
                coeffsM = mSelectedSpin == 1 ? &mUHFWfn->C_Alpha : &mUHFWfn->C_Beta;
                break;
            case SCFMethod::KS_DFT:
                coeffsM = &mKSWfn->C;
                break;
        }

        std::shared_ptr<ScalarField> field = std::make_shared<ScalarField>();
        field->Initilise(mIsosurfaceResolution,
                         GenerateMolecularWavefunction<float>(coeffsM->col(mSelectedMO), mMoleculeBasisSet),
                         mIsosurfaceVolume);
        field->Evaluate(true);

        mPosTrisurf.Initilise(field, mIsosurfaceVolume, mIsosurfaceVal);
        mNegTrisurf.Initilise(field, mIsosurfaceVolume, -mIsosurfaceVal);

        mPosTrisurf.Evaluate(true);
        mNegTrisurf.Evaluate(true);
    }

    constexpr const char *PROG_SCREEN_NAME = "Calculation Progress";

    void CalculationLayer::SelectMolecule() {
        Molecule (*moleculeFuncs[])() =
        {
            &Molecule::GenerateWaterMolecule,
            &Molecule::GenerateAmmoniaMolecule,
            &Molecule::GenerateCOMolecule,
            &Molecule::GenerateEtheneMolecule,
            &Molecule::GenerateMethaneMolecule,
            &Molecule::GenerateLithiumHydrideMolecule,
            &Molecule::GenerateHydrogenMolecule
        };

        mMolecule = std::make_unique<Molecule>(moleculeFuncs[mSelectedMolecule]());
        UploadMolecule();
    }

    void CalculationLayer::UploadMolecule() {
        if (mMolecule == nullptr)
            return;

        mInvalidateLastSCF = true;

        std::vector<float> data(mMolecule->GetNoAtoms() * MOLECULE_INSTANCE_FLOATS);
        for (size_t i = 0; i < mMolecule->GetNoAtoms(); i++) {
            const Atom &a = mMolecule->GetAtom(i);
            const AtomStyle &style = AtomStyle::GetStyle(a.GetAtomSymbol());

            // Atom scale
            data[i * 7 + 0] = mMolRenderStyle == MoleculeStyle::SpaceFill ? style.mRadius : 0.6f;
            // Atom Translation
            data[i * 7 + 1] = (float) a.GetPos().x();
            data[i * 7 + 2] = (float) a.GetPos().y();
            data[i * 7 + 3] = (float) a.GetPos().z();
            // Atom colour
            data[i * 7 + 4] = style.mColour.x;
            data[i * 7 + 5] = style.mColour.y;
            data[i * 7 + 6] = style.mColour.z;
        }

        mMoleculeBuffer.Queue(data.data(), data.size());
    }

    bool CalculationLayer::SelectMoleculeFromFile() {
        nfdchar_t *outPath = nullptr;
        nfdresult_t result = NFD_OpenDialog("xyz", nullptr, &outPath);

        if (result != NFD_OKAY) return false;

        try {
            mMolecule = std::make_unique<Molecule>(Molecule::LoadFromXYZ((const char *) outPath));
        } catch (const std::exception &e) {
            free(outPath);
            mMolecule = nullptr;
            return false;
        }
        free(outPath);

        mMolecule->ConstructBonding(2.0);
        UploadMolecule();
        return true;
    }

    void CalculationLayer::ExecuteCalculation() {
        if (mMolecule == nullptr || mBasisSets.empty())
            return;

        mPosTrisurf.Reset();
        mNegTrisurf.Reset();

        mOrbitalsOrdered = false;
        mSelectedMO = -1;

        mCalcBasisSet = mBasisSets[mSelectedBasisSet].name;

        //Correlated methods are all built on an RHF reference
        if (mCorrelationMethod != CorrelationMethod::None)
            mSCFMethod = SCFMethod::RHF;

        mCalcCorrelationMethod = mCorrelationMethod;

        //F12 methods build their own OBS + CABS integrals and reference, so they produce no SCF result (orbitals)
        if (IsF12Method(mCalcCorrelationMethod)) {
            ExecuteF12Calculation();
            return;
        }

        mMoleculeBasisSet = BasisSet(mCalcBasisSet, *mMolecule);
        mCalcSCFMethod = mSCFMethod;

        //Decide on the UI thread whether the RHF reference can be reused; it must exist and be for this
        //molecule and basis set
        const bool runRHF = !mUseLastSCF || mInvalidateLastSCF || !mRHFWfn.has_value();

        mConsoleOutput.clear();
        mConsoleBuffer.Clear(); {
            std::lock_guard<std::mutex> lock(mQMCGraphMutex);
            mQMCGraphData.clear();
        }

        auto onQMCReportCycle = [this](size_t reportCycle, double projectedEnergy) {
            std::lock_guard<std::mutex> lock(mQMCGraphMutex);
            mQMCGraphData.push_back(glm::vec2{(float) reportCycle, (float) projectedEnergy});
        };

        mCalculationThread.Start([this, runRHF, onQMCReportCycle]() {
            constexpr double SCF_CONVERGENCE = 1E-10;

            //Exceptions must not escape the worker thread (e.g. an odd number of electrons for RHF)
            bool failed = false;
            try {
                if (mCalcCorrelationMethod != CorrelationMethod::None) {
                    //Run SCF as RHF if a SCF determinant is not available
                    if (runRHF)
                        mRHFWfn = HF::RHF(
                            *mMolecule,
                            mMoleculeBasisSet,
                            SCF_CONVERGENCE,
                            100,
                            &mints,
                            &mCalculationThread,
                            &mConsoleOutput);

                    if (!mCalculationThread.ShouldStop())
                        switch (mCalcCorrelationMethod) {
                            case CorrelationMethod::MP2:
                                mMP2Energy = HF::MP2_RHF_ClosedShellAO(*mRHFWfn,
                                                                       mints.eri,
                                                                       &mCalculationThread,
                                                                       &mConsoleOutput);
                                break;
                            case CorrelationMethod::CCSD: {
                                mCCSDWfn = CC::RCCSD_DIIS(
                                    *mRHFWfn, mints, &mCalculationThread, &mConsoleOutput);
                                break;
                            }

                            case CorrelationMethod::CCSD_T: {
                                mCCSD_TWfn = SCESS::CC::RCCSD_T(
                                    *mRHFWfn, mints, &mCalculationThread, &mConsoleOutput);
                                break;
                            }

                            case CorrelationMethod::FCI: {
                                mFCIWfn = CI::CI(*mRHFWfn,
                                                 mints,
                                                 mMolecule->CalculateVnn(),
                                                 -1,
                                                 &mCalculationThread,
                                                 &mConsoleOutput);
                                break;
                            }

                            case CorrelationMethod::FCIQMC: {
                                mQMCEnergy = QMC::FCI::FCIQMC(*mRHFWfn,
                                                                 mints,
                                                                 mMolecule->CalculateVnn(),
                                                                 &mCalculationThread,
                                                                 &mConsoleOutput,
                                                                 onQMCReportCycle);
                                break;
                            }

                            case CorrelationMethod::CCQMC: {
                                mQMCEnergy = QMC::CC::CCQMC(*mRHFWfn,
                                                               mints,
                                                               &mCalculationThread,
                                                               onQMCReportCycle);
                                break;
                            }

                            default:
                                break;
                        }
                } else {
                    switch (mCalcSCFMethod) {
                        default:
                            break;
                        case SCFMethod::RHF:
                            mRHFWfn = HF::RHF(
                                *mMolecule,
                                mMoleculeBasisSet,
                                SCF_CONVERGENCE,
                                100,
                                &mints,
                                &mCalculationThread,
                                &mConsoleOutput);
                            break;
                        case SCFMethod::UHF:
                            mUHFWfn = HF::UHF(*mMolecule,
                                              mMoleculeBasisSet,
                                              SCF_CONVERGENCE,
                                              100,
                                              &mints,
                                              &mCalculationThread,
                                              &mConsoleOutput);
                            break;
                        case SCFMethod::KS_DFT:
                            mKSWfn = DFT::KS_DFT(*mMolecule,
                                                 mMoleculeBasisSet,
                                                 DFTGridResToRes(mDFTGridRes),
                                                 &mCalculationThread,
                                                 &mConsoleOutput);
                            break;
                    }
                }
            } catch (const std::exception &e) {
                failed = true;
                mConsoleOutput << "Calculation failed: " << e.what() << '\n';
            }

            const bool cancelled = mCalculationThread.ShouldStop();
            if (cancelled)
                mConsoleOutput << "Calculation cancelled\n";

            std::cout << mConsoleBuffer.Snapshot();
            //The RHF reference is valid for the current molecule and basis set only if it was computed (and not
            //cancelled) in this run
            const bool ranRHF = mCalcCorrelationMethod != CorrelationMethod::None
                                    ? runRHF
                                    : mCalcSCFMethod == SCFMethod::RHF;
            if (ranRHF)
                mInvalidateLastSCF = cancelled || failed;
            mCalculationThread.MarkAsFinished();
        });
    }
    //
    void CalculationLayer::ExecuteF12Calculation() {
        mCalcSCFMethod = SCFMethod::None;
        const std::string cabsName = mCABSSets[mSelectedCABS].name;

        mConsoleOutput.clear();
        mConsoleBuffer.Clear();
        mF12Result.reset();

        mCalculationThread.Start([this, cabsName]() {
            try {
                const F12::F12Calculation calc = F12::PrepareF12Calculation(mCalcBasisSet, cabsName, *mMolecule,
                                                                            &mCalculationThread, &mConsoleOutput);
                F12Result result;
                result.method = mCalcCorrelationMethod;
                result.cabs = cabsName;
                result.E_HF = calc.hfWfn.E_HF;
                result.E_CABSSingles = calc.E_CABSSingles;

                bool succeeded = !mCalculationThread.ShouldStop();
                if (succeeded)
                    switch (mCalcCorrelationMethod) {
                        case CorrelationMethod::MP2_F12_Fixed:
                        case CorrelationMethod::MP2_F12_Opt: {
                            const F12::MP2F12Result mp2 =
                                    mCalcCorrelationMethod == CorrelationMethod::MP2_F12_Fixed
                                        ? F12::MP2_F12_Fixed(calc.hfWfn, calc.ints, calc.idxs)
                                        : F12::MP2_F12_Optimise(calc.hfWfn, calc.ints, calc.idxs);
                            result.E_MP2 = mp2.E_MP2;
                            result.E_F12 = mp2.E_F12_MP2;
                            mConsoleOutput << std::setprecision(12)
                                    << "E(MP2) = " << mp2.E_MP2 << " hartrees\n"
                                    << "E(F12 MP2) = " << mp2.E_F12_MP2 << " hartrees\n"
                                    << "E(Total) = " << result.Total() << " hartrees\n";
                            break;
                        }
                        case CorrelationMethod::CCSD_F12b: {
                            const F12::CCSD_F12bResult ccsd = F12::CCSD_F12b(calc, &mCalculationThread,
                                                                             &mConsoleOutput);
                            succeeded = ccsd.converged;
                            result.E_MP2 = ccsd.E_MP2;
                            result.E_CCSD = ccsd.E_CCSD;
                            result.E_Coupling = ccsd.E_Coupling;
                            result.E_F12 = ccsd.E_F12;
                            break;
                        }
                        default:
                            succeeded = false;
                            break;
                    }

                if (succeeded && !mCalculationThread.ShouldStop())
                    mF12Result = result;
            } catch (const std::exception &e) {
                mConsoleOutput << "Calculation failed: " << e.what() << '\n';
            }

            if (mCalculationThread.ShouldStop())
                mConsoleOutput << "Calculation cancelled\n";

            std::cout << mConsoleBuffer.Snapshot();
            mCalculationThread.MarkAsFinished();
        });
    }
    //
    void CalculationLayer::StopCalculation() {
        //Do not join here; the worker notices ShouldStop(), unwinds and marks itself finished, at which point the
        //UI re-enables and a new calculation may be started
        mCalculationThread.RequestStop();
    }

    // void CalculationLayer::DrawProgressScreen() {
    //     float dpiFactor = ImGui::GetWindowDpiScale();
    //     ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Once, ImVec2(0.5f, 0.5f));
    //     ImGui::SetNextWindowSize({600.0f * dpiFactor, 300.0f * dpiFactor}, ImGuiCond_Once);
    //     if (ImGui::BeginPopupModal(PROG_SCREEN_NAME)) {
    //         ImSpinner::SpinnerAtom("##atomSpinner", 20.0f * dpiFactor, 2.5f);
    //         ImGui::SameLine();
    //         ImGui::Text("Running %s single-point energy. Please wait...", [this]() {
    //             if (mCorrelationMethod != CorrelationMethod::None)
    //                 switch (mCorrelationMethod) {
    //                     case CorrelationMethod::MP2:
    //                         return "MP2";
    //                     case CorrelationMethod::CCSD:
    //                         return "CCSD";
    //                     case CorrelationMethod::CCSD_T:
    //                         return "CCSD(T)";
    //                     case CorrelationMethod::FCI:
    //                         return "FCI";
    //                     default:
    //                         return "";
    //                 };
    //
    //             switch (mSCFMethod) {
    //                 case SCFMethod::RHF:
    //                     return "RHF";
    //                 case SCFMethod::UHF:
    //                     return "UHF";
    //                 case SCFMethod::KS_DFT:
    //                     return "DFT";
    //                 default:
    //                     return "";
    //             };
    //         }());
    //
    //         const std::string &consoleStr = mConsoleOutput.str();
    //         ImGui::InputTextMultiline("##consoleOutput", (char *) consoleStr.c_str(), consoleStr.length(),
    //                                   {-FLT_MIN, -FLT_MIN}, ImGuiInputTextFlags_ReadOnly);
    //
    //         if (mCalculationThread.HasFinished())
    //             ImGui::CloseCurrentPopup();
    //
    //         ImGui::EndPopup();
    //     }
    // }

    void CalculationLayer::OrderOrbitals() {
        auto sortCoeffs = [this](Eigen::VectorXd &Ep_, Eigen::MatrixXd &Coeffs,
                                 const Eigen::VectorXi &Occupation) {
            //Sort orbital indices by energy; degenerate orbitals keep their relative order
            std::vector<Eigen::Index> order(Ep_.rows());
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(),
                             [&Ep_](Eigen::Index a, Eigen::Index b) { return Ep_[a] < Ep_[b]; });

            Eigen::VectorXd sortedEp(Ep_.rows());
            Eigen::MatrixXd sortedCoeffs(Coeffs.rows(), Coeffs.cols());

            mMOLabels.reserve(mMOLabels.size() + Ep_.rows());
            for (Eigen::Index i = 0; i < Ep_.rows(); i++) {
                sortedEp[i] = Ep_[order[i]];
                sortedCoeffs.col(i) = Coeffs.col(order[i]);

                std::stringstream ss;
                ss.precision(3);
                ss << "MO " << i + 1 << "; E=" << sortedEp[i] << " hartrees. ";

                if (i + 1 < Occupation.rows())
                    if (Occupation(i + 1) == 0 && Occupation(i) > 0)
                        ss << "HOMO";
                if (i != 0)
                    if (Occupation(i - 1) > 0 && Occupation(i) == 0)
                        ss << "LUMO";

                mMOLabels.push_back(ss.str());
            }
            Ep_ = sortedEp;
            Coeffs = sortedCoeffs;
        };

        if (!HaveSCFResult())
            return;

        mMOLabels.clear();

        switch (mCalcSCFMethod) {
            default:
                return;
            case SCFMethod::RHF:
                sortCoeffs(mRHFWfn->Ep_Alpha, mRHFWfn->C_Alpha, mRHFWfn->Occ_Alpha);
                break;
            case SCFMethod::UHF:
                sortCoeffs(mUHFWfn->Ep_Alpha, mUHFWfn->C_Alpha, mUHFWfn->Occ_Alpha);
                sortCoeffs(mUHFWfn->Ep_Beta, mUHFWfn->C_Beta, mUHFWfn->Occ_Beta);
                break;
            case SCFMethod::KS_DFT:
                sortCoeffs(mKSWfn->Ep, mKSWfn->C, mKSWfn->Occ);
                break;
        }
    }

    void CalculationLayer::DrawQuantumCalcPanel() {
        if (ImGui::Begin("Quantum Calculation")) {
            if (mMolecule != nullptr) {
                ImGui::BeginDisabled(mCalculationThread.IsRunning());

                ImGui::SeparatorText("SCF Methods");

                ImGui::Indent(100.0);
                ImGui::SeparatorText("Hartree-Fock");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("Restricted Hartree Fock", reinterpret_cast<int *>(&mSCFMethod), 0)) {
                    mSCFMethod = SCFMethod::RHF;
                    mInvalidateLastSCF = true;
                }

                if (ImGui::RadioButton("Unrestricted Hartree Fock", reinterpret_cast<int *>(&mSCFMethod), 1)) {
                    mCorrelationMethod = CorrelationMethod::None;
                    mSCFMethod = SCFMethod::UHF;
                    mInvalidateLastSCF = true;
                }
                ImGui::Unindent(100.0);

                ImGui::SeparatorText("Kohn-Sham Density Functional Theory");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("Slater-X + VWM Corr", reinterpret_cast<int *>(&mSCFMethod), 2)) {
                    mCorrelationMethod = CorrelationMethod::None;
                    mSCFMethod = SCFMethod::KS_DFT;
                    mInvalidateLastSCF = true;
                }
                ImGui::Unindent(200.0);

                ImGui::SeparatorText("Electron Correlation Methods");
                ImGui::Indent(200.0);

                if (ImGui::RadioButton("None", reinterpret_cast<int *>(&mCorrelationMethod), -1))
                    mCorrelationMethod = CorrelationMethod::None;

                ImGui::Unindent(100.0);
                ImGui::SeparatorText("Moller-Plesset");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("MP2", reinterpret_cast<int *>(&mCorrelationMethod), 0)) {
                    mCorrelationMethod = CorrelationMethod::MP2;
                    mSCFMethod = SCFMethod::RHF;
                }
                ImGui::Unindent(100.0);

                ImGui::SeparatorText("Coupled Cluster");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("CCSD", reinterpret_cast<int *>(&mCorrelationMethod), 1)) {
                    mCorrelationMethod = CorrelationMethod::CCSD;
                    mSCFMethod = SCFMethod::RHF;
                }
                if (ImGui::RadioButton("CCSD(T)", reinterpret_cast<int *>(&mCorrelationMethod), 2)) {
                    mCorrelationMethod = CorrelationMethod::CCSD_T;
                    mSCFMethod = SCFMethod::RHF;
                }
                ImGui::Unindent(100.0);

                ImGui::SeparatorText("Configuration Interaction");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("Full CI", reinterpret_cast<int *>(&mCorrelationMethod), 3)) {
                    mCorrelationMethod = CorrelationMethod::FCI;
                    mSCFMethod = SCFMethod::RHF;
                }
                ImGui::Unindent(100.0);

                ImGui::SeparatorText("Quantum Monte Carlo");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("FCI", reinterpret_cast<int *>(&mCorrelationMethod), 4)) {
                    mCorrelationMethod = CorrelationMethod::FCIQMC;
                    mSCFMethod = SCFMethod::RHF;
                    std::lock_guard<std::mutex> lock(mQMCGraphMutex);
                    mQMCGraphData.clear();
                }
                if (ImGui::RadioButton("CC", reinterpret_cast<int *>(&mCorrelationMethod), 5)) {
                    mCorrelationMethod = CorrelationMethod::CCQMC;
                    mSCFMethod = SCFMethod::RHF;
                    std::lock_guard<std::mutex> lock(mQMCGraphMutex);
                    mQMCGraphData.clear();
                }
                ImGui::Unindent(100.0);

                ImGui::SeparatorText("Explicitly Correlated (F12)");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("MP2-F12 (fixed amplitudes)", reinterpret_cast<int *>(&mCorrelationMethod), 6)) {
                    mCorrelationMethod = CorrelationMethod::MP2_F12_Fixed;
                    mSCFMethod = SCFMethod::RHF;
                }
                if (ImGui::RadioButton("MP2-F12 (optimised amplitudes)", reinterpret_cast<int *>(&mCorrelationMethod), 7)) {
                    mCorrelationMethod = CorrelationMethod::MP2_F12_Opt;
                    mSCFMethod = SCFMethod::RHF;
                }
                if (ImGui::RadioButton("CCSD-F12b", reinterpret_cast<int *>(&mCorrelationMethod), 8)) {
                    mCorrelationMethod = CorrelationMethod::CCSD_F12b;
                    mSCFMethod = SCFMethod::RHF;
                }
                ImGui::Unindent(200.0);

                DrawBasisSetSelection();
                if (IsF12Method(mCorrelationMethod))
                    DrawCABSSelection();

                if (mSCFMethod == SCFMethod::KS_DFT) {
                    ImGui::SeparatorText("DFT Molecular Grid");
                    ImGui::Indent(100.0);
                    const char *molGridLabels[] = {"Low", "Medium", "High", "Ultrafine"};
                    ImGui::Combo("##molGridRes", reinterpret_cast<int *>(&mDFTGridRes), molGridLabels,
                                 IM_ARRAYSIZE(molGridLabels));
                    ImGui::Unindent(100.0);
                }

                if (!mCalculationThread.IsRunning() && mRHFWfn.has_value() && !IsF12Method(mCorrelationMethod)) {
                    ImGui::Checkbox("Use Last SCF Result", &mUseLastSCF);
                    ImGui::SameLine();
                }
                if (mInvalidateLastSCF && !IsF12Method(mCorrelationMethod)) {
                    ImGui::TextDisabled("SCF will be re-run");
                }

                ImGui::EndDisabled();

                ImGui::NewLine();

                if (!mCalculationThread.IsRunning()) {
                    ImGui::PushStyleColor(ImGuiCol_Button, {0.0f, 0.4f, 0.0f, 1.0f});
                    if (ImGui::Button("Execute Calculation", ImVec2{-FLT_MIN, 0.0f})) {
                        ExecuteCalculation();
                    }
                } else {
                    ImGui::PushStyleColor(ImGuiCol_Button, {0.8f, 0.1f, 0.0f, 1.0f});

                    const bool stopping = mCalculationThread.IsStopping();
                    ImGui::Text(stopping ? "Stopping %s single point energy..." : "Running %s single point energy...",
                                mCalcCorrelationMethod == CorrelationMethod::None
                                    ? SCFMethosToStr(mCalcSCFMethod)
                                    : CorrelationMethodToStr(mCalcCorrelationMethod));
                    ImGui::SameLine();
                    ImSpinner::SpinnerAtom("##atomSpinner", 12.0f * ImGui::GetWindowDpiScale(), 1.5f);

                    ImGui::BeginDisabled(stopping);
                    if (ImGui::Button(stopping ? "Stopping Calculation..." : "Stop Calculation",
                                      ImVec2{-FLT_MIN, 0.0f})) {
                        StopCalculation();
                    }
                    ImGui::EndDisabled();
                }

                ImGui::PopStyleColor();
            }
        }
        ImGui::End();
    }

    void CalculationLayer::DrawBasisSetSelection() {
        ImGui::SeparatorText("Basis Sets");
        if (mBasisSets.empty()) {
            ImGui::TextDisabled("No basis sets were found in Assets/Basis");
            return;
        }
        ImGui::Indent(100.0);
        //One collapsible group per family; the group holding the selected basis set starts open
        for (const std::string &family: BasisLibrary::ListFamilies()) {
            if (family == BasisLibrary::CABS_FAMILY)
                continue;

            ImGui::SetNextItemOpen(mBasisSets[mSelectedBasisSet].family == family, ImGuiCond_Appearing);
            if (!ImGui::TreeNode(family.c_str()))
                continue;
            for (int i = 0; i < (int) mBasisSets.size(); i++)
                if (mBasisSets[i].family == family &&
                    ImGui::RadioButton(mBasisSets[i].name.c_str(), &mSelectedBasisSet, i)) {
                    mInvalidateLastSCF = true;
                    SelectDefaultCABS();
                }
            ImGui::TreePop();
        }
        ImGui::Unindent(100.0);
    }

    void CalculationLayer::DrawCABSSelection() {
        ImGui::SeparatorText("Complementary Auxiliary Basis Set (CABS)");
        if (mCABSSets.empty() || mBasisSets.empty())
            return;
        ImGui::Indent(100.0);
        if (ImGui::BeginCombo("##cabsSelection", mCABSSets[mSelectedCABS].name.c_str())) {
            std::string lastFamily;
            for (int i = 0; i < (int) mCABSSets.size(); i++) {
                const BasisSetInfo &info = mCABSSets[i];
                if (info.family != lastFamily) {
                    ImGui::SeparatorText(info.family == BasisLibrary::CABS_FAMILY
                                             ? "CABS (OptRI)"
                                             : info.family.c_str());
                    lastFamily = info.family;
                }
                const bool selected = i == mSelectedCABS;
                if (ImGui::Selectable(info.name.c_str(), selected)) {
                    mSelectedCABS = i;
                    mCABSFollowsOBS = false;
                }
                if (selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        ImGui::SameLine();
        if (ImGui::Checkbox("Match basis set", &mCABSFollowsOBS))
            SelectDefaultCABS();
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("Use the OptRI set recommended for the orbital basis set, e.g.\n"
                              "cc-pVDZ-F12 -> cc-pVDZ-F12-OptRI, cc-pVDZ & aug-cc-pVDZ -> aug-cc-pVDZ-OptRI");
        if (mCABSFollowsOBS && BasisLibrary::DefaultCABS(mBasisSets[mSelectedBasisSet].name).empty())
            ImGui::TextDisabled("%s has no matching CABS; choose one above", mBasisSets[mSelectedBasisSet].name.c_str());
        ImGui::Unindent(100.0);
    }

    void CalculationLayer::DrawResultsPanel() {
        if (ImGui::Begin("Results")) {
            if (mCalculationThread.HasFinished()) {
                if (mCalcCorrelationMethod != CorrelationMethod::None && HaveCorrelationResult() &&
                    ImGui::CollapsingHeader(
                        "Electron Correlation Result", ImGuiTreeNodeFlags_DefaultOpen)) {
                    auto DrawMP2 = [this]() {
                        if (!mMP2Energy.has_value())return;
                        ImGui::SeparatorText("Moller-Plesset Order 2");
                        ImGui::Indent();
                        ImGui::Text("E(MP2) = %.8f hartrees", mMP2Energy.value());
                        ImGui::TextColored({0.0, 1.0, 0.0, 1.0},
                                           "E(Total) = E(Hartree-Fock) + E(MP2) = %.8f hartrees",
                                           mMP2Energy.value() + mRHFWfn->E_HF);
                        ImGui::Unindent();
                    };

                    auto DrawCCSD = [this]() {
                        if (!mCCSDWfn)return;
                        ImGui::SeparatorText("Coupled Cluster Singles & Doubles");
                        ImGui::Indent();
                        ImGui::Text("E(CCSD) = %.8f hartrees", mCCSDWfn->E_CCSD);
                        ImGui::TextColored({0.0, 1.0, 0.0, 1.0},
                                           "E(Total) = E(Hartree-Fock) + E(CCSD) = %.8f hartrees",
                                           mCCSDWfn->E_CCSD + mRHFWfn->E_HF);
                        ImGui::Unindent();
                    };

                    auto DrawCCSD_T = [this]() {
                        if (!mCCSD_TWfn)return;
                        ImGui::SeparatorText("Coupled Cluster Singles & Doubles with Pertubative Triples");
                        ImGui::Indent();
                        ImGui::Text("E(CCSD) = %.8f hartrees", mCCSD_TWfn->E_CCSD);
                        ImGui::Text("E(T) = %.8f hartrees", mCCSD_TWfn->E_T);
                        ImGui::TextColored({0.0, 1.0, 0.0, 1.0},
                                           "E(Total) = E(Hartree-Fock) + E(CCSD) + E(T) = %.8f hartrees",
                                           mCCSD_TWfn->E_CCSD + mRHFWfn->E_HF + mCCSD_TWfn->
                                           E_T);
                        ImGui::Unindent();
                    };

                    auto DrawFCI = [this]() {
                        if (!mFCIWfn)return;
                        ImGui::SeparatorText("Full Configuration Interaction");
                        ImGui::Indent();
                        ImGui::Text("E(FCI Correlation) = %.8f hartrees", mFCIWfn->E_CI(0) - mRHFWfn->E_HF);
                        ImGui::TextColored({0.0, 1.0, 0.0, 1.0},
                                           "E(Total) = E(Hartree Fock) + E(FCI Correlation) = %.8f hartrees",
                                           mFCIWfn->E_CI(0));
                        ImGui::Unindent();
                    };

                    auto DrawF12 = [this]() {
                        if (!mF12Result)return;
                        const F12Result &r = *mF12Result;
                        ImGui::SeparatorText(CorrelationMethodToStr(r.method));
                        ImGui::Indent();
                        ImGui::BulletText("OBS: %s, CABS: %s", mCalcBasisSet.c_str(), r.cabs.c_str());
                        ImGui::Text("E(HF) = %.8f hartrees", r.E_HF);
                        ImGui::Text("E(CABS singles correction) = %.8f hartrees", r.E_CABSSingles);
                        ImGui::Text("E(MP2) = %.8f hartrees", r.E_MP2);
                        if (r.method == CorrelationMethod::CCSD_F12b) {
                            ImGui::Text("E(CCSD, with geminal coupling) = %.8f hartrees", r.E_CCSD);
                            ImGui::Text("    of which E(coupling) = %.8f hartrees", r.E_Coupling);
                            ImGui::Text("E(F12 geminals) = %.8f hartrees", r.E_F12);
                        } else
                            ImGui::Text("E(F12 MP2) = %.8f hartrees", r.E_F12);
                        ImGui::Text("E(%s correlation) = %.8f hartrees", CorrelationMethodToStr(r.method),
                                    r.Correlation());
                        ImGui::TextColored({0.0, 1.0, 0.0, 1.0},
                                           "E(Total) = E(HF) + E(CABS singles) + E(correlation) = %.8f hartrees",
                                           r.Total());
                        ImGui::Unindent();
                    };

                    switch (mCalcCorrelationMethod) {
                        case CorrelationMethod::MP2_F12_Fixed:
                        case CorrelationMethod::MP2_F12_Opt:
                        case CorrelationMethod::CCSD_F12b: {
                            DrawF12();
                            break;
                        }

                        case CorrelationMethod::CCSD_T: {
                            DrawCCSD_T();
                            DrawMP2();
                            break;
                        }

                        case CorrelationMethod::CCSD: {
                            DrawCCSD();
                            DrawMP2();
                            break;
                        }

                        case CorrelationMethod::MP2: {
                            DrawMP2();
                            break;
                        }

                        case CorrelationMethod::FCI: {
                            DrawFCI();
                            break;
                        }

                        case CorrelationMethod::None:
                        default:
                            break;
                    }
                }

                if (HaveSCFResult() && ImGui::CollapsingHeader("SCF Result", ImGuiTreeNodeFlags_DefaultOpen)) {
                    if (mCalcSCFMethod != SCFMethod::None) {
                        if (!mOrbitalsOrdered) {
                            OrderOrbitals();
                            mOrbitalsOrdered = true;
                        }

                        ImGui::BulletText("Basis Set: %s", mCalcBasisSet.c_str());

                        if (mCalcSCFMethod == SCFMethod::KS_DFT)
                            ImGui::BulletText("Total Kohn-Sham DFT Electronic Energy : %.10f hartrees",
                                              mKSWfn->E_KS);
                        else
                            ImGui::BulletText("Total Hartree Fock Energy: %.10f hartrees",
                                              mCalcSCFMethod == SCFMethod::RHF ? mRHFWfn->E_HF : mUHFWfn->E_HF);


                        ImGui::Indent();
                        if (ImGui::CollapsingHeader("Orbitals", ImGuiTreeNodeFlags_DefaultOpen)) {
                            if (mSelectedMO != -1) {
                                if (mCalcSCFMethod == SCFMethod::RHF || mCalcSCFMethod == SCFMethod::KS_DFT) {
                                    Eigen::VectorXi &Occ = mCalcSCFMethod == SCFMethod::RHF
                                                               ? mRHFWfn->Occ_Alpha
                                                               : mKSWfn->Occ;
                                    Eigen::VectorXd &epsilon = mCalcSCFMethod == SCFMethod::RHF
                                                                   ? mRHFWfn->Ep_Alpha
                                                                   : mKSWfn->Ep;

                                    ImGui::TextColored({1.0f, 1.0f, 0.0f, 1.0f}, "Previewing MO %i",
                                                       mSelectedMO + 1);
                                    if (mSelectedMO + 1 < Occ.rows())
                                        if (Occ(mSelectedMO + 1) == 0 && Occ(
                                                mSelectedMO) > 0)
                                            ImGui::TextColored({1.0f, 0.0f, 0.0f, 1.0f}, "HOMO");
                                    if (mSelectedMO - 1 >= 0)
                                        if (Occ(mSelectedMO - 1) > 0 && Occ(
                                                mSelectedMO) == 0)
                                            ImGui::TextColored({0.0f, 0.5f, 0.5f, 1.0f}, "LUMO");
                                    ImGui::BulletText("Energy: %.3f hartrees", epsilon(mSelectedMO));
                                    ImGui::BulletText("Occupation: %i e-", Occ(mSelectedMO));
                                } else if (mCalcSCFMethod == SCFMethod::UHF) {
                                    const auto &Ep_ =
                                            mSelectedSpin == 1 ? mUHFWfn->Ep_Alpha : mUHFWfn->Ep_Beta;
                                    const auto &occupation = mSelectedSpin == 1
                                                                 ? mUHFWfn->Occ_Alpha
                                                                 : mUHFWfn->Occ_Beta;
                                    ImGui::TextColored({1.0f, 1.0f, 0.0f, 1.0f}, "Previewing MO %i",
                                                       mSelectedMO + 1);
                                    if (mSelectedMO + 1 < occupation.rows())
                                        if (occupation(mSelectedMO + 1) == 0 && occupation(mSelectedMO) > 0)
                                            ImGui::TextColored({1.0f, 0.0f, 0.0f, 1.0f}, "HOMO");
                                    if (mSelectedMO - 1 >= 0)
                                        if (occupation(mSelectedMO - 1) > 0 && occupation(mSelectedMO) == 0)
                                            ImGui::TextColored({0.0f, 0.5f, 0.5f, 1.0f}, "LUMO");
                                    ImGui::BulletText(
                                        "Spin: %s", mSelectedSpin == 1 ? "Alpha (up)" : "Beta (down)");
                                    ImGui::BulletText("Energy: %.3f hartrees", Ep_(mSelectedMO));
                                    ImGui::BulletText("Occupation: %i e-", occupation(mSelectedMO));
                                }
                                ImGui::Separator();
                            }

                            int selectedOption = -1;
                            int selectedOption2 = -1;
                            if (!(mPosTrisurf.IsRunning() || mNegTrisurf.IsRunning())) {
                                const char **orbitalLabels = (const char **) alloca(
                                    sizeof(const char **) * mMOLabels.size());

                                for (uint32_t i = 0; i < mMOLabels.size(); i++)
                                    orbitalLabels[i] = mMOLabels[i].c_str();

                                ImGui::Text("Select an MO to generate an isosurface");
                                if (mCalcSCFMethod == SCFMethod::RHF || mCalcSCFMethod == SCFMethod::KS_DFT)
                                    ImGui::ListBox("##listbox1", &selectedOption, &orbitalLabels[0],
                                                   (int) mMOLabels.size());
                                else if (mCalcSCFMethod == SCFMethod::UHF) {
                                    float width = ImGui::GetColumnWidth() * 0.5f;

                                    ImGui::Text("Alpha (spin up) Orbitals");
                                    ImGui::SameLine(width, 10.0f);
                                    ImGui::Text("Beta (spin down) Orbitals");
                                    ImGui::PushItemWidth(width);
                                    ImGui::ListBox("##listbox1", &selectedOption, &orbitalLabels[0],
                                                   (int) mUHFWfn->Ep_Alpha.rows());
                                    ImGui::SameLine();
                                    ImGui::ListBox("##listbox2", &selectedOption2,
                                                   &orbitalLabels[mUHFWfn->Ep_Alpha.rows()],
                                                   (int) mUHFWfn->Ep_Beta.rows());
                                    ImGui::PopItemWidth();
                                }
                            }

                            if (selectedOption != -1 || selectedOption2 != -1) {
                                if (mCalcSCFMethod == SCFMethod::RHF || mCalcSCFMethod == SCFMethod::KS_DFT)
                                    mSelectedMO = selectedOption;
                                else if (mCalcSCFMethod == SCFMethod::UHF) {
                                    if (selectedOption != -1) {
                                        mSelectedMO = selectedOption;
                                        mSelectedSpin = 1;
                                    } else {
                                        mSelectedMO = selectedOption2;
                                        mSelectedSpin = -1;
                                    }
                                }

                                mIsosurfaceUploaded = false;
                                mIsoWasNull = false;
                                RenderMOIsosurface();
                            }

                            if ((mPosTrisurf.HasFinished() && mNegTrisurf.HasFinished()) && !mIsosurfaceUploaded) {
                                // if (mPosTrisurf.triangles.size() != 0)
                                UploadTriSurface(mPosTrisurfBuffer, mPosTrisurf);
                                // if (mNegTrisurf.triangles.size() != 0)
                                UploadTriSurface(mNegTrisurfBuffer, mNegTrisurf);

                                mIsoWasNull = mPosTrisurf.GetTriangles().empty() && mNegTrisurf.GetTriangles().empty();

                                mPosTrisurf.Reset();
                                mNegTrisurf.Reset();
                                mIsosurfaceUploaded = true;
                            }

                            if (mPosTrisurf.IsRunning() || mNegTrisurf.IsRunning()) {
                                ImSpinner::SpinnerTwinBall("##spinner1", 10.0f * ImGui::GetWindowDpiScale(),
                                                           10.0f * ImGui::GetWindowDpiScale(), 3.0, 6.0,
                                                           {1.0f, 0.0f, 0.0f, 1.0f},
                                                           ImSpinner::white, 3.0f, 1);
                                ImGui::SameLine();
                                ImGui::Text("Generating isosurface...");
                                ImGui::SameLine();
                                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4{0.6f, 0.0f, 0.0f, 1.0f});
                                if (ImGui::Button("Cancel", {100, 0})) {
                                    mPosTrisurf.Reset();
                                    mNegTrisurf.Reset();
                                }
                                ImGui::PopStyleColor();
                            } else {
                                if (mIsosurfaceUploaded) {
                                    if (!mIsoWasNull)
                                        ImGui::TextColored({0.0f, 1.0f, 0.0f, 1.0f}, "Isosurface generated.");
                                    else
                                        ImGui::TextColored({1.0f, 0.0f, 0.0f, 1.0f},
                                                           "Wavefunction has no valid isosurface in the given parameters.");
                                }

                                if (ImGui::CollapsingHeader("Isosurface Options", ImGuiTreeNodeFlags_DefaultOpen)) {
                                    ImGui::Indent();

                                    float xPos = ImGui::GetCursorPosX() + 200.0f;

                                    ImGui::Text("Value: ");
                                    ImGui::SameLine();
                                    ImGui::SetCursorPosX(xPos);
                                    ImGui::DragFloat("##isosurfVal", &mIsosurfaceVal, 0.001f, 1E-10f, 1.0f);

                                    ImGui::Text("Resolution: ");
                                    ImGui::SameLine();
                                    ImGui::SetCursorPosX(xPos);
                                    ImGui::DragInt("##isoResolutionVal", &mIsosurfaceResolution, 1.0f, 10,
                                                   MAX_ISOSURFACE_RESOLUTION, "%d", ImGuiSliderFlags_AlwaysClamp);

                                    ImGui::Text("Volume: ");
                                    ImGui::SameLine();
                                    ImGui::SetCursorPosX(xPos);
                                    ImGui::DragFloat("bohr##isoVolumeVal", &mIsosurfaceVolume, 0.5f, 0.01f, 100.0f);

                                    ImGui::Text("Opacity: ");
                                    ImGui::SameLine();
                                    ImGui::SetCursorPosX(xPos);
                                    ImGui::DragFloat("##isoOpacity", &mIsoOpacity, 0.0f, 0.01f, 1.0f);

                                    ImGui::Unindent();
                                }
                            }
                        }
                    }
                }
            }
        }
        ImGui::End();
    }

    void CalculationLayer::DrawMolSelectPanel() {
        static const char *moleculeNames[] =
        {
            "Water",
            "Ammonia",
            "Carbon Monoxide",
            "Ethene",
            "Methane",
            "Lithium Hydride",
            "Hydrogen"
        };

        if (ImGui::Begin("Select Molecule")) {
            ImGui::BeginDisabled(mCalculationThread.IsRunning());

            ImGui::Indent(100.0f);
            int lastSelectedMol = mSelectedMolecule;
            for (int i = 0; i < IM_ARRAYSIZE(moleculeNames); i++)
                if (ImGui::RadioButton(moleculeNames[i], &mSelectedMolecule, i))
                    mInvalidateLastSCF = true;

            if (ImGui::RadioButton("Load from XYZ File...", &mSelectedMolecule, IM_ARRAYSIZE(moleculeNames) + 1)) {
                if (!SelectMoleculeFromFile())
                    mSelectedMolecule = lastSelectedMol;
            } else if (mSelectedMolecule != lastSelectedMol && mSelectedMolecule != -1)
                SelectMolecule();

            ImGui::Unindent(100.0f);
            ImGui::NewLine();
            ImGui::SeparatorText("Selected Molecule");
            ImGui::Indent(100.0f);

            ImGui::EndDisabled();


            ImGui::Text("Draw Molecule: ");
            float xPos = ImGui::GetCursorPosX() + ImGui::GetItemRectSize().x + 20.0f;
            ImGui::SameLine();
            ImGui::SetCursorPosX(xPos);
            ImGui::Checkbox("##drawMolCheck", &mRenderMolecule);
            if (mRenderMolecule) {
                ImGui::Text("Opacity: ");
                ImGui::SameLine();
                ImGui::SetCursorPosX(xPos);
                ImGui::SetCursorPosX(xPos);
                ImGui::SliderFloat("##molOpacity", &mMoleculeOpacity, 0.0f, 1.0f, "%.2f");

                ImGui::Text("Render Style: ");
                ImGui::SameLine();
                ImGui::SetCursorPosX(xPos);
                static const char *renderStyles[] = {"Spacefill", "Ball & Stick"};

                MoleculeStyle lastStyle = mMolRenderStyle;
                ImGui::Combo("##molRenderStyle", (int *) &mMolRenderStyle, renderStyles, IM_ARRAYSIZE(renderStyles));

                if (lastStyle != mMolRenderStyle)
                    UploadMolecule();
            }
        }
        ImGui::End();
    }

    void CalculationLayer::DrawConsolePanel() {
        if (ImGui::Begin("Console")) {
            mConsoleSnapshot = mConsoleBuffer.Snapshot();
            ImGui::InputTextMultiline("##consoleOutput", mConsoleSnapshot.data(), mConsoleSnapshot.length() + 1,
                                      {-FLT_MIN, -FLT_MIN}, ImGuiInputTextFlags_ReadOnly);
        }
        ImGui::End();
    }

    void CalculationLayer::DrawQMCGraph() {
        if (ImGui::Begin("Quantum Monte Carlo")
            &&
            (mCorrelationMethod == CorrelationMethod::CCQMC || mCorrelationMethod == CorrelationMethod::FCIQMC)) {
            ImGui::Text("%s", [this]() {
                switch (mCorrelationMethod) {
                    default:
                        return "";
                    case CorrelationMethod::CCQMC:
                        return "Coupled cluster quantum monte carlo";
                    case CorrelationMethod::FCIQMC:
                        return "Full CI quantum monte carlo";
                }
            }());
            if (mCalculationThread.HasFinished() && mRHFWfn.has_value()) {
                ImGui::Text("Final QMC electronic energy: %.8f", mQMCEnergy);
                ImGui::Text("Final QMC correlation energy: %.8f", mQMCEnergy - mRHFWfn->E_HF);
            }

            //Copy the points so the calculation thread is not blocked while plotting
            std::vector<glm::vec2> graphData; {
                std::lock_guard<std::mutex> lock(mQMCGraphMutex);
                graphData.assign(mQMCGraphData.begin(), mQMCGraphData.end());
            }

            if (!graphData.empty() && ImPlot::BeginPlot("Projected Energy", ImVec2(-1, -1))) {
                ImPlot::SetupAxes("Cycle No.", "E(Proj) / hartrees");
                ImPlot::SetupAxisLimits(ImAxis_X1, graphData.back().x, graphData.back().x + 250,
                                        ImGuiCond_Once);
                ImPlot::PushStyleColor(ImPlotCol_Line, ImVec4(1.0f, 0.0f, 0.0f, 1.0f));
                ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 1.0f);
                ImPlot::PlotLine<float>("E(Proj)",
                                        &graphData.data()->x,
                                        &graphData.data()->y,
                                        static_cast<int>(graphData.size()),
                                        0, 0, sizeof(glm::vec2));
                ImPlot::PopStyleColor();
                ImPlot::PopStyleVar();

                ImPlot::EndPlot();
            }
        }
        ImGui::End();
    }

    void CalculationLayer::OnImGuiRender() {
        DrawMolSelectPanel();

        //DrawProgressScreen();
        if (mCalculationThread.IsRunning())
            ImGui::OpenPopup(PROG_SCREEN_NAME);

        DrawQuantumCalcPanel();
        DrawResultsPanel();
        DrawConsolePanel();
        DrawQMCGraph();
    }

    void CalculationLayer::OnEvent(Event &e) {
    }
}
