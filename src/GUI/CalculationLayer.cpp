#include <pch.h>
#include "CalculationLayer.h"

#include <ImGui/imgui.h>
#include <nfd/nfd.h>

#include "GraphicsEngine/Scene/Camera3D.h"
#include "GraphicsEngine/3D/VertexArray.h"
#include "GraphicsEngine/3D/ShaderProgram.h"
#include "GraphicsEngine/Maths/Geometry.h"

#include "Chemistry/Quantum/Basis/BasisUtils.h"
#include "Chemistry/Quantum/Molecule.h"
#include "Chemistry/Quantum/Basis/Gaussian.h"
#include "Chemistry/Quantum/HF/HartreeFock.h"
#include "Chemistry/Quantum/DFT/DFT.h"
#include "Maths/Isosurfaces.h"

#include "Application/Application.h"
#include "Chemistry/Quantum/CC/CoupledCluster.h"

#include "GUI/ImSpinner.h"

namespace SCESS {
    using namespace GEng;

    void CalculationLayer::UploadTriSurface(VertexArray &vArray, const TriSurf &tri) {
#ifdef _DEBUG
        if (!tri.HasFinished()) {
            throw std::runtime_error("Trisurface generation thread has not completed.");
        }
#endif

        if (vArray.GetNoOfVBuffers() != 0) {
            const sRef<VertexBuffer> &vBuffer = vArray.GetVBuffer(0);
            vBuffer->StoreData((const float *) (tri.GetTriangles().data()),
                               tri.GetTriangles().size() * 3 * 3);
            vArray.UpdateElementLength(tri.GetTriangles().size() * 3);
        }
    }

    void CalculationLayer::OnInit() {
        mDiffuseProg = ShaderProgram::CreateFromFile(
            Application::Get().GetResourcePath(R"(Shaders/ColourDiffuse.glsl)").string());

        mMoleculeDiffuseProg = ShaderProgram::CreateFromFile(
            Application::Get().GetResourcePath(R"(Shaders/MoleculeColourDiffuse.glsl)").string());

        mSolidColourProg = ShaderProgram::CreateFromFile(
            Application::Get().GetResourcePath(R"(Shaders/SolidColour.glsl)").string());

        mPosTrisurfVAO = VertexArray::Create();
        mPosTrisurfVAO->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));
        mNegTrisurfVAO = VertexArray::Create();
        mNegTrisurfVAO->AddVertexBuffer(VertexBuffer::Create(BufferUsage::DynamicUpdate));

        mMoleculeVAO = VertexArray::Create();
        Mesh_F32 sphereMesh = Geometry::CreateQuadsphere(1.0f, 30);
        mMoleculeVAO->AddVertexBuffer(VertexBuffer::Create(sphereMesh.GetVertexArray().data(),
                                                           sphereMesh.GetVertexFloatCount()));
        mMoleculeVAO->AddIndexBuffer(IndexBuffer::Create(sphereMesh.GetIndicesArray().data(),
                                                         sphereMesh.GetIndexIntCount()));

        // All molecule data is packed in the following format:
        //(Atom Scale; float) (Atom Translation; vec3) (Atom Colour; vec3) --
        mMoleculeBuffer = VertexBuffer::Create(BufferUsage::DynamicUpdate);
        // Atom Scale
        mMoleculeVAO->AddVertexBuffer(mMoleculeBuffer, VertexSpecification{1, 7, 0, 1});
        // Atom Translation
        mMoleculeVAO->AddVertexBuffer(mMoleculeBuffer, VertexSpecification{3, 7, 1, 1});
        // Atom Colour
        mMoleculeVAO->AddVertexBuffer(mMoleculeBuffer, VertexSpecification{3, 7, 4, 1});

        mMoleculeBondVAO = VertexArray::Create(Geometry::CreateCylinder(0.2f, 0.7f, 20));

        mBasisSets = BasisUtils::ListBasisSets();
        mSelectedBasisSet = std::max(0,
                                     (int) std::distance(mBasisSets.begin(),
                                                         std::find(mBasisSets.begin(), mBasisSets.end(), "sto-3g")));
    }

    void CalculationLayer::OnDestroy() {
    }

    void CalculationLayer::OnUpdate(const sRef<SCESS::Framebuffer> &screenBuffer, const Camera3D &cam,
                                    const glm::mat4 &projectionMatrix, const glm::vec2 &nearFarPlanes) {
        GPUContext *ctx = GPUContext::Get();
        ctx->SetPolygonMode(GPUContext::PolygonMode::Solid);
        //ctx->SetFaceCulling(false);

        glm::mat4 pvMatrix = projectionMatrix * cam.GetViewMatrix();

        if (mMolecule != nullptr && mRenderMolecule) {
            mMoleculeDiffuseProg->Bind();
            mMoleculeDiffuseProg->UploadMatF("uMPVMatrix", pvMatrix);
            mMoleculeDiffuseProg->UploadVecF("uLightPos", cam.GetActualPosition());
            mMoleculeDiffuseProg->UploadF("uOpacity", mMoleculeOpacity);

            mMoleculeVAO->Bind();
            mMoleculeVAO->SubmitArrays(Primitive::Triangles, mMolecule->GetNoAtoms());
            mMoleculeVAO->Unbind();

            mMoleculeDiffuseProg->Unbind();
        }

        // if (mRenderMolecule)
        //{

        //	if (mMolecule != nullptr)
        //	{
        //		mSphereVAO->Bind();
        //		for (const Atom& a : mMolecule->GetAtoms())
        //		{
        //			const AtomStyle& style = AtomStyle::GetStyle(a.GetAtomSymbol());
        //			mDiffuseProg->UploadMatF("uMPVMatrix", pvMatrix * glm::translate(glm::vec3{ a.GetPos().x(), a.GetPos().y() ,a.GetPos().z() }) *
        //				glm::scale(glm::vec3{ style.mRadius }));
        //			mDiffuseProg->UploadVecF("uColour", style.mColour);

        //			mSphereVAO->SubmitArrays();
        //		}
        //		mSphereVAO->Unbind();
        //	}
        //}

        if (mMolRenderStyle == MoleculeStyle::BallStick && mMolecule != nullptr && mRenderMolecule) {
            mSolidColourProg->Bind();
            // mSolidColourProg->UploadVecF("uLightPos", cam.GetActualPosition());
            mSolidColourProg->UploadMatF("uPVMatrix", pvMatrix);
            mSolidColourProg->UploadVecF("uColour", glm::vec3(0.1f));
            mSolidColourProg->UploadF("uOpacity", mMoleculeOpacity);

            for (const BondPair &bond: mMolecule->GetBonds()) {
                const Atom &a1 = mMolecule->GetAtom(bond.atom1);
                const Atom &a2 = mMolecule->GetAtom(bond.atom2);
                glm::vec3 pos1 = glm::vec3(a1.GetPos().x(), a1.GetPos().y(), a1.GetPos().z());
                glm::vec3 pos2 = glm::vec3(a2.GetPos().x(), a2.GetPos().y(), a2.GetPos().z());
                glm::vec3 dirVec = pos1 - pos2;
                glm::vec3 midpoint = dirVec * 0.5f + pos2;

                glm::mat3 rotationMat = GenerateBondRotMatrix(glm::normalize(dirVec), {0.0f, 0.0f, 1.0f});

                glm::mat4 transformMat =
                        glm::translate(midpoint) *
                        glm::mat4(rotationMat) * glm::scale(glm::vec3(1.0f, 1.0f, dirVec.length()));

                mSolidColourProg->UploadMatF("uTransformMatrix", transformMat);

                mMoleculeBondVAO->Bind();
                mMoleculeBondVAO->SubmitArrays();
            }
        }

        ctx->SetPolygonMode(mPolygonMode);

        mDiffuseProg->Bind();
        mDiffuseProg->UploadVecF("uLightPos", cam.GetActualPosition());
        mDiffuseProg->UploadMatF("uPVMatrix", pvMatrix);
        mDiffuseProg->UploadMatF("uTransformMatrix", glm::mat4(1.0f));
        mDiffuseProg->UploadF("uOpacity", mIsoOpacity);

        mDiffuseProg->UploadVecF("uColour", {1.0f, 0.0f, 0.0f});
        mPosTrisurfVAO->Bind();
        mPosTrisurfVAO->SubmitArrays();
        mPosTrisurfVAO->Unbind();

        mDiffuseProg->UploadVecF("uColour", {0.0f, 0.0f, 1.0f});
        mNegTrisurfVAO->Bind();
        mNegTrisurfVAO->SubmitArrays();
        mNegTrisurfVAO->Unbind();

        mDiffuseProg->Unbind();
    }

    void CalculationLayer::RenderMOIsosurface() {
        const Eigen::MatrixXd *coeffsM = nullptr;

        mPosTrisurf.Stop();
        mNegTrisurf.Stop();

        if (mSCFMethod == SCFMethod::RHF) {
            coeffsM = &mRHFResult->Coeffs;
        } else {
            coeffsM = mSelectedSpin == 1 ? &mUHFResult->CoeffsAlpha : &mUHFResult->CoeffsBeta;
        }

        std::shared_ptr<ScalarField> field = std::make_shared<ScalarField>();
        field->Initilise(mIsosurfaceResolution,
                         GenerateMolecularWavefunction<float>(coeffsM->col(mSelectedMO), *mMolecule),
                         mIsosurfaceVolume);
        field->Evaluate(true);

        mPosTrisurf.Initilise(field, mIsosurfaceVolume, mIsosurfaceVal);
        mNegTrisurf.Initilise(field, mIsosurfaceVolume, -mIsosurfaceVal);

        mPosTrisurf.Evaluate(true);
        mNegTrisurf.Evaluate(true);
    }

    constexpr const char *PROG_SCREEN_NAME = "Calculation Progress";

    void CalculationLayer::SelectMolecule() {
        static SCESS::Molecule (*moleculeFuncs[])() =
        {
            {&GenerateWaterMolecule},
            {&GenerateAmmoniaMolecule},
            {&GenerateCOMolecule},
            {&GenerateEtheneMolecule},
            {&GenerateMethaneMolecule},
            {&GenerateLithiumHydrideMolecule},
            {&GenerateHydrogenMolecule}
        };

        mMolecule = std::make_unique<SCESS::Molecule>(moleculeFuncs[mSelectedMolecule]());
        UploadMolecule();
    }

    void CalculationLayer::UploadMolecule() {
        if (mMolecule == nullptr)
            return;

        float *data = new float[mMolecule->GetNoAtoms() * 7];
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

        mMoleculeBuffer->StoreData(data, mMolecule->GetAtoms().size() * 7);

        delete[] data;
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

    glm::mat3 CalculationLayer::GenerateBondRotMatrix(const glm::vec3 &bondDirNorm, const glm::vec3 &startingDirNorm) {
        glm::vec3 v = glm::cross(bondDirNorm, startingDirNorm);
        float c = glm::dot(bondDirNorm, startingDirNorm);

        glm::mat3 vSkew = {
            0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0
        };

        return glm::mat3(1.0f) + vSkew + vSkew * vSkew / (1 + c);
    }

    void CalculationLayer::ExecuteCalculation() {
        if (mMolecule == nullptr)
            return;

        mPosTrisurf.Reset();
        mNegTrisurf.Reset();

        mOrbitalsOrdered = false;
        mSelectedMO = -1;

        mCalcBasiSet = mBasisSets[mSelectedBasisSet];
        mMolecule->MakeOrbitals(mCalcBasiSet);

        mConsoleOutput.clear();
        mConsoleOutput.str("");

        mCalculationThread.Start([this]() {
            if (mCorrelationMethod != CorrelationMethod::None) {
                mSCFMethod = SCFMethod::RHF;

                //Run SCF as RHF if a SCF determinant is not available
                if (mRHFResult == nullptr || !mUseLastSCF)
                    mRHFResult = std::make_shared<RHFResult>(
                        std::move(SCESS::HartreeFock::Restricted(*mMolecule, 1E-5, 100, &mConsoleOutput,
                                                                 mCalcBasiSet)));

                switch (mCorrelationMethod) {
                    case CorrelationMethod::MP2:
                        mMP2Energy = SCESS::HartreeFock::MP2ClosedShell(*mRHFResult, &mConsoleOutput);
                        break;
                    case CorrelationMethod::CCSD: {
                        auto ccResult = SCESS::CoupledCluster::CCSD(*mRHFResult, &mConsoleOutput);
                        if (ccResult.has_value()) {
                            mCCSDResult = std::make_shared<CCSDResult>(ccResult.value());
                            mMP2Energy = ccResult->mp2CorrelationEnergy;
                        }
                        break;
                    }

                    default:
                        break;
                }
            } else {
                switch (mSCFMethod) {
                    case SCFMethod::RHF:
                        mRHFResult = std::make_shared<RHFResult>(
                            std::move(SCESS::HartreeFock::Restricted(*mMolecule, 1E-5, 100, &mConsoleOutput,
                                                                     mCalcBasiSet)));
                        break;
                    case SCFMethod::UHF:
                        mUHFResult = std::make_shared<UHFResult>(
                            std::move(SCESS::HartreeFock::Unrestricted(*mMolecule, 1E-5, 100, &mConsoleOutput,
                                                                       mCalcBasiSet)));
                        break;
                    case SCFMethod::DFT:
                        mDFTEnergy = SCESS::DFT::RunDFT(*mMolecule, &mConsoleOutput, mCalcBasiSet);
                        break;
                }
            }

            std::cout << mConsoleOutput.str();
            mUseLastSCF = true;
            mCalculationThread.MarkAsFinished();
        });
    }

    void CalculationLayer::DrawProgressScreen() {
        float dpiFactor = ImGui::GetWindowDpiScale();
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Once, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize({600.0f * dpiFactor, 300.0f * dpiFactor}, ImGuiCond_Once);
        if (ImGui::BeginPopupModal(PROG_SCREEN_NAME)) {
            ImSpinner::SpinnerAtom("##atomSpinner", 20.0f * dpiFactor, 2.5f);
            ImGui::SameLine();
            ImGui::Text("Running %s single-point energy. Please wait...", [this]() {
                if (mCorrelationMethod != CorrelationMethod::None)
                    switch (mCorrelationMethod) {
                        case CorrelationMethod::MP2:
                            return "MP2";
                        case CorrelationMethod::CCSD:
                            return "CCSD";
                        default:
                            return "";
                    };

                switch (mSCFMethod) {
                    case SCFMethod::RHF:
                        return "RHF";
                    case SCFMethod::UHF:
                        return "UHF";
                    case SCFMethod::DFT:
                        return "DFT";
                    default:
                        return "";
                };
            }());

            const std::string &consoleStr = mConsoleOutput.str();
            ImGui::InputTextMultiline("##consoleOutput", (char *) consoleStr.c_str(), consoleStr.length(),
                                      {0, 0}, ImGuiInputTextFlags_ReadOnly);

            if (mCalculationThread.HasFinished())
                ImGui::CloseCurrentPopup();

            ImGui::EndPopup();
        }
    }

    void CalculationLayer::OrderOrbitals() {
        auto sortCoeffs = [this](Eigen::VectorXd &Energies, Eigen::MatrixXd &Coeffs,
                                 const Eigen::VectorXi &Occupation) {
            std::unordered_map<double, uint32_t> sortMap;
            for (Eigen::Index i = 0; i < Energies.rows(); i++) {
                sortMap.insert({Energies[i], i});
            }
            std::sort(Energies.begin(), Energies.end());

            Eigen::MatrixXd sortedCoeffs = Eigen::MatrixXd::Zero(Coeffs.rows(), Coeffs.rows());

            mMOLabels.reserve(Energies.rows());
            for (Eigen::Index i = 0; i < Energies.rows(); i++) {
                uint32_t index = sortMap.at(Energies[i]);
                sortedCoeffs.col(i) = Coeffs.col(index);

                std::stringstream ss;
                ss.precision(3);
                ss << "MO " << i + 1 << "; E=" << Energies[i] << " hartrees. ";

                if (i + 1 < Occupation.rows())
                    if (Occupation(i + 1) == 0 && Occupation(i) > 0)
                        ss << "HOMO";
                if (i != 0)
                    if (Occupation(i - 1) > 0 && Occupation(i) == 0)
                        ss << "LUMO";

                mMOLabels.push_back(ss.str());
            }
            Coeffs = sortedCoeffs;
        };

        mMOLabels.clear();
        if (mSCFMethod == SCFMethod::RHF)
            sortCoeffs(mRHFResult->Energies, mRHFResult->Coeffs, mRHFResult->Occupation);
        else {
            sortCoeffs(mUHFResult->EnergiesAlpha, mUHFResult->CoeffsAlpha, mUHFResult->OccupationAlpha);
            sortCoeffs(mUHFResult->EnergiesBeta, mUHFResult->CoeffsBeta, mUHFResult->OccupationBeta);
        }
    }

    void CalculationLayer::DrawQuantumCalcPanel() {
        if (ImGui::Begin("Quantum Calculation")) {
            if (mMolecule != nullptr) {
                ImGui::SeparatorText("SCF Methods");

                ImGui::Indent(100.0);
                ImGui::SeparatorText("Hartree-Fock");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("Restricted Hartree Fock", reinterpret_cast<int *>(&mSCFMethod), 0)) {
                    mSCFMethod = SCFMethod::RHF;
                }

                if (ImGui::RadioButton("Unrestricted Hartree Fock", reinterpret_cast<int *>(&mSCFMethod), 1)) {
                    mCorrelationMethod = CorrelationMethod::None;
                    mSCFMethod = SCFMethod::UHF;
                }
                ImGui::Unindent(100.0);

                ImGui::SeparatorText("Density Functional Theory");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("LYP", reinterpret_cast<int *>(&mSCFMethod), 2)) {
                    mCorrelationMethod = CorrelationMethod::None;
                    mSCFMethod = SCFMethod::DFT;
                }
                ImGui::Unindent(200.0);


                ImGui::SeparatorText("Electron Correlation Methods");
                ImGui::Indent(200.0);

                if (ImGui::RadioButton("None", reinterpret_cast<int *>(&mCorrelationMethod), 0))
                    mCorrelationMethod = CorrelationMethod::None;

                ImGui::Unindent(100.0);
                ImGui::SeparatorText("Moller-Plesset");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("MP2", reinterpret_cast<int *>(&mCorrelationMethod), 1)) {
                    mCorrelationMethod = CorrelationMethod::MP2;
                    mSCFMethod = SCFMethod::RHF;
                }
                ImGui::Unindent(100.0);

                ImGui::SeparatorText("Coupled Cluster");
                ImGui::Indent(100.0);
                if (ImGui::RadioButton("CCSD", reinterpret_cast<int *>(&mCorrelationMethod), 2)) {
                    mCorrelationMethod = CorrelationMethod::CCSD;
                    mSCFMethod = SCFMethod::RHF;
                }
                ImGui::Unindent(200.0);

                ImGui::SeparatorText("Basis Sets");
                ImGui::Indent(100.0);
                int lastSelectedBasisSet = mSelectedBasisSet;
                for (int i = 0; i < mBasisSets.size(); i++)
                    ImGui::RadioButton(mBasisSets[i].c_str(), &mSelectedBasisSet, i);

                mUseLastSCF = lastSelectedBasisSet != mSelectedBasisSet;
                ImGui::Unindent(100.0);

                ImGui::NewLine();
                ImGui::PushStyleColor(ImGuiCol_Button, {0.0f, 0.4f, 0.0f, 1.0f});
                if (ImGui::Button("Execute Calculation", ImVec2{-FLT_MIN, 0.0f})) {
                    ExecuteCalculation();
                }
                ImGui::PopStyleColor();
            }
        }
        ImGui::End();
    }

    void CalculationLayer::DrawResultsPanel() {
        if (ImGui::Begin("Results")) {
            if (mCalculationThread.HasFinished()) {
                if (mCorrelationMethod != CorrelationMethod::None && ImGui::CollapsingHeader(
                        "Electron Correlation Result", ImGuiTreeNodeFlags_DefaultOpen)) {
                    switch (mCorrelationMethod) {
                        case CorrelationMethod::CCSD: {
                            if (!mCCSDResult)break;
                            ImGui::SeparatorText("Coupled Cluster Singles & Doubles");
                            ImGui::Indent();
                            ImGui::Text("E(CCSD) = %.8f hartrees", mCCSDResult->ccsdCorrelationEnergy);
                            ImGui::TextColored({0.0, 1.0, 0.0, 1.0},
                                               "E(Total) = E(Hartree-Fock) + E(CCSD) = %.8f hartrees",
                                               mCCSDResult->ccsdCorrelationEnergy + mRHFResult->Energy);
                            ImGui::Unindent();
                        }

                        case CorrelationMethod::MP2: {
                            if (!mMP2Energy.has_value())break;
                            ImGui::SeparatorText("Moller-Plesset Order 2");
                            ImGui::Indent();
                            ImGui::Text("E(MP2) = %.8f hartrees", mMP2Energy.value());
                            ImGui::TextColored({0.0, 1.0, 0.0, 1.0},
                                               "E(Total) = E(Hartree-Fock) + E(MP2) = %.8f hartrees",
                                               mMP2Energy.value() + mRHFResult->Energy);
                            ImGui::Unindent();
                            break;
                        }


                        case CorrelationMethod::None:
                        default:
                            break;
                    }
                }

                if (ImGui::CollapsingHeader("SCF Result", ImGuiTreeNodeFlags_DefaultOpen)) {
                    if (mSCFMethod == SCFMethod::RHF || mSCFMethod == SCFMethod::UHF) {
                        if (!mOrbitalsOrdered) {
                            OrderOrbitals();
                            mOrbitalsOrdered = true;
                        }

                        ImGui::BulletText("Basis Set: %s", mCalcBasiSet.c_str());
                        ImGui::BulletText("Total Hartree Fock Energy: %.4f hartrees",
                                          mSCFMethod == SCFMethod::RHF ? mRHFResult->Energy : mUHFResult->Energy);

                        ImGui::Indent();
                        if (ImGui::CollapsingHeader("Orbitals", ImGuiTreeNodeFlags_DefaultOpen)) {
                            if (mSelectedMO != -1) {
                                if (mSCFMethod == SCFMethod::RHF) {
                                    ImGui::TextColored({1.0f, 1.0f, 0.0f, 1.0f}, "Previewing MO %i", mSelectedMO + 1);
                                    if (mSelectedMO + 1 < mRHFResult->Occupation.rows())
                                        if (mRHFResult->Occupation(mSelectedMO + 1) == 0 && mRHFResult->Occupation(
                                                mSelectedMO) > 0)
                                            ImGui::TextColored({1.0f, 0.0f, 0.0f, 1.0f}, "HOMO");
                                    if (mSelectedMO - 1 > 0)
                                        if (mRHFResult->Occupation(mSelectedMO - 1) > 0 && mRHFResult->Occupation(
                                                mSelectedMO) == 0)
                                            ImGui::TextColored({0.0f, 0.5f, 0.5f, 1.0f}, "LUMO");
                                    ImGui::BulletText("Energy: %.3f hartrees", mRHFResult->Energies(mSelectedMO));
                                    ImGui::BulletText("Occupation: %i e-", mRHFResult->Occupation(mSelectedMO));
                                } else {
                                    const auto &energies =
                                            mSelectedSpin == 1 ? mUHFResult->EnergiesAlpha : mUHFResult->EnergiesBeta;
                                    const auto &occupation = mSelectedSpin == 1
                                                                 ? mUHFResult->OccupationAlpha
                                                                 : mUHFResult->OccupationBeta;
                                    ImGui::TextColored({1.0f, 1.0f, 0.0f, 1.0f}, "Previewing MO %i", mSelectedMO + 1);
                                    if (mSelectedMO + 1 < occupation.rows())
                                        if (occupation(mSelectedMO + 1) == 0 && occupation(mSelectedMO) > 0)
                                            ImGui::TextColored({1.0f, 0.0f, 0.0f, 1.0f}, "HOMO");
                                    if (mSelectedMO - 1 > 0)
                                        if (occupation(mSelectedMO - 1) > 0 && occupation(mSelectedMO) == 0)
                                            ImGui::TextColored({0.0f, 0.5f, 0.5f, 1.0f}, "LUMO");
                                    ImGui::BulletText("Spin: %s", mSelectedSpin == 1 ? "Alpha (up)" : "Beta (down)");
                                    ImGui::BulletText("Energy: %.3f hartrees", energies(mSelectedMO));
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
                                if (mSCFMethod == SCFMethod::RHF)
                                    ImGui::ListBox("##listbox1", &selectedOption, &orbitalLabels[0],
                                                   (int) mMOLabels.size());
                                else {
                                    float width = ImGui::GetColumnWidth() * 0.5f;

                                    ImGui::Text("Alpha (spin up) Orbitals");
                                    ImGui::SameLine(width, 10.0f);
                                    ImGui::Text("Beta (spin down) Orbitals");
                                    ImGui::PushItemWidth(width);
                                    ImGui::ListBox("##listbox1", &selectedOption, &orbitalLabels[0],
                                                   (int) mUHFResult->EnergiesAlpha.rows());
                                    ImGui::SameLine();
                                    ImGui::ListBox("##listbox2", &selectedOption2,
                                                   &orbitalLabels[mUHFResult->EnergiesAlpha.rows()],
                                                   (int) mUHFResult->EnergiesBeta.rows());
                                    ImGui::PopItemWidth();
                                }
                            }

                            if (selectedOption != -1 || selectedOption2 != -1) {
                                if (mSCFMethod == SCFMethod::RHF)
                                    mSelectedMO = selectedOption;
                                else {
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
                                UploadTriSurface(*mPosTrisurfVAO, mPosTrisurf);
                                // if (mNegTrisurf.triangles.size() != 0)
                                UploadTriSurface(*mNegTrisurfVAO, mNegTrisurf);

                                mIsoWasNull = mPosTrisurf.GetTriangles().size() == 0 && mNegTrisurf.GetTriangles().
                                              size() == 0;

                                mPosTrisurf.Reset();
                                mNegTrisurf.Reset();
                                mIsosurfaceUploaded = true;
                            }

                            if (mPosTrisurf.IsRunning() || mNegTrisurf.IsRunning()) {
                                ImSpinner::SpinnerTwinBall("spinner1", 16.0, 16.0f, 3.0, 6.0, {1.0f, 0.0f, 0.0f, 1.0f},
                                                           ImSpinner::white, 3.0f, 1);
                                ImGui::SameLine();
                                ImGui::Text("Generating isosurface...");
                                if (ImGui::Button("Cancel", {100, 0})) {
                                    mPosTrisurf.Reset();
                                    mNegTrisurf.Reset();
                                }
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
                                    ImGui::DragInt("##isoResolutionVal", &mIsosurfaceResolution, 1.0f, 10, 1000);

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
                    } else if (mSCFMethod == SCFMethod::DFT) {
                        ImGui::TextColored({0.0, 1.0, 0.0, 1.0}, "Density Functional calculation finished!");
                        ImGui::BulletText("Total DFT Energy : %.4f hartrees", mDFTEnergy);
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
            ImGui::Indent(100.0f);
            int lastSelectedMol = mSelectedMolecule;
            for (int i = 0; i < IM_ARRAYSIZE(moleculeNames); i++)
                ImGui::RadioButton(moleculeNames[i], &mSelectedMolecule, i);

            if (ImGui::RadioButton("Load from XYZ File...", &mSelectedMolecule, IM_ARRAYSIZE(moleculeNames) + 1)) {
                if (!SelectMoleculeFromFile())
                    mSelectedMolecule = lastSelectedMol;
            } else if (mSelectedMolecule != lastSelectedMol && mSelectedMolecule != -1)
                SelectMolecule();

            ImGui::Unindent(100.0f);
            ImGui::NewLine();
            ImGui::SeparatorText("Selected Molecule");
            ImGui::Indent(100.0f);


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
            const std::string &consoleStr = mConsoleOutput.str();
            ImGui::InputTextMultiline("##consoleOutput", (char *) consoleStr.c_str(), consoleStr.length(),
                                      {-FLT_MIN, -FLT_MIN}, ImGuiInputTextFlags_ReadOnly);
        }
        ImGui::End();
    }

    void CalculationLayer::OnImGuiRender() {
        DrawMolSelectPanel();

        DrawProgressScreen();
        if (mCalculationThread.IsRunning())
            ImGui::OpenPopup(PROG_SCREEN_NAME);

        DrawQuantumCalcPanel();
        DrawResultsPanel();
        DrawConsolePanel();
    }

    void CalculationLayer::OnEvent(Event &e) {
    }
}
