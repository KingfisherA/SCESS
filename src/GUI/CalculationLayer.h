#pragma once

#include"GraphicsEngine/Layers/Layer.h"
#include"GraphicsEngine/3D/VertexArray.h"
#include"GraphicsEngine/3D/ShaderProgram.h"
#include"GraphicsEngine/Utilities/Threading.h"

#include"Chemistry/Quantum/Molecule.h"
#include "Maths/Isosurfaces.h"


namespace GEng {
    class Camera3D;
    class ShaderProgram;
    class VertexArray;
    class VertexBuffer;
}

namespace SCESS {
    class GaussianCollection;
    struct RHFResult;
    struct UHFResult;
    class TriSurf;

    enum class CalcMethod {
        RHF = 0,
        UHF,
        DFT
    };

    enum class MoleculeStyle:int {
        SpaceFill = 0,
        BallStick = 1
    };

    class CalculationLayer : public GEng::Layer {
    public:
        CalculationLayer() : Layer("Graph Layer") {
        }

        //virtual ~CalculationLayer();

        virtual void OnAttach() override;
        virtual void OnDetach() override;
        virtual void OnUpdate(const sRef<GEng::Framebuffer> &screenBuffer, const GEng::Camera3D &cam,
                              const glm::mat4 &projectionMatrix, const glm::vec2 &nearFarPlanes) override;
        virtual void OnImGuiRender() override;
        virtual void OnEvent(GEng::Event &e) override;

    private:
        void UploadTriSurface(GEng::VertexArray &vArray, const TriSurf &tri);
        void RenderMOIsosurface();
        void SelectMolecule();
        void UploadMolecule();
        void SelectMoleculeFromFile();

        glm::mat3 GenerateBondRotMatrix(const glm::vec3 &bondDirNorm, const glm::vec3 &startingDirNorm);

        void ExecuteCalculation();
        void DrawProgressScreen();
        void OrderOrbitals();

        std::vector<std::string> mBasisSets;
        int mSelectedBasisSet = 0;

        int mSelectedMolecule = -1;
        std::string mMoleculeName;
        uRef<Molecule> mMolecule;
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
