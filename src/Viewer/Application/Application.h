#pragma once

#include "Dahlia/Application/GraphicsApp.h"
#include "Dahlia/Utilities/Memory.h"
#include "Dahlia/GUI/Display.h"
#include "Dahlia/Scene/Camera3D.h"
#include "Dahlia/Events/WindowEvents.h"
#include "Dahlia/Layers/LayerStack.h"
#include "Dahlia/GUI/ImGuiLayer.h"

#include<boost/filesystem.hpp>

namespace DH {
    class KeyPressedEvent;
    class KeyReleasedEvent;
    class MouseMoveEvent;
    class MouseButtonPressedEvent;
    class MouseButtonReleasedEvent;
    class MouseScrollEvent;
    class CalculationLayer;
    class Framebuffer;
    class Texture2D;
}

namespace SCESS::Viewer {
    class CalculationLayer;

    class Application : public DH::GraphicsApp {
    public:
        Application();
        ~Application() override = default;


        [[nodiscard]] static Application &Get() { return *sInstance; }

        [[nodiscard]] const DH::Camera3D &GetCamera() const { return mCam; }
        [[nodiscard]] DH::Camera3D &GetCamera() { return mCam; }

        void UpdateProjectionMatrix();
        [[nodiscard]] glm::mat4 &GetProjectionMatrix() { return mProjectionMatrix; }
        [[nodiscard]] const glm::mat4 &GetProjectionMatrix() const { return mProjectionMatrix; }

    private:
        void OnViewportResize() override;
        void OnRender() override;
        void OnImGuiRender() override;
        bool OnInit() override;
        void OnDestroy() override;

        bool OnMouseMove(DH::MouseMoveEvent &e) override;
        bool OnMouseScroll(DH::MouseScrollEvent &e) override;

        void RenderMenuBar();
        void RenderSplashScreen();
        void RenderCameraPanel();
        void RenderSceneInformation();

        CalculationLayer *mCalcLayer = nullptr;

        DH::Camera3D mCam;

        static Application *sInstance;
        glm::mat4 mProjectionMatrix;

        sRef<DH::Texture2D> mAppIcon;

        bool mSplashScreenOpen = true;
    };
}
