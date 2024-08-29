#pragma once

#include<pch.h>

#include "GraphicsEngine/Application/GraphicsApp.h"
#include"GraphicsEngine/Utilities/Memory.h"
#include"GraphicsEngine/GUI/Display.h"
#include"GraphicsEngine/Scene/Camera3D.h"
#include"GraphicsEngine/Events/WindowEvents.h"
#include"GraphicsEngine/Layers/LayerStack.h"
#include"GraphicsEngine/GUI/ImGuiLayer.h"

#include<boost/filesystem.hpp>

namespace GEng {
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

namespace SCESS {
    class CalculationLayer;

    class Application : public GEng::GraphicsApp {
    public:
        Application();
        ~Application() override = default;


        [[nodiscard]] static Application &Get() { return *sInstance; }

        [[nodiscard]] const GEng::Camera3D &GetCamera() const { return mCam; }
        [[nodiscard]] GEng::Camera3D &GetCamera() { return mCam; }

        void UpdateProjectionMatrix();
        [[nodiscard]] glm::mat4 &GetProjectionMatrix() { return mProjectionMatrix; }
        [[nodiscard]] const glm::mat4 &GetProjectionMatrix() const { return mProjectionMatrix; }

        [[nodiscard]] const boost::filesystem::path &GetAssetsPath() const { return mAssetsPath; }
        [[nodiscard]] boost::filesystem::path GetResourcePath(const std::string& resource) {
            return boost::filesystem::path(mAssetsPath).append(resource);
        };

    private:
        void OnViewportResize() override;
        void OnRender() override;
        void OnImGuiRender() override;
        bool OnInit() override;
        void OnDestroy() override;

        bool OnMouseMove(GEng::MouseMoveEvent &e) override;
        bool OnMouseScroll(GEng::MouseScrollEvent &e) override;

        void RenderMenuBar();
        void RenderSplashScreen();
        void RenderCameraPanel();
        void RenderSceneInformation();

        CalculationLayer *mCalcLayer = nullptr;

        GEng::Camera3D mCam;

        static Application *sInstance;
        glm::mat4 mProjectionMatrix;

        sRef<GEng::Texture2D> mAppIcon;

        boost::filesystem::path mAssetsPath = "../../../Assets";

        bool mSplashScreenOpen = true;
    };
}
