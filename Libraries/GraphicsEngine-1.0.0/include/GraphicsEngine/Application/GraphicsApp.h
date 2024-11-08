#pragma once

#include "GraphicsEngine/GEng_API.h"
#include <string>
#include <ImGui/imgui.h>
#include <glm/glm.hpp>

#include "GraphicsEngine/Utilities/Memory.h"
#include "GraphicsEngine/Layers/LayerStack.h"
#include "GraphicsEngine/GUI/Display.h"
#include "GraphicsEngine/3D/Framebuffer.h"
#include "GraphicsEngine/GUI/ImGuiLayer.h"
#include "GraphicsEngine/GUI/Text.h"
#include "GraphicsEngine/Resources/Resources.h"

namespace GEng {
    class Event;
    class Layer;
    class WindowCloseEvent;
    class KeyPressedEvent;
    class KeyReleasedEvent;
    class MouseButtonPressedEvent;
    class MouseButtonReleasedEvent;
    class MouseMoveEvent;
    class MouseScrollEvent;
    class WindowResizeEvent;
    class Framebuffer;
    class ImGuiLayer;

    using size_t = std::size_t;

#define GRAPHICS_APP_INIT_IMGUI_CTX ImGui::SetCurrentContext(m_ImGuiLayer->GetContext().ctx);\
    ImGui::SetAllocatorFunctions(m_ImGuiLayer->GetContext().allocFunc, m_ImGuiLayer->GetContext().freeFunc, m_ImGuiLayer->GetContext().userData);

    GENG_API class GraphicsApp {
    public:
        GraphicsApp(const char *displayTitle,
                    size_t width, size_t height, const FramebufferSpecification &viewportSpec,
                    const Font &font = Res::FontRes::GetDefaultFont());
        virtual ~GraphicsApp() = default;

        bool Init();
        void Run();

        void OnEvent(Event &e);

        void PushLayer(Layer *layer);
        void PushOverlay(Layer *overlay);

        Framebuffer &GetViewportFramebuffer() { return *mViewportFramebuffer; }
        const inline Display &GetDisplay() const { return *m_Window; }
        constexpr double GetDeltaTime() const { return m_DeltaTime; }
        constexpr double GetTime() const { return m_Time; }

    protected:
        virtual void OnViewportResize() {
        }

        virtual void OnRender() {
        }

        virtual void OnImGuiRender() {
        }

        virtual bool OnInit() { return true; }
        virtual void OnDestroy();

        virtual bool OnWindowResized(WindowResizeEvent &e);
        virtual bool OnWindowClose(WindowCloseEvent &e);
        virtual bool OnKeyPressed(KeyPressedEvent &e);
        virtual bool OnKeyReleased(KeyReleasedEvent &e);
        virtual bool OnMouseButtonPressed(MouseButtonPressedEvent &e);
        virtual bool OnMouseButtonReleased(MouseButtonReleasedEvent &e);
        virtual bool OnMouseMove(MouseMoveEvent &e);
        virtual bool OnMouseScroll(MouseScrollEvent &e);
        bool mIsViewportFocused = false;
        bool mRenderViewport = true;

        void RenderViewport();
        void RenderCameraPanel();
        void RenderSceneInformation();

        std::string m_VersionStr, m_VendorStr, m_ShadingLanguageStr;

        bool m_IsRunning = true;
        uRef<Display> m_Window;
        sRef<Framebuffer> mViewportFramebuffer = nullptr;
        ImVec2 m_viewPortPanelSize = {0, 0};
        ImVec2 m_lastViewportSize = {0, 0};

        ImGuiLayer *m_ImGuiLayer = nullptr;
        LayerStack m_LayerStack;

        double m_Time = 0.0;
        double m_DeltaTime = 0.0;

        glm::dvec2 mCurrentMousePos;
        glm::dvec2 mLastMousePos;
        glm::dvec2 mDeltaMousePos;
        int m_MouseButtonDown = -1;
        int m_KeyDown = -1;
        int m_KeyMods = 0;

        double mDeltaScroll = 0.0;
        double mLastScrollTime = 0.0;

        std::string mViewportName;

        Font mFont;

    private:
        GENG_API struct InitParams {
            const char *displayTitle = nullptr;
            size_t width = 0, height = 0;
            FramebufferSpecification viewportSpec;
        };

        uRef<InitParams> mInitParams = nullptr;
    };
}
