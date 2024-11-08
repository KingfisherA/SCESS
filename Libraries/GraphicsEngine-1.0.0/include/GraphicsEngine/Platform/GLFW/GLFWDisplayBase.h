#pragma once

#include <GLFW/glfw3.h>
#include "GraphicsEngine/GUI/Display.h"

#include<string>

namespace GEng {
    class GLFWDisplayBase : public Display {
    public:
        ~GLFWDisplayBase() override;

        void Show() const override { glfwShowWindow(mWindow); }
        void Hide() const override { glfwHideWindow(mWindow); }
        void Focus() const override { glfwFocusWindow(mWindow); }

        void Refresh() override;
        [[nodiscard]] size_t GetWidth() const override { return mWidth; }
        [[nodiscard]] size_t GetHeight() const override { return mHeight; }
        [[nodiscard]] glm::vec2 GetContentScaling() const override;

        void BindEventCallback(const EventCallbackFnc &callback) override {
            mEventCallback = callback;
        }

        [[nodiscard]] bool CloseRequested() const override { return glfwWindowShouldClose(mWindow); }
        [[nodiscard]] double GetTime() const override { return glfwGetTime(); }

        [[nodiscard]] void *GetAPIDisplay() const override { return (void *) mWindow; };

    protected:
        GLFWDisplayBase() = default;

        void CreateDisplay(const DisplayProps &props);

        static void OnWindowResize(GLFWwindow *window, int width, int height);
        static void OnWindowClose(GLFWwindow *window);
        static void OnKeyPressed(GLFWwindow *window, int key, int scancode, int action, int mods);
        static void OnCharPressed(GLFWwindow *window, unsigned int keyCode);
        static void OnMousePressed(GLFWwindow *window, int button, int action, int mods);
        static void OnScroll(GLFWwindow *window, double xoffset, double yoffset);
        static void OnCursorPos(GLFWwindow *window, double xPos, double yPos);

        GLFWwindow *mWindow = nullptr;

        size_t mWidth = 0, mHeight = 0;
        std::string mTitle;

        EventCallbackFnc mEventCallback;
    };
} // GEng
