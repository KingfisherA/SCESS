#pragma once
#include<pch.h>
#include"imgui.h"
#include"Utilities/Memory.h"
#include"Display.h"
#include"Scene/Camera.hpp"
#include"Events/WindowEvents.h"

#include"Layers/LayerStack.h"
#include"GUI/ImGuiLayer.h"
#include"GUI/GraphLayer.h"

namespace SCESS
{
	class KeyPressedEvent;
	class KeyReleasedEvent;
	class MouseMoveEvent;
	class MouseButtonPressedEvent;
	class MouseButtonReleasedEvent;
	class MouseScrollEvent;

	class Framebuffer;

	class Application
	{
	public:
		Application();
		~Application();

		int Run();

		void OnEvent(Event& e);

		void PushLayer(Layer* layer);
		void PushOverlay(Layer* overlay);

		void UpdateProjectionMatrix();

		inline glm::mat4& GetProjectionMatrix() { return projectionMatrix; }

		static inline Application& Get() { return *s_Instance; }
		const inline Display& GetDisplay() const { return *m_Window; }
		constexpr double GetDeltaTime() const { return m_DeltaTime; }
		constexpr double GetTime() const { return m_Time; }
		const Camera& GetCamera() const { return mCam; }
		Camera& GetCamera() { return mCam; }

	private:
		bool OnWindowClose(WindowCloseEvent& e);
		bool OnKeyPressed(KeyPressedEvent& e);
		bool OnKeyReleased(KeyReleasedEvent& e);
		bool OnMouseButtonPressed(MouseButtonPressedEvent& e);
		bool OnMouseButtonReleased(MouseButtonReleasedEvent& e);
		bool OnMouseMove(MouseMoveEvent& e);
		bool OnMouseScroll(MouseScrollEvent& e);
		bool m_isViewportFocused = false;

		void RenderViewport();
		void RenderCameraPanel();
		void RenderSceneInformation();

		std::string m_VersionStr, m_VendorStr, m_ShadingLanguageStr;

		bool m_IsRunning = true;
		uRef<Display> m_Window;
		sRef<Framebuffer> mainFramebuffer = nullptr;
		ImVec2 m_viewPortPanelSize = { 0,0 };
		ImVec2 m_lastViewportSize = { 0,0 };

		SCESS::ImGuiLayer* m_ImGuiLayer = nullptr;
		SCESS::GraphLayer* mGraphLayer = nullptr;
		LayerStack m_LayerStack;

		double m_Time = 0.0;
		double m_DeltaTime = 0.0;
		glm::vec2 m_MouseLastPos;
		int m_MouseButtonDown = -1;
		int m_KeyDown = -1;

		Camera mCam;

		static Application* s_Instance;
		static glm::mat4 projectionMatrix;
	};
}
