#pragma once
#include<pch.h>
#include"Display.h"


namespace SCESS
{
	constexpr size_t OPENGL_MAJOR_VERSION = 4;
	constexpr size_t OPENGL_MINOR_VERSION = 6;

	class GLFWDisplay : public SCESS::Display
	{
	public:
		GLFWDisplay(const DisplayProps& properties);
		virtual ~GLFWDisplay();

		virtual void Show() const { glfwShowWindow(m_WindowID); }
		virtual void Hide() const { glfwHideWindow(m_WindowID); }
		virtual void Focus() const { glfwFocusWindow(m_WindowID); }

		virtual void Refresh() override;
		virtual size_t GetWidth() const override { return m_Data.width; }
		virtual size_t GetHeight() const override { return m_Data.height; }

		virtual inline void BindEventCallback(const EventCallbackFnc& callback) override { m_Data.EventCallback = callback; }

		virtual bool CloseRequested() const override { return glfwWindowShouldClose(m_WindowID); }
		virtual double GetTime() const override { return glfwGetTime(); }

		virtual void* GetAPIDisplay() const override { return (void*)m_WindowID; };

	private:
		void Initilise(const DisplayProps& properties,
			const size_t majorVersion = OPENGL_MAJOR_VERSION,
			const size_t minorVersion = OPENGL_MINOR_VERSION);
		static void PrintErrorCallback(int errorCode, const char* message);

		GLFWwindow* m_WindowID = nullptr;

		struct GLFWData
		{
			std::string title;
			size_t width = 0, height = 0;

			EventCallbackFnc EventCallback;
		};

		GLFWData m_Data;
	};
}
