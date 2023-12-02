#include<pch.h>
#include "GLFWDisplay.h"
#include"Events/Events.h"
#include"Events/WindowEvents.h"
#include"Events/KeyEvents.h"
#include"Events/MouseEvents.h"
//
//

static bool isGLFWInitilised = false;

//
SCESS::GLFWDisplay::GLFWDisplay(const DisplayProps& properties)
{
	this->Initilise(properties);
}
//
SCESS::GLFWDisplay::~GLFWDisplay()
{
	glfwDestroyWindow(m_WindowID);
}
//
void SCESS::GLFWDisplay::Refresh()
{
	
	glfwSwapBuffers(m_WindowID);
	glfwPollEvents();
}
//
void SCESS::GLFWDisplay::Initilise(const DisplayProps& properties, const size_t majorVersion, const size_t minorVersion)
{
	if (!isGLFWInitilised)
	{
		glfwSetErrorCallback(&SCESS::GLFWDisplay::PrintErrorCallback);
		int errorCode = glfwInit();
		if (!errorCode)
		{
			glfwTerminate();
			return;
		}

		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, majorVersion);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, minorVersion);
	}

	m_Data.width = properties.width;
	m_Data.height = properties.height;
	m_Data.title = properties.title;

	m_WindowID = glfwCreateWindow(properties.width, properties.height, properties.title.c_str(), nullptr, nullptr);
	glfwSetWindowUserPointer(m_WindowID, &this->m_Data);
	glfwMakeContextCurrent(m_WindowID);
	glfwShowWindow(m_WindowID);

	//Set GLFW callbacks
	glfwSetWindowSizeCallback(m_WindowID, [](GLFWwindow* window, int width, int height)
		{
			GLFWData& data = *(GLFWData*)glfwGetWindowUserPointer(window);
			data.width = width;
			data.height = height;

			WindowResizeEvent resizeEvt(width, height);
			data.EventCallback(resizeEvt);
		});
	glfwSetWindowCloseCallback(m_WindowID, [](GLFWwindow* window) {
		GLFWData& data = *(GLFWData*)glfwGetWindowUserPointer(window);

		WindowCloseEvent closeEvent;
		data.EventCallback(closeEvent);
		});
	glfwSetKeyCallback(m_WindowID, [](GLFWwindow* window, int key, int scancode, int action, int mods) {
			GLFWData& data = *(GLFWData*)glfwGetWindowUserPointer(window);
			switch (action)
			{
				case GLFW_PRESS:
				{
					KeyPressedEvent event(key, 0);
					data.EventCallback(event);
					break;
				}

				case GLFW_RELEASE:
				{
					KeyReleasedEvent event(key);
					data.EventCallback(event);
					break;
				}

				case GLFW_REPEAT:
				{
					KeyPressedEvent event(key, 1);
					data.EventCallback(event);
					break;
				}
			}

		});
	glfwSetCharCallback(m_WindowID, [](GLFWwindow* window, unsigned int keyCode) {
		GLFWData& data = *(GLFWData*)glfwGetWindowUserPointer(window);
		KeyTypedEvent event(keyCode);
		data.EventCallback(event);
		});
	glfwSetMouseButtonCallback(m_WindowID, [](GLFWwindow* window, int button, int action, int mods) {
			GLFWData& data = *(GLFWData*)glfwGetWindowUserPointer(window);

			auto GLFWToMouseButton = [](int button)
			{
				switch (button)
				{
				case GLFW_MOUSE_BUTTON_LEFT:
					return MouseButton::Left;
				case GLFW_MOUSE_BUTTON_MIDDLE:
					return MouseButton::Centre;
				case GLFW_MOUSE_BUTTON_RIGHT:
					return MouseButton::Right;
				default:
					return MouseButton::None;
				}
			};
			
			switch (action)
			{
				case GLFW_PRESS:
				{
					MouseButtonPressedEvent event(GLFWToMouseButton(button));
					data.EventCallback(event);
					break;
				}
				case GLFW_RELEASE:
				{
					MouseButtonReleasedEvent event(GLFWToMouseButton(button));
					data.EventCallback(event);
					break;
				}
			}
		});
	glfwSetScrollCallback(m_WindowID, [](GLFWwindow* window, double xoffset, double yoffset)
		{
			GLFWData& data = *(GLFWData*)glfwGetWindowUserPointer(window);

			MouseScrollEvent event((float)xoffset, (float)yoffset);
			data.EventCallback(event);
		});
	glfwSetCursorPosCallback(m_WindowID, [](GLFWwindow* window, double xPos, double yPos)
		{
			GLFWData& data = *(GLFWData*)glfwGetWindowUserPointer(window);

			MouseMoveEvent event({ (float)xPos, (float)yPos });
			data.EventCallback(event);
		});

}
//
void SCESS::GLFWDisplay::PrintErrorCallback(int errorCode, const char* message)
{
	std::cerr << "[GLFW Error]: Code: " << errorCode << " Message: " << message << '\r' << '\n';
}
//
