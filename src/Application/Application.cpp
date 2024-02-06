#include<pch.h>
#include"Application.h"
#include"imgui.h"

#include"Events/Events.h"
#include"Events/KeyEvents.h"
#include"Events/MouseEvents.h"
#include"3D/VertexArray.h"
#include"Platform/Windows/Console.h"
#include"Platform/API Core.h"
#include"Utilities/Memory.h"
#include"Display.h"
#include"3D/Texture.h"

#include"GUI/ImGuiLayer.h"
#include"GUI/CalculationLayer.h"

#include"3D/Framebuffer.h"

#include"Utilities/IO.h"
#include<filesystem>
#include "Utilities/PathUtils.h"

//#define IMSPINNER_DEMO
//#include "GUI/ImSpinner.h"


#pragma region Entry Point

int main(int argc, char** argv)
{
#ifdef PLATFORM_WIN
	SCESS::Console::SetHandle(GetStdHandle(STD_OUTPUT_HANDLE));
#endif // PLATFORM_WIN

	SCESS::PathUtils::Init();

	if (!std::filesystem::exists(SCESS::PathUtils::GetAssetPath().string()))
	{
		std::cout << "Required assets directory not found relative to executable.\n";
		system("PAUSE");
		return -1;
	}

	SCESS::Application smApp;
	return smApp.Run();
}


#pragma endregion

namespace SCESS
{
	//Constants
	constexpr float NEAR_PLANE = 0.1f;
	constexpr float FAR_PLANE = 100.0f;

	//Static Variable Initilisation
	Application* Application::s_Instance = nullptr;
	glm::mat4 Application::projectionMatrix(1.0f);


#define BIND_EVENT_FNC(x) std::bind(&Application::x, this, std::placeholders::_1)

	Application::Application() : m_MouseLastPos(0.0f)
	{
		std::cout << "Schrodinger's Cat's Electronic System Solver\n Meow!\n";
		UNI_ASSERT(s_Instance != nullptr, Application::Application, "There can only exist one instantiation of \"Application\"");
		s_Instance = this;
		m_Window.reset(Display::Create({ "Schrodinger's Cat's Electronic System Solver", 800, 600 }));
		m_Window->BindEventCallback(BIND_EVENT_FNC(OnEvent));

		glDisable(GL_CULL_FACE);
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);

		m_VersionStr = (const char*)glGetString(GL_VERSION);
		m_VersionStr.insert(0, "OpenGL ");
		m_VendorStr = (const char*)glGetString(GL_VENDOR);
		m_ShadingLanguageStr = (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION);
		m_ShadingLanguageStr.insert(0, "GLSL ");

		APICore::Initilise();

		m_ImGuiLayer = new SCESS::ImGuiLayer();
		PushOverlay(m_ImGuiLayer);

		mGraphLayer = new CalculationLayer();
		PushLayer(mGraphLayer);

		mCam.MoveTo(glm::vec3(0.0f, 0.0f, 0.0f), 10.0f);

		mainFramebuffer = Framebuffer::Create(FramebufferSpecification{ 800, 600, 1, {
			FramebufferAttachment(PixelFormat::RGBA8),
			FramebufferAttachment(PixelFormat::DEPTH24_STENCIL8)
			} });

		UpdateProjectionMatrix();

		TextureLoader loader(PathUtils::GetAssetPath(R"(Icon\SCESS Small.png)").string());
		mAppIcon = Texture2D::Create(loader);
	}

	Application::~Application()
	{
	}

	int Application::Run()
	{
		while (m_IsRunning)
		{
			{
				double newTime = m_Window->GetTime();
				m_DeltaTime = newTime - m_Time;
				m_Time = newTime;
				mCam.UpdateCamera(m_DeltaTime * 0.1f, m_KeyDown);
			}

			if (m_viewPortPanelSize.x != this->m_lastViewportSize.x ||
				m_viewPortPanelSize.y != this->m_lastViewportSize.y)
			{
				mainFramebuffer->Resize(m_viewPortPanelSize.x, m_viewPortPanelSize.y);

				for (Layer* layer : m_LayerStack)
					layer->Resize(m_viewPortPanelSize.x, m_viewPortPanelSize.y);
				//m_PlanetLayer->Resize(m_viewPortPanelSize.x, m_viewPortPanelSize.y);
				this->UpdateProjectionMatrix();
				m_lastViewportSize = m_viewPortPanelSize;
			}

			mainFramebuffer->Bind();
			glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

			//m_PlanetLayer->OnUpdate(mainFramebuffer, mCam, projectionMatrix, { NEAR_PLANE, FAR_PLANE });
			for (Layer* layer : m_LayerStack)
			{
				layer->OnUpdate(mainFramebuffer, mCam, projectionMatrix, { NEAR_PLANE, FAR_PLANE });
			}

			mainFramebuffer->Unbind();

			m_ImGuiLayer->Begin();
			RenderMenuBar();
			RenderSplashScreen();

			for (Layer* layer : m_LayerStack)
				layer->OnImGuiRender();

			RenderViewport();
			RenderCameraPanel();
			RenderSceneInformation();

			m_ImGuiLayer->End();

			m_Window->Refresh();
		}

		return 0;
	}

	void Application::RenderMenuBar()
	{
		if (ImGui::BeginMenuBar())
		{
			if (ImGui::BeginMenu("About"))
			{
				ImGui::MenuItem("Splash Screen", nullptr, &mSplashScreenOpen);
				ImGui::EndMenu();
			}
		}
		ImGui::EndMenuBar();
	}

	void Application::RenderSplashScreen()
	{
		ImVec2 center = ImGui::GetMainViewport()->GetCenter();
		ImGui::SetNextWindowPos(center, 0, ImVec2(0.5f, 0.5f));
		ImGui::SetNextWindowSize({ 600,600 });
		if (ImGui::BeginPopup("splashScreenPopup"))
		{
			constexpr float imgWidth = 350;
			ImGui::SetCursorPosX((ImGui::GetWindowWidth() - imgWidth) * 0.5f);
			ImGui::Image((ImTextureID)mAppIcon->GetRendererID(), { (float)imgWidth,
				(float)imgWidth });
			ImGui::PushFont(ImGui::GetIO().Fonts->Fonts[1]);
			ImGui::TextCentered("SCESS Version 1");
			ImGui::PopFont();
			ImGui::Text("Schrodinger's Cat's Electronic System Solver");
			ImGui::Text("By Charlie Wheelwright; cspw2");
			ImGui::Text("(Meow)");
			ImGui::EndPopup();
		}

		if (mSplashScreenOpen)
		{
			ImGui::OpenPopup("splashScreenPopup");
			mSplashScreenOpen = false;
		}
	}

	void Application::RenderViewport()
	{
		//ImGui::ShowDemoWindow(nullptr);
		ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2{ 0,0 });
		if (ImGui::Begin("3D Viewport", nullptr))
		{
			m_isViewportFocused = ImGui::IsWindowFocused();
			m_viewPortPanelSize = ImGui::GetContentRegionAvail();
			ImGui::Image(mainFramebuffer->QueryTextureAttachment()->GetRendererID(), ImVec2{ (float)mainFramebuffer->GetWidth(),
				(float)mainFramebuffer->GetHeight() });
		}
		ImGui::End();
		ImGui::PopStyleVar();
	}

	void Application::RenderCameraPanel()
	{
		if (!ImGui::Begin("Camera"))
		{
			ImGui::End();
			return;
		}

		if (ImGui::CollapsingHeader("Camera Properties", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Indent();
			static glm::vec3 newPosition{ 0.0f };
			static float newZoom = 0.0f;
			const glm::vec3& centrePos = mCam.GetCentrePosition();
			const glm::vec3& actualPos = mCam.GetActualPosition();

			newPosition = centrePos;
			newZoom = mCam.GetZoom();

			if (ImGui::Button("Reset Camera##cam"))
			{
				newPosition = { 0.0f, 0.0f, 0.0f };
				newZoom = 10.0f;
			}

			ImGui::Text("Camera Centre Position:");
			ImGui::DragFloat3("##cam_centre_pos", &newPosition[0]);
			ImGui::SameLine();
			if (ImGui::Button("Reset##cam_pos"))
				newPosition = { 0.0f, 0.0f, 0.0f };
			if (newPosition != centrePos)
				mCam.SetCentrePosition(newPosition);

			ImGui::Text("Camera Zoom:");
			ImGui::DragFloat("##cam_zoom", &newZoom, 0.5f);
			ImGui::SameLine();
			if (ImGui::Button("Reset##cam_zoom"))
				newZoom = 10.0f;
			if (newZoom != mCam.GetZoom())
				mCam.SetZoom(newZoom);

			ImGui::Text("Camera Actual Position: X=%.2f, Y=%.2f, Z=%.2f", actualPos.x, actualPos.y, actualPos.z);
			ImGui::Text("Pitch=%.2f, Yaw=%.2f, Roll=%.2f", mCam.GetPitch(), mCam.GetYaw(), mCam.GetRoll());
			ImGui::Unindent();
		}

		if (ImGui::CollapsingHeader("Render Style", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Indent();
			const char* items[] = { "Solid", "Wireframe" };
			static int selectedItem = 0;
			ImGui::Combo("Render Style##combo", &selectedItem, items, 2);

			for (Layer* layer : m_LayerStack)
				layer->SetRenderStyle(selectedItem == 0 ? GL_FILL : GL_LINE);

			ImGui::Unindent();
		}

		ImGui::End();
	}

	void Application::RenderSceneInformation()
	{
		if (!ImGui::Begin("Scene Information"))
		{
			ImGui::End();
			return;
		}

		if (ImGui::CollapsingHeader("Adaptor Information", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Text("Renderer API: %s", m_VersionStr.c_str());
			ImGui::Text("Vendor: %s", m_VendorStr.c_str());
			ImGui::Text("Shading Language Version: %s", m_ShadingLanguageStr.c_str());
		}

		if (ImGui::CollapsingHeader("Frame Information", ImGuiTreeNodeFlags_DefaultOpen))
		{
			const ImGuiIO& myIO = ImGui::GetIO();
			ImGui::Text("Frame rate: %.0f fps", myIO.Framerate);
			ImGui::Text("Frame time: %.1f ms", myIO.DeltaTime * 1000.0f);
		}

		ImGui::End();
	}

	void Application::PushLayer(Layer* layer)
	{
		m_LayerStack.PushLayer(layer);
	}

	void Application::PushOverlay(Layer* overlay)
	{
		m_LayerStack.PushOverlay(overlay);
	}

	void Application::UpdateProjectionMatrix()
	{
		projectionMatrix = glm::perspective<float>(glm::radians(45.0f),
			(float)mainFramebuffer->GetSpecification().Width / (float)mainFramebuffer->GetSpecification().Height, NEAR_PLANE, FAR_PLANE);
	}

	void Application::OnEvent(Event& e)
	{
		EventDispatcher dispatcher(e);

		auto DispatchApplicationEvents = [this, &dispatcher]() {
			if (dispatcher.GetEvent().IsHandled()) return;
			dispatcher.Dispatch<WindowCloseEvent>(BIND_EVENT_FNC(OnWindowClose));
			dispatcher.Dispatch<WindowResizeEvent>([this](WindowResizeEvent& e) -> bool
				{
					glViewport(0, 0, e.GetWidth(), e.GetHeight());
					this->UpdateProjectionMatrix();

					return true;
				});
			dispatcher.Dispatch<KeyPressedEvent>(BIND_EVENT_FNC(OnKeyPressed));
			dispatcher.Dispatch<KeyReleasedEvent>(BIND_EVENT_FNC(OnKeyReleased));
			dispatcher.Dispatch<MouseMoveEvent>(BIND_EVENT_FNC(OnMouseMove));
			dispatcher.Dispatch<MouseButtonPressedEvent>(BIND_EVENT_FNC(OnMouseButtonPressed));
			dispatcher.Dispatch<MouseButtonReleasedEvent>(BIND_EVENT_FNC(OnMouseButtonReleased));
			dispatcher.Dispatch<MouseScrollEvent>(BIND_EVENT_FNC(OnMouseScroll));
			};

		DispatchApplicationEvents();
		for (auto it = m_LayerStack.rbegin(); it != m_LayerStack.rend(); ++it)
		{
			if (e.IsHandled())
				break;
			(*it)->OnEvent(e);
		}
	}

	bool Application::OnWindowClose(WindowCloseEvent& e)
	{
		m_IsRunning = false;
		return true;
	}

	bool Application::OnKeyPressed(KeyPressedEvent& e)
	{
		m_KeyDown = e.GetKeyCode();
		m_KeyMods = e.GetModifiers();
		return false;
	}

	bool Application::OnKeyReleased(KeyReleasedEvent& e)
	{
		m_KeyDown = -1;
		m_KeyMods = 0;
		return false;
	}

	bool Application::OnMouseButtonPressed(MouseButtonPressedEvent& e)
	{
		m_MouseButtonDown = (int)e.GetMouseButton();
		return true;
	}

	bool Application::OnMouseButtonReleased(MouseButtonReleasedEvent& e)
	{
		m_MouseButtonDown = -1;
		return true;
	}

	bool Application::OnMouseMove(MouseMoveEvent& e)
	{
		const glm::dvec2& mousePos = e.GetMousePos();
		glm::dvec2 deltaMouse(mousePos.x - m_MouseLastPos.x, mousePos.y - m_MouseLastPos.y);
		if (m_isViewportFocused) mCam.UpdateMouse((float)deltaMouse.x, (float)deltaMouse.y, m_MouseButtonDown, m_KeyMods);
		m_MouseLastPos = mousePos;
		return true;
	}

	bool Application::OnMouseScroll(MouseScrollEvent& e)
	{
		double currentScrollTime = m_Window->GetTime();
		double deltaScroll = glm::max(1.0 / (currentScrollTime - mLastScrollTime), 1.0) * glm::sign(e.GetYOffset());
		mLastScrollTime = currentScrollTime;

		if (m_isViewportFocused) mCam.UpdateZoom(deltaScroll);
		return false;
	}

}