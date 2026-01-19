#include<Quantum/Q_PCH.h>
#include "Viewer/Application/Application.h"
#include <ImGui/imgui.h>
#include <glm/gtc/matrix_transform.hpp>

#include "Dahlia/Events/MouseEvents.h"
#include "Dahlia/Platform/GPUContext.h"
#include "Dahlia/GUI/Display.h"
#include "Dahlia/3D/Texture.h"
#include "Dahlia/3D/Framebuffer.h"
#include "Dahlia/GUI/ImGuiLayer.h"
#include "Dahlia/Utilities/Assert.h"

#include "Viewer/GUI/CalculationLayer.h"
#include <filesystem>

#include "xc.h"
#include "Quantum/CC/CC.h"
#include "Quantum/Integrals/IntegralsLibInt.h"
#include "Resources/AssetFolder.h"

#pragma region Entry Point

int main(int argc, char **argv) {
#ifdef PLATFORM_WIN
	SCESS::Console::SetHandle(GetStdHandle(STD_OUTPUT_HANDLE));
#endif // PLATFORM_WIN
    std::cout << "Schrodinger's Cat's Electronic System Solver\n    By Charlie Wheelwright (cspw2) \n    Meow!\n";

    std::cout << "---Libraries---\n";
    int vmajor, vminor, vmicro;
    xc_version(&vmajor, &vminor, &vmicro);
    printf("    Using LibXC Version: %i.%i.%i\n", vmajor, vminor, vmicro);

    using namespace SCESS;

    Viewer::Application app;

    //Initilise libint2
    Integrals::LibInt::InitLibInt();

    if (!Res::AssetsFolder::Get().Init() || !app.Init(DH::GraphicsAPI::OpenGL, 4, 0)) {
        std::cout << "Failed to initialise SCESS application. Program will now close. \n";
        std::cout << "Press any key to continue...\n";

        //Close connection to libint2
        Integrals::LibInt::CloseLibInt();

        std::cin.get();
        return 1;
    }

    app.Run();

    //Close connection to libint2
    Integrals::LibInt::CloseLibInt();

    return 0;
}

#pragma endregion

namespace SCESS::Viewer {
    using namespace DH;
    // Constants
    constexpr float NEAR_PLANE = 0.1f;
    constexpr float FAR_PLANE = 100.0f;

    // Static Variable Initilisation
    Application *Application::sInstance = nullptr;

#define BIND_EVENT_FNC(x) std::bind(&Application::x, this, std::placeholders::_1)

    Application::Application() : GraphicsApp("Schrodinger's Cat's Electronic System Solver",
                                             800, 600,
                                             FramebufferSpecification{
                                                 800, 600, 1,
                                                 {
                                                     FramebufferAttachment(PixelFormat::RGBA8),
                                                     FramebufferAttachment(PixelFormat::DEPTH24_STENCIL8)
                                                 }
                                             }),
                                 mProjectionMatrix(1.0f) {
        DH_ASSERT(sInstance != nullptr, Application::Application,
                    "There can only exist one instantiation of \"Application\"");
        sInstance = this;
    }

    bool Application::OnInit() {
        this->mRenderViewport = true;
        this->mViewportName = "Molecule View";

        float dpiFactor = mWindow->GetContentScaling().x * 0.1f;
        // mWindow->SetSize(size_t((float) mWindow->GetWidth() * dpiFactor),
        //                  size_t((float) mWindow->GetHeight() * dpiFactor));
        mWindow->SetSize(800, 600);

        GPUContext *ctx = GPUContext::Get();
        ctx->SetDepthTest(false);
        ctx->SetFaceCulling(false);
        ctx->SetBlending(true, GPUContext::Blending::Alpha);

        mCalcLayer = new CalculationLayer();
        this->PushLayer((Layer*)mCalcLayer);

        mCam.SetCentrePosition({0.0f, 0.0f, 0.0f});
        mCam.SetZoom(10.0f);
        UpdateProjectionMatrix();

        TextureLoader loader(Res::AssetsFolder::Get().GetResourcePath("Icon/SCESS Small.png").string());
        mAppIcon = Texture2D::Create(loader);

        return true;
    }

    void Application::OnDestroy() {
        mAppIcon.reset();
    }


    void Application::OnViewportResize() {
        GraphicsApp::OnViewportResize();
        UpdateProjectionMatrix();
    }

    void Application::OnRender() {
        mCam.UpdateCamera((float) mDeltaTime * 0.1f);

        mViewportFramebuffer->Bind();
        mViewportFramebuffer->Clear({glm::vec3(0.01f), 1.0f});

        for (Layer *layer: mLayerStack)
            layer->OnUpdate(mViewportFramebuffer, mCam, mProjectionMatrix, {NEAR_PLANE, FAR_PLANE});
    }

    void Application::OnImGuiRender() {
        RenderSplashScreen();
        RenderMenuBar();
        RenderCameraPanel();
        RenderSceneInformation();
    }

    void Application::RenderMenuBar() {
        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("About")) {
                ImGui::MenuItem("Splash Screen", nullptr, &mSplashScreenOpen);
                ImGui::EndMenu();
            }
        }
        ImGui::EndMenuBar();
    }

    void Application::RenderSplashScreen() {
        ImVec2 center = ImGui::GetMainViewport()->GetCenter();
        ImGui::SetNextWindowPos(center, 0, ImVec2(0.5f, 0.5f));

        float dpiFactor = mWindow->GetContentScaling().x * 0.5f;
        ImGui::SetNextWindowSize({600.0f * dpiFactor, 600.0f * dpiFactor});
        if (ImGui::BeginPopup("splashScreenPopup")) {
            float imgWidth = 350.0f * dpiFactor;
            ImGui::SetCursorPosX((ImGui::GetWindowWidth() - imgWidth) * 0.5f);
            ImGui::Image((ImTextureID) (intptr_t) mAppIcon->GetRendererID(), {
                             (float) imgWidth,
                             (float) imgWidth
                         });
            ImGui::PushFont(ImGui::GetIO().Fonts->Fonts[1]);
            ImGui::TextCentered("SCESS Version 1");
            ImGui::PopFont();
            ImGui::Text("Schrodinger's Cat's Electronic System Solver");
            ImGui::Text("By Charlie Wheelwright; cspw2");
            ImGui::Text("(Meow)");
            ImGui::EndPopup();
        }

        if (mSplashScreenOpen) {
            ImGui::OpenPopup("splashScreenPopup");
            mSplashScreenOpen = false;
        }
    }

    void Application::RenderCameraPanel() {
        if (!ImGui::Begin("Camera")) {
            ImGui::End();
            return;
        }

        if (ImGui::CollapsingHeader("Camera Properties", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Indent();
            static glm::vec3 newPosition{0.0f};
            static float newZoom = 0.0f;
            const glm::vec3 &centrePos = mCam.GetCentrePosition();
            const glm::vec3 &actualPos = mCam.GetActualPosition();

            newPosition = centrePos;
            newZoom = mCam.GetZoom();

            if (ImGui::Button("Reset Camera##cam")) {
                newPosition = {0.0f, 0.0f, 0.0f};
                newZoom = 10.0f;
            }

            ImGui::Text("Camera Centre Position:");
            ImGui::DragFloat3("##cam_centre_pos", &newPosition[0]);
            ImGui::SameLine();
            if (ImGui::Button("Reset##cam_pos"))
                newPosition = {0.0f, 0.0f, 0.0f};
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

        if (ImGui::CollapsingHeader("Render Style", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Indent();
            const char *items[] = {"Solid", "Wireframe"};
            static int selectedItem = 0;
            ImGui::Combo("Render Style##combo", &selectedItem, items, 2);

            for (Layer *layer: mLayerStack)
                layer->SetRenderStyle((GPUContext::PolygonMode) selectedItem);

            ImGui::Unindent();
        }

        ImGui::End();
    }

    void Application::RenderSceneInformation() {
        if (!ImGui::Begin("Scene Information")) {
            ImGui::End();
            return;
        }

        if (ImGui::CollapsingHeader("Adaptor Information", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Text("Renderer API: %s", mVersionStr.c_str());
            ImGui::Text("Vendor: %s", mVendorStr.c_str());
            ImGui::Text("Shading Language Version: %s", mShadingLanguageStr.c_str());
        }

        if (ImGui::CollapsingHeader("Frame Information", ImGuiTreeNodeFlags_DefaultOpen)) {
            const ImGuiIO &myIO = ImGui::GetIO();
            ImGui::Text("Frame rate: %.0f fps", myIO.Framerate);
            ImGui::Text("Frame time: %.1f ms", myIO.DeltaTime * 1000.0f);
        }

        ImGui::End();
    }


    void Application::UpdateProjectionMatrix() {
        mProjectionMatrix = glm::perspective<float>(glm::radians(45.0f),
                                                    (float) mViewportFramebuffer->GetSpecification().Width / (float)
                                                    mViewportFramebuffer->GetSpecification().Height, NEAR_PLANE,
                                                    FAR_PLANE);
    }


    bool Application::OnMouseMove(MouseMoveEvent &e) {
        GraphicsApp::OnMouseMove(e);
        if (mIsViewportFocused) {
            mDeltaMousePos /= glm::vec2(mWindow->GetWidth(), mWindow->GetHeight());
            mCam.UpdateMovement(mDeltaMousePos, mMouseButtonDown, mKeyDown,
                                mKeyMods);
        }
        return true;
    }

    bool Application::OnMouseScroll(MouseScrollEvent &e) {
        GraphicsApp::OnMouseScroll(e);
        // if (mIsViewportFocused)
        //     mCam.UpdateZoom(mDeltaScroll);
        return true;
    }
}
