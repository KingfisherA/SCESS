#pragma once

#include <ImGui/imgui.h>

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Layers/Layer.h"
#include "GraphicsEngine/Events/Events.h"
#include "GraphicsEngine/GUI/Text.h"


namespace ImGui {
    GENG_API void TextCentered(const char *txt);
}

namespace GEng {
    class Display;
    class Font;

    GENG_API class ImGuiLayer : public Layer {
    public:
        class Context {
        public:
            // Context() {
            //     ctx = ImGui::CreateContext();
            //     ImGui::GetAllocatorFunctions(&allocFunc, &freeFunc, &userData);
            // }
            Context() = default;
            Context(const Context &) = default;
            ImGuiContext *ctx = nullptr;
            ImGuiMemAllocFunc allocFunc = nullptr;
            ImGuiMemFreeFunc freeFunc = nullptr;
            void *userData = nullptr;
        };

        ImGuiLayer(const Display *mainDisplay, const Font &font);
        ~ImGuiLayer() override = default;

        void OnInit() override;
        void OnDestroy() override;
        void OnImGuiRender() override;

        void OnEvent(Event &e) override;

        void RenderDockspace();


        [[nodiscard]] Context &GetContext() { return mContext; }

        void Begin();
        void End();

    private:
        void SetDarkThemeColors();

        Font mFont;
        Context mContext;

        const Display *mDisplay;
    };
}
