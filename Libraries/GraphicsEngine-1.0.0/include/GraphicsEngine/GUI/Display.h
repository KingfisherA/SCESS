#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Events/Events.h"
#include "GraphicsEngine/Platform/GPUContext.h"

#include <glm/vec2.hpp>

namespace GEng {
    GENG_API struct DisplayProps {
        std::string title;
        size_t width = 0, height = 0;
    };

    GENG_API class Display {
    public:
        using EventCallbackFnc = std::function<void(Event &)>;

        virtual ~Display() = default;

        virtual void Refresh() = 0;
        virtual size_t GetWidth() const = 0;
        virtual size_t GetHeight() const = 0;
        virtual glm::vec2 GetContentScaling() const =0;

        virtual void Show() const = 0;
        virtual void Hide() const = 0;
        virtual void Focus() const = 0;

        virtual bool CloseRequested() const = 0;
        virtual double GetTime() const = 0;

        virtual void BindEventCallback(const EventCallbackFnc &callback) = 0;

        virtual void *GetAPIDisplay() const = 0;

        static Display *Create(GPUContext *ctx, const DisplayProps &properties);
    };
}
