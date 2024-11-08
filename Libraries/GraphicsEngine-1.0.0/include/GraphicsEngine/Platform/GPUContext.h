#pragma once

#include"GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Utilities/Assert.h"

#include<glm/vec4.hpp>
#include <string>


namespace GEng {
    GENG_API enum class API {
        OpenGL = 0
    };

    GENG_API class GPUContext {
    public:
        virtual ~GPUContext() = default;

        GENG_API enum class PolygonMode {
            Solid = 0,
            Lines,
            Points
        };

        GENG_API enum class Blending {
            Alpha = 0,
            Additive
        };

        static GPUContext *Create(const API api, uint32_t major, uint32_t minor);

        template<typename T = GPUContext>
        static T *Get() { return reinterpret_cast<T *>(sCurrentContext); }
        static API GetAPI() { return sCurrentContext->GetContextAPI(); }

        virtual API GetContextAPI() const = 0;
        virtual void Reset() =0;

        virtual void Clear(const glm::vec4 &clearColour) =0;

        virtual std::string GetVersionStr() =0;
        virtual std::string GetVendorStr() =0;
        virtual std::string GetShadingVersionStr() =0;

        virtual void SetFaceCulling(bool activate) =0;
        virtual void SetDepthTest(bool activate) =0;
        virtual void SetPolygonMode(const PolygonMode mode) =0;
        virtual void SetBlending(bool activate, const Blending mode = Blending::Alpha) =0;

        virtual void SetViewport(const glm::ivec4 &bounds) =0;

        virtual void ActivateDefaultFramebuffer() =0;

        virtual void DestroyContext() =0;

    private:
        static GPUContext *sCurrentContext;
    };
}
