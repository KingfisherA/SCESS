#pragma once

#include "../../GEng_API.h"
#include "GraphicsEngine/Platform/GLFW/GLFWDisplayBase.h"

#include <GLFW/glfw3.h>

namespace GEng {
    class GPUContext;

    GENG_API class OpenGLDisplay : public GLFWDisplayBase {
    public:
        OpenGLDisplay(GPUContext *ctx, const DisplayProps &properties);

        void Refresh() override;
    };
}
