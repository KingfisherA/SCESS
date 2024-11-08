#pragma once

#include "GraphicsEngine/GEng_API.h"
#include <map>
#include "GraphicsEngine/3D/Framebuffer.h"
#include "GraphicsEngine/Platform/OpenGL/OpenGLFramebuffer.h"
#include "GraphicsEngine/Platform/OpenGL/OpenGLTexture.h"

namespace GEng {
    GENG_API class OpenGLFramebuffer : public Framebuffer {
    public:
        OpenGLFramebuffer(const FramebufferSpecification &spec);

    private:
        void CreateFramebuffer();

    public:
        virtual ~OpenGLFramebuffer();

        virtual void Bind() override;
        virtual void Unbind() override;

        virtual void Resize(size_t newWidth, size_t newHeight) override;

        virtual FramebufferAttachment *QueryAttachment(size_t index = 0) override;
        virtual Texture2D *QueryTextureAttachment(size_t index = 0) override;

        virtual const FramebufferSpecification &GetSpecification() const override { return specification; }
        virtual FramebufferSpecification &GetSpecification() override { return specification; }

        void Clear(const glm::vec4 &clearColour) override;

        virtual unsigned int GetRenderID() const { return m_FBOID; }

    private:
        void ClearAttachments();

        bool mHasDepth = false;

        unsigned int m_FBOID = 0;
        FramebufferSpecification specification;

        std::map<unsigned int, OpenGLTexture2D> textureAttachments;
        std::map<unsigned int, GLuint> renderbufferAttachments;
    };
}
