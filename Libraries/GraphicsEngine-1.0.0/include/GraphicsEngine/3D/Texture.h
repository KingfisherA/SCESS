#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Utilities/Memory.h"
#include <string_view>
#include <glm/glm.hpp>

#include <GL/glew.h>

namespace GEng {
    using Dimension2D = glm::vec<2, size_t>;

    GENG_API struct PixelFormat {
    public:
        enum PixelEnum : unsigned char {
            None = 0,

            // Colour
            R,
            RGB,
            RGBA8,
            RG,
            RED_INTEGER,

            // Depth/stencil
            DEPTH24_STENCIL8
        };

        PixelFormat() noexcept = default;

        PixelFormat(PixelEnum pixEnum) noexcept : mValue(pixEnum) {
        }

        PixelFormat(const PixelFormat &format) noexcept : mValue(format.mValue) {
        }

        PixelFormat(PixelFormat &&format) noexcept : mValue(format.mValue) {
        }

        inline operator PixelEnum() noexcept { return mValue; }

        inline PixelFormat &operator=(const PixelFormat &format) noexcept {
            mValue = format.mValue;
            return *this;
        }

        inline PixelFormat &operator=(PixelFormat &&format) noexcept {
            mValue = format.mValue;
            return *this;
        }

        static bool IsDepth(PixelFormat format) noexcept {
            switch (format) {
                case GEng::PixelFormat::DEPTH24_STENCIL8:
                    return true;

                default:
                    return false;
            }
        }

        bool IsDepth() { return IsDepth(*this); }

        static GLenum FormatToGL(PixelFormat format) noexcept {
            switch (format) {
                case GEng::PixelFormat::None:
                    return GL_NONE;

                case GEng::PixelFormat::RGB:
                    return GL_RGB;

                case GEng::PixelFormat::RGBA8:
                    return GL_RGBA;

                case GEng::PixelFormat::RG:
                    return GL_RG;

                case GEng::PixelFormat::R:
                case GEng::PixelFormat::RED_INTEGER:
                    return GL_RED;

                case GEng::PixelFormat::DEPTH24_STENCIL8:
                    return GL_DEPTH24_STENCIL8;
            }
        }

        GLenum FormatToGL() { return FormatToGL(*this); }

        PixelEnum mValue = PixelEnum::None;
    };

    GENG_API class TextureLoader {
    public:
        using byte = unsigned char;

        explicit TextureLoader(const std::string_view filePath);

        ~TextureLoader();

        byte *GetData() { return m_Data; }
        const byte *GetData() const { return m_Data; }

        const Dimension2D &GetDimensions() const { return m_imageDimensions; }

        const PixelFormat GetFormat() const { return m_imageFormat; }

    private:
        byte *m_Data = nullptr;
        Dimension2D m_imageDimensions;
        PixelFormat m_imageFormat;
    };

    GENG_API enum class MipmapFiltering : unsigned char {
        Bilinear,
        Trilinear
    };

    GENG_API class Texture {
    public:
        virtual ~Texture() = default;

        virtual void Bind(unsigned short unit) = 0;

        virtual void Unbind() = 0;

        virtual unsigned int GetRendererID() const noexcept = 0;
    };

    GENG_API class Texture2D : public Texture {
    public:
        using byte = unsigned char;

        virtual void GenerateMipMaps(const MipmapFiltering filtering) = 0;

        virtual void DeleteMipMaps() = 0;

        virtual const Dimension2D &GetDimensions() const = 0;

        const size_t GetWidth() const { return this->GetDimensions().x; }
        const size_t GetHeight() const { return this->GetDimensions().y; }
        const double GetAspect() const { return (double) this->GetDimensions().x / (double) this->GetDimensions().y; }

        const double GetInvAspect() const {
            return (double) this->GetDimensions().y / (double) this->GetDimensions().x;
        }

        virtual void UpdateData(const Dimension2D offset, const Dimension2D size, const byte *data) = 0;

        static sRef<Texture2D> Create(const std::string_view filePath);

        static sRef<Texture2D> Create(const TextureLoader &imageData);

        static sRef<Texture2D> Create(const Dimension2D &dimensions, const byte *data, const PixelFormat format);
    };
}
