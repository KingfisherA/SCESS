#pragma once
#include<pch.h>

#include"Utilities/Memory.h"
#include<string_view>

namespace SCESS
{
	using Dimension2D = glm::vec<2, size_t>;

	struct PixelFormat
	{
	public:
		enum PixelEnum : unsigned char
		{
			None = 0,

			// Colour
			RGB,
			RGBA8,
			RG,
			RED_INTEGER,

			//Depth/stencil
			DEPTH24_STENCIL8
		};

		PixelFormat() noexcept = default;
		PixelFormat(PixelEnum pixEnum) noexcept : value(pixEnum) {}
		PixelFormat(const PixelFormat& format) noexcept : value(format.value) {}
		PixelFormat(PixelFormat&& format) noexcept : value(format.value) {}

		inline operator PixelEnum() noexcept { return value; }
		inline PixelFormat& operator=(const PixelFormat& format) noexcept
		{
			value = format.value;
			return *this;
		}
		inline PixelFormat& operator=(PixelFormat&& format) noexcept
		{
			value = format.value;
			return *this;
		}

		static bool IsDepth(PixelFormat format) noexcept
		{
			switch (format)
			{
			case SCESS::PixelFormat::None:
				return false;

			case SCESS::PixelFormat::RGBA8:
				return false;

			case SCESS::PixelFormat::RED_INTEGER:
				return false;

			case SCESS::PixelFormat::DEPTH24_STENCIL8:
				return true;

			}
		}
		bool IsDepth() { return IsDepth(*this); }
		static GLenum FormatToGL(PixelFormat format) noexcept
		{
			switch (format)
			{
			case SCESS::PixelFormat::None:
				return GL_NONE;
				
			case SCESS::PixelFormat::RGB:
				return GL_RGB;

			case SCESS::PixelFormat::RGBA8:
				return GL_RGB8;
				
			case SCESS::PixelFormat::RG:
				return GL_RG;

			case SCESS::PixelFormat::RED_INTEGER:
				return GL_RED;
				 
			case SCESS::PixelFormat::DEPTH24_STENCIL8:
				return GL_DEPTH24_STENCIL8;
				 
			}
		}
		GLenum FormatToGL() { return FormatToGL(*this); }

		PixelEnum value = PixelEnum::None;
	};

	class TextureLoader
	{
	public:
		using byte = unsigned char;

		explicit TextureLoader(const std::string_view filePath);
		~TextureLoader();

		byte* GetData() { return m_Data; }
		const byte* GetData() const { return m_Data; }

		const Dimension2D& GetDimensions() const { return m_imageDimensions; }

		const PixelFormat GetFormat() const { return m_imageFormat; }

	private:
		byte* m_Data = nullptr;
		Dimension2D m_imageDimensions;
		PixelFormat m_imageFormat;
	};

	enum class MipmapFiltering : unsigned char
	{
		Bilinear,
		Trilinear
	};

	class Texture
	{
	public:
		virtual ~Texture() = default;

		virtual void Bind(unsigned short unit) = 0;
		virtual void Unbind() = 0;

		virtual void* GetRendererID() const noexcept = 0;
	};

	class Texture2D : public Texture
	{
	public:
		using byte = unsigned char;

		virtual void GenerateMipMaps(const MipmapFiltering filtering) = 0;
		virtual void DeleteMipMaps() = 0;

		virtual const Dimension2D& GetDimensions() const = 0;
		const size_t GetWidth() const { return this->GetDimensions().x; }
		const size_t GetHeight() const { return this->GetDimensions().y; }
		const double GetAspect() const { return (double)this->GetDimensions().x / (double)this->GetDimensions().y; }
		const double GetInvAspect() const { return (double)this->GetDimensions().y / (double)this->GetDimensions().x; }

		virtual void UpdateData(const Dimension2D offset, const Dimension2D size, const byte* data) = 0;

		static sRef<Texture2D> Create(const std::string_view filePath);
		static sRef<Texture2D> Create(const TextureLoader& imageData);
		static sRef<Texture2D> Create(const Dimension2D& dimensions, const byte* data, const PixelFormat format);
	};

}