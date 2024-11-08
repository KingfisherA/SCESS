#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/3D/Texture.h"

namespace GEng
{
	GENG_API class OpenGLTexture2D : public Texture2D
	{
	public:
		using byte = unsigned char;

		explicit OpenGLTexture2D(const Dimension2D &dimensions, const byte *textureData,
								 PixelFormat texFormat);
		explicit OpenGLTexture2D(GLenum texID, const Dimension2D &dimensions,
								 PixelFormat texFormat);
		OpenGLTexture2D(const OpenGLTexture2D &) = delete;
		OpenGLTexture2D(OpenGLTexture2D &&tex);

	private:
		void CreateTexture(const Dimension2D &dimensions, const byte *textureData,
						   PixelFormat texFormat);

	public:
		virtual ~OpenGLTexture2D();

		virtual void Bind(unsigned short unit) override;
		virtual void Unbind() override;

		virtual void GenerateMipMaps(const MipmapFiltering filtering) override;
		virtual void DeleteMipMaps() override;

		virtual const Dimension2D &GetDimensions() const override { return m_imageDimensions; }

		virtual void UpdateData(const Dimension2D offset, const Dimension2D size, const byte *data) override;

		virtual unsigned int GetRendererID() const noexcept override { return m_texID; }

	private:
		Dimension2D m_imageDimensions;

		GLenum m_texID = 0;
		PixelFormat m_pixelFormat;
		bool m_mipMapsGenerated = false;

		unsigned short m_texUnit = 0;
	};
}