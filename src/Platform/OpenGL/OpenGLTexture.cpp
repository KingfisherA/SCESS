#include<pch.h>

#include"Platform/OpenGL/OpenGLCore.h"
#include"Platform/OpenGL/OpenGLTexture.h"

//
//

namespace SCESS
{
	OpenGLTexture2D::OpenGLTexture2D(const Dimension2D& dimensions, 
		const byte* textureData, PixelFormat texFormat) : 
		m_imageDimensions(dimensions),
		m_pixelFormat(texFormat)
	{
		CreateTexture(dimensions, textureData, texFormat);
	}

	OpenGLTexture2D::OpenGLTexture2D(GLenum texID, const Dimension2D& dimensions,
		PixelFormat texFormat) :
		m_texID(texID),
		m_imageDimensions(dimensions),
		m_pixelFormat(texFormat)
	{
	}

	OpenGLTexture2D::OpenGLTexture2D(OpenGLTexture2D&& tex) :
		m_imageDimensions(tex.m_imageDimensions),
		m_mipMapsGenerated(tex.m_mipMapsGenerated),
		m_pixelFormat(tex.m_pixelFormat),
		m_texID(tex.m_texID),
		m_texUnit(tex.m_texUnit)
	{
		tex.m_texID = 0;
	}

	void OpenGLTexture2D::CreateTexture(const Dimension2D& dimensions, const byte* textureData,
		PixelFormat texFormat)
	{
		if (m_texID != 0)
		{
			glDeleteTextures(1, &m_texID);
		}

		glGenTextures(1, &m_texID);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, m_texID);
		glTexImage2D(GL_TEXTURE_2D, 0, texFormat.FormatToGL(), dimensions.x, dimensions.y, 
			0, texFormat.FormatToGL(), GL_UNSIGNED_BYTE, (const void*)textureData);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	OpenGLTexture2D::~OpenGLTexture2D()
	{
		glDeleteTextures(1, &m_texID);
		m_texID = 0;
	}

	void OpenGLTexture2D::Bind(unsigned short unit)
	{
		m_texUnit = unit;

		glActiveTexture(GL_TEXTURE0 + unit);
		glBindTexture(GL_TEXTURE_2D, m_texID);
	}

	void OpenGLTexture2D::Unbind()
	{
		glActiveTexture(GL_TEXTURE0 + m_texUnit);
		glBindTexture(GL_TEXTURE_2D, 0);
	}

	void OpenGLTexture2D::GenerateMipMaps(const MipmapFiltering filtering)
	{
		Bind(0);
		if(!m_mipMapsGenerated) glGenerateMipmap(GL_TEXTURE_2D);
		glTextureParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTextureParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
			[filtering]()
			{
				switch (filtering)
				{
				default:
				case MipmapFiltering::Bilinear:
					return GL_LINEAR_MIPMAP_NEAREST;
				case MipmapFiltering::Trilinear:
					return GL_LINEAR_MIPMAP_LINEAR;
				}
			}());
		Unbind();

		m_mipMapsGenerated = true;
	}

	void OpenGLTexture2D::DeleteMipMaps()
	{
		Bind(0);
		glTextureParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTextureParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		Unbind();
	}

	void OpenGLTexture2D::UpdateData(const Dimension2D offset, const Dimension2D size, const byte* data)
	{
		Bind(0);
		glTexSubImage2D(GL_TEXTURE_2D, 0, offset.x, offset.y, size.x, size.y, m_pixelFormat.FormatToGL(), 
			GL_UNSIGNED_BYTE, (const void*)data);
		Unbind();
	}
}