#include<pch.h>

#include"Platform/API Core.h"
#include"Utilities\Assert.h"

#include"Texture.h"
#include"Platform/OpenGL/OpenGLTexture.h"
#include"stb/stb_image.h"

//
//

namespace SCESS
{

#pragma region Texture Data

	TextureLoader::TextureLoader(const std::string_view filePath)
	{
		stbi_set_flip_vertically_on_load(0);
		int width, height, channels;
		m_Data = (byte*)stbi_load(filePath.data(), &width, &height, &channels, 0);

		UNI_ASSERT(m_Data == nullptr, SCESS::TextureLoader::TextureData, "Image file cannot be loaded into memory.");

		m_imageDimensions.x = width;
		m_imageDimensions.y = height;

		m_imageFormat = [channels]() {
			switch (channels)
			{
			default:
				return PixelFormat::None;
			case 1:
				return PixelFormat::RED_INTEGER;
			case 2:
				return PixelFormat::RG;
			case 3:
				return PixelFormat::RGB;
			case 4:
				return PixelFormat::RGBA8;
			}
		}();

	}

	TextureLoader::~TextureLoader()
	{
		free(m_Data);
		m_Data = nullptr;
	}

#pragma endregion

#pragma region Texture2D

	sRef<Texture2D> Texture2D::Create(const std::string_view filePath)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			{
				TextureLoader imageData(filePath);
				return make<OpenGLTexture2D>(imageData.GetDimensions(), imageData.GetData(), 
					imageData.GetFormat());
			}
		}
	}

	sRef<Texture2D> Texture2D::Create(const TextureLoader& imageData)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			{
				return make<OpenGLTexture2D>(imageData.GetDimensions(), imageData.GetData(),
					imageData.GetFormat());
			}
		}
	}

	sRef<Texture2D> Texture2D::Create(const Dimension2D& dimensions, const byte* data, const PixelFormat format)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			return make<OpenGLTexture2D> (dimensions, data, format);
		}
	}

#pragma endregion
}