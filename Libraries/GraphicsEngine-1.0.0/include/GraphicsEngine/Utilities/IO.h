#pragma once

#include"GraphicsEngine/GEng_API.h"
#include <string>
#include <string_view>

namespace GEng::IO
{
	GENG_API std::string ReadTextFile(const std::string_view filePath);

	namespace Dir
	{
		GENG_API std::string_view GetDirectoryFromPath(const std::string_view path);
	}
}