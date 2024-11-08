#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "boost/filesystem.hpp"

namespace GEng
{
	GENG_API class PathUtils
	{
	public:
		static void Init()
		{
			boost::filesystem::initial_path();
		}

		static boost::filesystem::path GetAbsolute(const std::string_view relative)
		{
			boost::filesystem::path exeDir = boost::filesystem::initial_path();
			exeDir.append(relative);
			return exeDir;
		}

		static boost::filesystem::path GetAssetPath()
		{
#ifdef OS_LINUX
			return boost::filesystem::path("../../../Assets");
#elif OS_WIN
			return boost::filesystem::path("D:/dev/C++/SCESS/Assets");
#endif // OS
		}

		static boost::filesystem::path GetAssetPath(const std::string_view relativeAsset)
		{
			return GetAssetPath().append(relativeAsset);
		}
	};
}
