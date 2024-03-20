#pragma once
#include "boost/filesystem.hpp"

namespace SCESS
{
	class PathUtils
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
			//return boost::filesystem::initial_path().append("\\Assets");
			return boost::filesystem::path("/media/charlie/Data/dev/C++/SCESS/Assets");
		}

		static boost::filesystem::path GetAssetPath(const std::string_view relativeAsset)
		{
			return GetAssetPath().append(relativeAsset);
		}
	};
}
