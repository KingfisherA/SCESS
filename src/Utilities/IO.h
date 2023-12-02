#pragma once
#include<string>
#include<string_view>

namespace SCESS::IO
{
	std::string ReadTextFile(const std::string_view filePath);

	namespace Dir
	{
		std::string_view GetDirectoryFromPath(const std::string_view path);
	}
}