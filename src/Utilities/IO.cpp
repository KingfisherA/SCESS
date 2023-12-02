#include<pch.h>
#include"IO.h"
#include"Utilities/Assert.h"
//
//

namespace SCESS::IO
{
	std::string ReadTextFile(const std::string_view filePath)
	{
		std::ifstream file(filePath.data());
		if (!file.good())
		{
			file.close();
			using namespace std::string_literals;
			UNI_FATAL(SCESS::IO::ReadFile, "Could not open file requested \""s + filePath.data() + "\"")
				return "";
		}
		std::string outputStr((std::istreambuf_iterator<char>(file)),
			std::istreambuf_iterator<char>());

		return outputStr;
	}
}

std::string_view SCESS::IO::Dir::GetDirectoryFromPath(const std::string_view path)
{
	size_t i = path.rfind('\\', path.length());
	if (i == path.npos)
		return "";

	return path.substr(0, i);
}
