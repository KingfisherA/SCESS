#pragma once

#include "GraphicsEngine/GEng_API.h"
#include <vector>
#include <string>
#include <string_view>

namespace GEng::StringOperations
{
	GENG_API bool Replace(std::string &str, const std::string_view from, const std::string_view to);

	GENG_API std::vector<std::string> Split(const std::string_view s,
								   const char seperator);

	GENG_API std::vector<std::string> Split(const std::string_view s,
								   const std::string_view seperator);

	GENG_API void FindAll(const std::string_view str, const std::string_view searchPhrase, std::vector<size_t> &results);
}