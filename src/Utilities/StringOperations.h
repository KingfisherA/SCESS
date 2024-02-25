#pragma once

#include<vector>
#include<string>
#include<string_view>

namespace SCESS::StringOperations
{
	bool Replace(std::string& str, const std::string_view from, const std::string_view to);

	std::vector<std::string> Split(const std::string_view s,
		const char seperator);

	std::vector<std::string> Split(const std::string_view s,
		const std::string_view seperator);

	void FindAll(const std::string_view str, const std::string_view searchPhrase, std::vector<size_t>& results);
}