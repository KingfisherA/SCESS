#include<pch.h>
#include"StringOperations.h"
//
//

namespace SCESS::StringOperations
{
	bool Replace(std::string& str, const std::string_view from, const std::string_view to)
	{
		size_t start_pos = str.find(from);
		if (start_pos == std::string::npos)
			return false;
		str.replace(start_pos, from.length(), to);
		return true;
	}

	std::vector<std::string> Split(const std::string_view s,
		const char seperator)
	{
		std::vector<std::string> output;
		std::string_view::size_type prev_pos = 0, pos = 0;

		while ((pos = s.find(seperator, pos)) != std::string::npos)
		{
			std::string substring(s.substr(prev_pos, pos - prev_pos));

			if (substring != "")
				output.push_back(substring);

			prev_pos = ++pos;
		}

		output.push_back(s.substr(prev_pos, pos - prev_pos).data()); // Last word
		return output;
	}

	std::vector<std::string> Split(const std::string_view s,
		const std::string_view seperator)
	{
		std::vector<std::string> output;
		std::string_view::size_type prev_pos = 0, pos = 0;

		while ((pos = s.find(seperator, pos)) != std::string::npos)
		{
			std::string substring(s.substr(prev_pos, pos - prev_pos));
			
			if (substring != "")
				output.push_back(substring);

			prev_pos = ++pos;
		}

		output.push_back(s.substr(prev_pos, pos - prev_pos).data()); // Last word
		return output;
	}

	void FindAll(const std::string_view str, const std::string_view searchPhrase, std::vector<size_t>& results)
	{
		size_t offset = 0;
		while (offset = str.find(searchPhrase, offset) != str.npos)
		{
			results.push_back(offset);
		}
	}
}