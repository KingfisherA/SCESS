#pragma once
#include <string_view>
#include <nlohmann/json.hpp>
#include <boost/locale/conversion.hpp>
#include <optional>
#include <iostream>

#include "Application/Application.h"

namespace SCESS
{
	class PeriodicTable
	{
	public:
		static nlohmann::json sPeriodicTable;

		static void OpenQuery()
		{
			if (sPeriodicTable == nlohmann::json(nullptr))
			{
				std::ifstream stream(Application::Get().GetResourcePath("Data/Periodic Table Data.json").string());
				sPeriodicTable = nlohmann::json::parse(stream);
				stream.close();
			}
		}

		template <typename T>
		static std::optional<T> QueryValueSym(const std::string_view eleSymbol, const std::string_view name)
		{
			for (nlohmann::json &obj : sPeriodicTable)
			{
				if (obj["Symbol"].get<std::string>() == eleSymbol)
				{
					return obj[name].get<T>();
				}
			}

			return {};
		}

		template <typename T>
		static std::optional<T> QueryValueName(std::string eleName, const std::string_view name)
		{
			std::string elementName = boost::locale::to_title<char>(eleName);
			for (nlohmann::json &obj : sPeriodicTable)
			{
				if (obj["Name"].get<std::string>() == elementName)
				{
					return obj[name].get<T>();
				}
			}

			return {};
		}

		static nlohmann::json QueryElementSym(std::string eleSymbol)
		{
			for (nlohmann::json &obj : sPeriodicTable)
			{
				if (obj["Symbol"].get<std::string>() == eleSymbol)
				{
					return obj;
				}
			}

			return {};
		}

		static nlohmann::json QueryElementName(std::string eleName)
		{
			std::string elementName = boost::locale::to_title<char>(eleName);
			for (nlohmann::json &obj : sPeriodicTable)
			{
				if (obj["Name"].get<std::string>() == elementName)
				{
					return obj;
				}
			}

			return {};
		}

		static void CloseQuery()
		{
			sPeriodicTable = nlohmann::json();
		}
	};
}