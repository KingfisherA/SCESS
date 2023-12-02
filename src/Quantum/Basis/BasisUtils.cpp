#include<pch.h>

#include "BasisUtils.h"


//
//
std::optional<nlohmann::json> SCESS::BasisUtils::LoadBasisSetJSON(
	const std::string_view jsonFilePath)
{
	if (!std::filesystem::exists(jsonFilePath))
	{
		return std::nullopt;
	}

	std::ifstream jsonFile(jsonFilePath.data());
	
	nlohmann::json inputtedJson;
	try
	{
		inputtedJson = nlohmann::json::parse(jsonFile);
	}
	catch (const std::exception& e)
	{
		std::cout << "Could not parse JSON. \n" << e.what()<<'\n';
		return std::nullopt;
	}

	jsonFile.close();

	return inputtedJson;
}