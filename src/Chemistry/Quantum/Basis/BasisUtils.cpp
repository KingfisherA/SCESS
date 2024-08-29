#include<pch.h>

#include "Chemistry/Quantum/Basis/BasisUtils.h"
#include"Application/Application.h"

#include<boost/filesystem.hpp>


//
//
namespace SCESS::BasisUtils
{
	std::vector<std::string> ListBasisSets()
	{
		namespace fs = boost::filesystem;
		std::vector<std::string> basisSets;
		fs::path basisSetPath = Application::Get().GetResourcePath("Basis");
		
		fs::recursive_directory_iterator itr(basisSetPath);

		fs::recursive_directory_iterator end;

		while (itr != end)
		{
			if (fs::is_regular_file(itr->path()))
			{
				basisSets.emplace_back(itr->path().stem().string());
			}

			boost::system::error_code ec;
			itr.increment(ec);
		}

		return basisSets;
	}

	std::optional<nlohmann::json> LoadBasisSetJSON(
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
			std::cout << "Could not parse JSON. \n" << e.what() << '\n';
			return std::nullopt;
		}

		jsonFile.close();

		return inputtedJson;
	}

}
