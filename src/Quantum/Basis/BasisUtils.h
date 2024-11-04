#pragma once

#include "Quantum/Basis/Gaussian.h"
#include<vector>
#include<nlohmann/json.hpp>
#include<string_view>
#include<optional>

namespace SCESS::BasisUtils 
{
	struct BasisSet
	{
	public:

	};

	std::vector<std::string> ListBasisSets();
	std::optional<nlohmann::json> LoadBasisSetJSON(const std::string_view jsonFilePath);

	std::vector<Eigen::Vector3i> GenerateAngularMomenta(const size_t angMom);
}