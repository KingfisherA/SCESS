#pragma once

#include "Quantum/Basis/Gaussian.h"
#include<vector>
#include<nlohmann/json.hpp>
#include<string_view>
#include<optional>

namespace SCESS::BasisUtils 
{

	std::vector<std::string> ListBasisSets();
	std::optional<nlohmann::json> LoadBasisSetJSON(std::string_view jsonFilePath);

	std::vector<Eigen::Vector3i> GenerateAngularMomenta(size_t angMom);
}