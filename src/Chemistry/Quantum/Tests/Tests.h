#pragma once

#include<glm/vec3.hpp>
#include<numbers>
#include <cmath>

namespace SCESS
{
	float OrbitalPZ(float x, float y, float z, const glm::vec3& pos) {
		x = x - pos.x;
		y = y - pos.y;
		z = z - pos.z;
		float r = std::sqrt(x * x + y * y + z * z);
		return z * std::exp(-r / 2) / std::sqrt(32 * std::numbers::pi);
	}

	float OrbitalDZ2(float x, float y, float z)
	{
		float r = std::sqrt(x * x + y * y + z * z);
		return (1.0 / (81 * std::sqrt(6 * std::numbers::pi))) * (3 * z * z - r * r) * std::exp(-r / 3);
	}
}