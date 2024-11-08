#pragma once

#include "GraphicsEngine/GEng_API.h"
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <cmath>
#include <math.h>

GENG_API inline bool operator<(const glm::vec2 &v1, const glm::vec2 &v2)
{
	return glm::length2(v1) < glm::length2(v2);
}
GENG_API inline bool operator>(const glm::vec2 &v1, const glm::vec2 &v2)
{
	return glm::length2(v1) > glm::length2(v2);
}
GENG_API inline bool operator<=(const glm::vec2 &v1, const glm::vec2 &v2)
{
	return glm::length2(v1) <= glm::length2(v2);
}
GENG_API inline bool operator>=(const glm::vec2 &v1, const glm::vec2 &v2)
{
	return glm::length2(v1) >= glm::length2(v2);
}

/// <summary>
/// A namespace for performing mathematical opearations
/// </summary>
namespace GEng::Maths
{
	constexpr double PI = 3.141592653589793238;
	constexpr double Tau = PI * 2;
	constexpr double DegreesPerRadian = 180.0 / PI;
	constexpr double RadiansPerDegree = PI / 180.0;
	/// <summary>
	/// Interpolate between two values using a linear function
	/// </summary>
	/// <param name="a">The first value to interpolate between</param>
	/// <param name="b">The second value to interpolate between</param>
	/// <param name="f">The amount to interpolate</param>
	/// <returns>The interpolated value</returns>
	GENG_API float LinearInterpolate(float a, float b, float f);
	/// <summary>
	/// Interpolate between two values using a quicker, less accurate linear function
	/// </summary>
	/// <param name="a">The first value to interpolate between</param>
	/// <param name="b">The second value to interpolate between</param>
	/// <param name="f">The amount to interpolate</param>
	/// <returns>The interpolated value</returns>
	GENG_API float LinearInterpolateSimplified(float a, float b, float f);
	///< summary>
	/// Interpolate between two values using a cosine function
	/// </summary>
	/// <param name="y1">The first value to interpolate between</param>
	/// <param name="y2">The second value to interpolate between</param>
	/// <param name="mu">The amount to interpolate</param>
	/// <returns>The interpolated value</returns>
	GENG_API double CosineInterpolate(double y1, double y2, double mu);
	/// <summary>
	/// Rotate a 2D std::vector
	/// </summary>
	/// <param name="angleIn">The angle, in radians, to rotate</param>
	/// <param name="vectorIn">The std::vector to rotate</param>
	/// <returns>The rotated std::vector</returns>
	GENG_API glm::vec2 Rotate2DVector(float angleIn, glm::vec2 vectorIn);
	/// <summary>
	/// Create a spherical billboard rotation matrix
	/// </summary>
	/// <param name="rotationMatrix">The rotation matrix to edit with the new transformations</param>
	/// <param name="cameraPosition">The current position of the camera or focus target (where the object should be facing towards)</param>
	/// <param name="objPos">The object's position in world space</param>
	/// <returns>The spherically billboarded view matrix</returns>
	GENG_API glm::mat4 SphericalBillboard(const glm::vec3 &cameraPosition, const glm::vec3 &objectPosition);
	/// <summary>
	/// Convert a value in degrees to one in radians
	/// </summary>
	/// <param name="degrees">The value, in degrees, you wish to convert</param>
	/// <returns>The value in radians</returns>
	template <class A>
	inline A ToRadians(A degrees)
	{
		return degrees * RadiansPerDegree;
	}
	/// <summary>
	///  Convert a value in radians to one in degrees
	/// </summary>
	/// <param name="radians">The value, in radians, you wish to convert</param>
	/// <returns>The value in degrees</returns>
	template <class A>
	inline A ToDegrees(A radians)
	{
		return radians * DegreesPerRadian;
	}
	/// <summary>
	/// Remap a range to a new range
	/// </summary>
	/// <param name="number">The number to remap</param>
	/// <param name="oldRangeMin">The old minimum range</param>
	/// <param name="oldRangeMax">The old maximum range</param>
	/// <param name="newRangeMin">The new minimum range</param>
	/// <param name="newRangeMax">The new maximum range</param>
	/// <returns>The remapped value</returns>
	GENG_API float Remap(float number, float oldRangeMin, float oldRangeMax, float newRangeMin, float newRangeMax);
	/// <summary>
	/// Get an exponential interpolated distance value
	/// </summary>
	/// <param name="distanceIn">The distance to get the value for</param>
	/// <param name="maxDistance">The maximum distance that the value should be 0 past; the end of the curve</param>
	/// <param name="fadeDensity">The density of the fade</param>
	/// <param name="fadeGradient">The gradient of the fade</param>
	/// <returns>The interpolated value</returns>
	GENG_API float GetFadeValue(float distanceIn, float maxDistance, float fadeDensity, float fadeGradient);
	/// <summary>
	/// Get a linear interpolated distance value
	/// </summary>
	/// <param name="distanceIn">The distance to interpolate</param>
	/// <param name="minDistance">The value which defines the start of the linear function</param>
	/// <param name="fadeSpeed">The gradient of the linear function; how fast should the fade occur</param>
	/// <returns>The interpolated value</returns>
	GENG_API float GetLinearFadeValue(float distanceIn, float minDistance, float fadeSpeed);
	//
	template <typename T>
	T IntegerPow(T base, T exponent)
	{
		if (exponent == 0)
			return 1;
		T result = base;

		for (size_t i = 1; i < exponent; i++)
		{
			result *= base;
		}

		return result;
	}
	//
	GENG_API inline glm::mat3 CalculateNormalMatrix(const glm::mat4 &transformMatrix)
	{
		return glm::mat3(glm::transpose(glm::inverse(transformMatrix)));
	}
	//
}

namespace SV::Maths::Gravity
{
	constexpr double G = 6.67408E-11;
}
