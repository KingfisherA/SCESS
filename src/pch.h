#pragma once

#include<iostream>
#include<format>
#include<fstream>
#include<functional>
#include<thread>
#include<algorithm>
#include<utility>
#include<exception>
#include<fstream>
#include<filesystem>
#include<functional>

#include<string>
#include<sstream>
#include<string_view>
#include<vector>
#include<unordered_map>
#include<unordered_set>
#include<numbers>

//Open-GL Includes
#include<GL/glew.h>
#include<GLFW/glfw3.h>

#define GLM_RADIANS
#include<glm/glm.hpp>
#include<glm/gtx/transform.hpp>
#include<glm/gtx/norm.hpp>
#include <glm/gtc/type_ptr.hpp>

using GLid = unsigned int;
using Index = size_t;
constexpr size_t IndexUndefined = 0xffffffffffffffff;

#ifdef PLATFORM_WIN
#define NOMINMAX
	#include<Windows.h>
#endif // PLATFORM_WIN

//Eigen includes
#include<Eigen/Core>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>

//JSON parsing support
#include<nlohmann/json.hpp>

//Boost
#include "boost/filesystem.hpp"
