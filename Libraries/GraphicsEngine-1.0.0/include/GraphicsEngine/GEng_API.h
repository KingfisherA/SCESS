#pragma once

//|| _WIN64 || _W64
#ifdef _WIN32

#ifdef GENG_STATIC
#define GENG_API
#else
#ifdef GENG_BUILD
#define GENG_API __declspec(dllexport)
#else
#define GENG_API __declspec(dllimport)
#endif
#endif

#else
#define GENG_API
#endif