#pragma once
#include <iostream>

#ifdef PLATFORM_WIN
#include "GraphicsEngine/Platform/Windows/Console.h"

#define GENG_ASSERT(condition, object, message)                                             \
    if (condition)                                                                         \
    {                                                                                      \
        GEng::Console::SetConsoleColour(GEng::ConsoleColour::Red);                         \
        std::cerr << "Assertion Failed! [" << #object << "]  " << message << '\n'; \
    }
#define GENG_FATAL(object, message)                                   \
    GEng::Console::SetConsoleColour(GEng::ConsoleColour::BrightRed); \
    std::cerr << "[" << #object << "]  " << message << '\n';
#define GENG_WARN(object, message)                                       \
    GEng::Console::SetConsoleColour(GEng::ConsoleColour::BrightYellow); \
    std::cerr << "[" << #object << "]  " << message << '\n';

#define GENG_LOG(col, object, message)                          \
    GEng::Console::SetConsoleColour(GEng::ConsoleColour::col); \
    std::cerr << '[' << #object << ']' << ' ' << message << '\n';

#else

#define GENG_ASSERT(condition, object, message)                                             \
    if (condition)                                                                         \
    {                                                                                      \
        std::cerr << "Assertion Failed! [" << #object << "] " << message << '\n'; \
        throw std::runtime_error("");\
    }
#define GENG_FATAL(object, message) std::cerr << '[' << #object << "] " << message << '\n';
#define GENG_WARN(object, message) std::cerr << '[' << #object << "] " << message << '\n';

#define GENG_LOG(col, object, message) std::cerr << '[' << #object << "] " << message << '\n';

#endif
