#pragma once
#include<pch.h>

#ifdef PLATFORM_WIN
#include"Platform/Windows/Console.h"

#define UNI_ASSERT(condition, object, message) if(condition){SCESS::Console::SetConsoleColour(SCESS::ConsoleColour::Red); std::cerr <<"Assertion Failed! [" << #object << "] Message: " << message << '\n'; }
#define UNI_FATAL(object, message) SCESS::Console::SetConsoleColour(SCESS::ConsoleColour::BrightRed); std::cerr <<"[" << #object << "] Message: " << message << '\n';
#define UNI_WARN(object, message) SCESS::Console::SetConsoleColour(SCESS::ConsoleColour::BrightYellow); std::cerr <<"[" << #object << "] Message: " << message << '\n';

#define UNI_LOG(col, object, message) SCESS::Console::SetConsoleColour(SCESS::ConsoleColour::col); std::cerr <<'[' << #object << ']' << ' ' << message << '\n';

#endif