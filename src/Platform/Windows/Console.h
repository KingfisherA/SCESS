#pragma once
#include<iostream>

#ifdef PLATFORM_WIN
namespace SCESS
{
	enum class ConsoleColour : unsigned char
	{
		Black = 0,
		Blue = 1,
		Green = 2,
		Aqua = 3,
		Red = 4,
		Purple = 5,
		Yellow = 6,
		White = 7,
		Grey = 8,
		BrightBlue = 9,
		BrightGreen = 10,
		BrightAqua = 11,
		BrightRed = 12,
		BrightPurple = 13,
		BrightYellow = 14
	};

	class Console
	{
	public:
		//
		static void SetHandle(void* handle);
		static void SetConsoleColour(ConsoleColour col);
		//
	private:
		static void* cHandle;
	};
}
#endif

