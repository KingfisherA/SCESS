#include<pch.h>

#ifdef PLATFORM_WIN
#include "Console.h"
//
//

void* SCESS::Console::cHandle = nullptr;

namespace SCESS
{
	//
	void Console::SetHandle(void* handle)
	{
		cHandle = handle;
	}
	//
	void Console::SetConsoleColour(ConsoleColour col)
	{
		SetConsoleTextAttribute(cHandle, WORD(col));
	}
	//
}
#endif
