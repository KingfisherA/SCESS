#include<pch.h>
#include "Display.h"
#include "Platform/GLFW/GLFWDisplay.h"
#include "Utilities/Memory.h"
//
//

namespace SCESS
{
	//
	Display *Display::Create(const DisplayProps &properties)
	{
		return new GLFWDisplay(properties);
	}
	//
}
