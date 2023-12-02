#pragma once
#include<pch.h>
#include"Platform/OpenGL/OpenGLCore.h"
#include"Utilities/Assert.h"

enum class API
{
	OpenGL = 0 
};

constexpr API PROGRAM_API = API::OpenGL;

class APICore abstract
{
public:
	static void Initilise()
	{
		using namespace std::literals::string_literals;
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			OpenGLCore::Initilise([](int err, const char* mes) {UNI_FATAL(OpenGLCore;  GLEW API Initilisation Error, "Unable to intilise GLEW. Error Message: "s + mes); });
			break;
		}
	}
};