#pragma once
#include<pch.h>


class OpenGLCore abstract
{
public: 
	enum class Blending
	{
		Alpha, Additive
	};

	static GLenum Initilise(void (*ErrorCallback)(int errorCode, const char* message));
	static void SetBlending(const Blending blend);

};

