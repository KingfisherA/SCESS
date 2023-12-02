#include<pch.h>
#include "OpenGLCore.h"
#include"Utilities/Assert.h"
//
//

static std::unordered_set<GLuint> gPreviousErrorId;
GLenum OpenGLCore::Initilise(void(*ErrorCallback)(int errorCode, const char* message))
{
	GLenum glewError = glewInit();
	if (glewError != GLEW_OK)
	{
		ErrorCallback(glewError, (const char*)glewGetErrorString(glewError));
		return glewError;
	}

	glEnable(GL_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	
	glEnable(GL_DEBUG_OUTPUT);
	glDebugMessageCallback([](
		GLenum source, 
		GLenum type, 
		GLuint id, 
		GLenum serverity, 
		GLsizei length, 
		const GLchar* message, 
		const void* userParam) {
			if (gPreviousErrorId.find(id) == gPreviousErrorId.end())
			{
				gPreviousErrorId.insert(id);
				if (type == GL_DEBUG_TYPE_ERROR)
				{
					UNI_FATAL(OpenGL Error, message);
				}
				else
				{
					UNI_LOG(Blue, OpenGL Information, message);
				}
			}
		}, nullptr);

	UNI_LOG(Green, OpenGL GLEW API, "GLEW Intilised OK!");

	return glewError;
}
//
void OpenGLCore::SetBlending(const Blending blend)
{
	switch (blend)
	{
	case Blending::Alpha:
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		break;
	case Blending::Additive:
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		break;
	}
}
//