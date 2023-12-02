#include<pch.h>
#include"Framebuffer.h"
#include"Platform/API Core.h"
#include"Platform/OpenGL/OpenGLFramebuffer.h"
//
//

namespace SCESS
{
	sRef<Framebuffer> Framebuffer::Create(const FramebufferSpecification& spec)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			return make<OpenGLFramebuffer>(spec);
			break;
		default:
			return nullptr;
		}
	}
}
