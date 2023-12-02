#include<pch.h>
#include"OpenGLFramebuffer.h"
#include"Platform/OpenGL/OpenGLCore.h"
#include"Utilities/Assert.h"
//
//

namespace SCESS
{
	OpenGLFramebuffer::OpenGLFramebuffer(const FramebufferSpecification& spec) :
		specification(spec)
	{
		this->CreateFramebuffer();
	}
	void OpenGLFramebuffer::CreateFramebuffer()
	{
		auto CreateTex2D = [](size_t width, size_t height, PixelFormat format) -> unsigned int
		{
			GLuint tex;
			glGenTextures(1, &tex);

			glBindTexture(GL_TEXTURE_2D, tex);
			bool isDepth = PixelFormat::IsDepth(format);
			glTexImage2D(GL_TEXTURE_2D, 0, PixelFormat::FormatToGL(format), width, height, 0,
				isDepth ? GL_DEPTH_STENCIL : GL_RGB, isDepth ? GL_UNSIGNED_INT_24_8 : GL_UNSIGNED_BYTE, nullptr);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glBindTexture(GL_TEXTURE_2D, 0);

			return tex;
		};
		auto CreateRenderbuffer = [](size_t width, size_t height, PixelFormat format) -> unsigned int
		{
			GLuint rboID;
			glGenRenderbuffers(1, &rboID);

			glBindRenderbuffer(GL_RENDERBUFFER, rboID);
			glRenderbufferStorage(GL_RENDERBUFFER, PixelFormat::FormatToGL(format), width, height);
			glBindRenderbuffer(GL_RENDERBUFFER, 0);

			return rboID;
		};
		auto CreateMultisampledTex2D = [](size_t width, size_t height, PixelFormat format,
			unsigned int samples) -> unsigned int
		{
			GLuint tex;
			glGenTextures(1, &tex);

			glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, tex);
			glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, samples, PixelFormat::FormatToGL(format), width, height, GL_TRUE);
			glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, 0);

			return tex;
		};
		auto CreateMultisampledRenderbuffer = [](size_t width, size_t height, PixelFormat format,
			unsigned int samples) -> unsigned int
		{
			GLuint rboID;
			glGenRenderbuffers(1, &rboID);

			glBindRenderbuffer(GL_RENDERBUFFER, rboID);
			glRenderbufferStorageMultisample(GL_RENDERBUFFER, samples, PixelFormat::FormatToGL(format), width, height);
			glBindRenderbuffer(GL_RENDERBUFFER, 0);

			return rboID;
		};

		if (m_FBOID != 0)
		{
			ClearAttachments();
			glDeleteFramebuffers(1, &m_FBOID);
		}

		glGenFramebuffers(1, &m_FBOID);
		glBindFramebuffer(GL_FRAMEBUFFER, m_FBOID);

		bool depthAttached = false;
		unsigned int colourIndex = 0;
		for (size_t i = 0; i < specification.Attachments.Attachments.size(); i++)
		{
			FramebufferAttachment currentAttach = specification.Attachments.Attachments[i];
			unsigned int fboAttachment = 0;
			if (PixelFormat::IsDepth(currentAttach.PixelFormat))
			{
				fboAttachment = GL_DEPTH_STENCIL_ATTACHMENT;
			}
			else
			{
				fboAttachment = GL_COLOR_ATTACHMENT0 + colourIndex;
				colourIndex++;
			}
			UNI_ASSERT(depthAttached && PixelFormat::IsDepth(currentAttach.PixelFormat), SCESS::OpenGLFramebuffer::CreateFramebuffer, "Cannot have more than one depth attachment to an OpenGL framebuffer.");

			if (specification.Samples > 1)
			{
				if (currentAttach.Renderbuffer)
				{
					GLuint rboID = CreateMultisampledRenderbuffer(specification.Width, specification.Height, currentAttach.PixelFormat,
						specification.Samples);

					glFramebufferRenderbuffer(GL_FRAMEBUFFER, fboAttachment, GL_RENDERBUFFER, rboID);
					renderbufferAttachments.insert({ (unsigned int)i, rboID });
				}
				else
				{
					GLuint texID = CreateMultisampledTex2D(specification.Width, specification.Height, currentAttach.PixelFormat, specification.Samples);

					glFramebufferTexture2D(GL_FRAMEBUFFER, fboAttachment, GL_TEXTURE_2D_MULTISAMPLE, texID, 0);
					//textureAttachments.insert({ (unsigned int)i,  texID });
					textureAttachments.emplace(std::piecewise_construct,
						std::forward_as_tuple((unsigned int)i), std::forward_as_tuple(texID, 
							Dimension2D{specification.Width, specification.Height}, currentAttach.PixelFormat));
				}
			}
			else
			{
				if (currentAttach.Renderbuffer)
				{
					GLuint rboID = CreateRenderbuffer(specification.Width, specification.Height, currentAttach.PixelFormat);

					glFramebufferRenderbuffer(GL_FRAMEBUFFER, fboAttachment, GL_RENDERBUFFER, rboID);
					renderbufferAttachments.insert({ (unsigned int)i, rboID });
				}
				else
				{
					GLuint texID = CreateTex2D(specification.Width, specification.Height, currentAttach.PixelFormat);

					glFramebufferTexture2D(GL_FRAMEBUFFER, fboAttachment, GL_TEXTURE_2D, texID, 0);
					//textureAttachments.insert({ (unsigned int)i, texID });
					textureAttachments.emplace(std::piecewise_construct,
						std::forward_as_tuple((unsigned int)i), std::forward_as_tuple(texID,
							Dimension2D{ specification.Width, specification.Height }, currentAttach.PixelFormat));
				}
			}

			depthAttached |= PixelFormat::IsDepth(currentAttach.PixelFormat);
		}

		GLenum error = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
		if (error != GL_FRAMEBUFFER_COMPLETE)
		{
			ClearAttachments();
			glDeleteFramebuffers(1, &m_FBOID);
			m_FBOID = 0;

			auto FramebufferErrorToStr = [](GLenum errorCode) -> const char*
			{
				switch (errorCode)
				{
				case GL_FRAMEBUFFER_UNDEFINED:
					return "Default framebuffer is not available.";
				case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
					return "An attachment point to the framebuffer is incomplete";
				case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
					return "The framebuffer has an inadequate amount of attachments to be functional. Ensure at least one attachment is added to the framebuffer before creating it.";
				case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER:
					return "Framebuffer violation; A GL_FRAMEBUFFER_ATTACHMENT_OBJECT_TYPE is GL_NONE for a color attachment point(s) named by GL_DRAW_BUFFERi";
				case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER:
					return "GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER";
				case GL_FRAMEBUFFER_UNSUPPORTED:
					return "Combination of attachment internal formats is disallowed on the current OpenGL driver.";
				case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE:
					return "The number of samples for each multisampled attachment is not the same. Please ensure that each multisampled attachment used the same number of samples.";
				case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS:
					return "GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS";
				default:
					return "Unknown framebuffer error occured";
				}
			};


			UNI_FATAL(SCESS::OpenGLFramebuffer::CreateFramebuffer, "Framebuffer specification is invalid and"
				"hence framebuffer creation has failed.Change specification and recreate framebuffer. Error Message: "
				+ std::string(FramebufferErrorToStr(error)));
			return;
		}
	}
	//
	OpenGLFramebuffer::~OpenGLFramebuffer()
	{
		ClearAttachments();

		glDeleteFramebuffers(1, &m_FBOID);
		m_FBOID = 0;
	}
	//
	void OpenGLFramebuffer::Bind()
	{
		glBindFramebuffer(GL_FRAMEBUFFER, m_FBOID);
		glViewport(0, 0, specification.Width, specification.Height);
	}
	void OpenGLFramebuffer::Unbind()
	{
		glBindFramebuffer(GL_FRAMEBUFFER, 0);
	}
	//
	void OpenGLFramebuffer::Resize(size_t newWidth, size_t newHeight)
	{
		specification.Width = newWidth;
		specification.Height = newHeight;

		this->CreateFramebuffer();
	}
	//
	FramebufferAttachment* OpenGLFramebuffer::QueryAttachment(size_t index)
	{
		if (index + 1 > specification.Attachments.Attachments.size())
			return nullptr;

		return &specification.Attachments.Attachments[index];
	}
	Texture2D* OpenGLFramebuffer::QueryTextureAttachment(size_t index)
	{
		if (auto texIt = textureAttachments.find(index); texIt != textureAttachments.end())
		{
			return (Texture2D*)&(*texIt).second;
		}

		return nullptr;
	}
	void OpenGLFramebuffer::ClearAttachments()
	{
		for (auto& kvPair : renderbufferAttachments)
		{
			glDeleteRenderbuffers(1, &kvPair.second);
		}

		textureAttachments.clear();
		renderbufferAttachments.clear();
	}
	//

}