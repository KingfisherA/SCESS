#pragma once
#include"Utilities/Memory.h"
#include"Texture.h"

namespace SCESS
{

	struct FramebufferAttachment
	{
	public:
		FramebufferAttachment() = default;
		FramebufferAttachment(PixelFormat format) : PixelFormat(format) { }
		FramebufferAttachment(PixelFormat format, const bool isRenderbuffer) : 
			PixelFormat(format), Renderbuffer(isRenderbuffer) { }

		PixelFormat PixelFormat = PixelFormat::None;
		bool Renderbuffer = false;
	};


	struct FramebufferAttachmentSpecification
	{
		FramebufferAttachmentSpecification() = default;
		FramebufferAttachmentSpecification(const std::initializer_list<FramebufferAttachment>& attachments)
			: Attachments(attachments) {}

		std::vector<FramebufferAttachment> Attachments;
	};

	struct FramebufferSpecification
	{
	public:
		size_t Width = 0, Height = 0, Samples = 1;
		FramebufferAttachmentSpecification Attachments;
	};

	class Framebuffer
	{
	public:
		virtual ~Framebuffer() = default;

		virtual void Bind() = 0;
		virtual void Unbind() = 0;

		virtual void Resize(size_t newWidth, size_t newHeight) = 0;

		virtual FramebufferAttachment* QueryAttachment(size_t index = 0) = 0;
		virtual Texture2D* QueryTextureAttachment(size_t index = 0) = 0;

		virtual const FramebufferSpecification& GetSpecification() const = 0;
		virtual FramebufferSpecification& GetSpecification() = 0;

		inline const size_t GetWidth() const { return GetSpecification().Width; }
		inline const size_t GetHeight() const { return GetSpecification().Height; }

		static sRef<Framebuffer> Create(const FramebufferSpecification& spec);
	};
}