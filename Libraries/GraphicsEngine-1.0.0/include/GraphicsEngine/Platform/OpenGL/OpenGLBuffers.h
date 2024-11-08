#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/3D/Buffers.h"
#include"GraphicsEngine/Platform/OpenGL/OpenGLContext.h"

namespace GEng
{
	GENG_API class OpenGLVertexBuffer : public VertexBuffer
	{
	public:
		OpenGLVertexBuffer(BufferUsage usageIn);
		OpenGLVertexBuffer(BufferUsage usageIn, const float *data, size_t floatCount);
		virtual ~OpenGLVertexBuffer();
		//
		virtual void StoreData(const float *data, size_t floatCount) override;
		virtual void UpdateData(const float *data, const size_t offset, const size_t floatCount) override;
		virtual void Bind() override;
		virtual void Unbind() override;
		virtual size_t GetBufferSize() const override;
		virtual size_t GetBufferByteSize() const override;
		//
		const inline GLid GetGLId() const { return m_VBOID; }
		//
	private:
		GLid m_VBOID = 0;
		size_t m_BufferFloatSize = 0;
		BufferUsage m_Usage;
	};

	GENG_API class OpenGLIndexBuffer : public IndexBuffer
	{
	public:
		OpenGLIndexBuffer(BufferUsage usageIn);
		OpenGLIndexBuffer(BufferUsage usageIn, const unsigned int *data, size_t intCount);
		virtual ~OpenGLIndexBuffer();
		//
		virtual void StoreData(const unsigned int *data, size_t intCount) override;
		virtual void Bind() override;
		virtual void Unbind() override;
		virtual size_t GetBufferSize() const override;
		virtual size_t GetBufferByteSize() const override;
		//
		const inline GLid GetGLId() const { return m_VBOID; }
		//
	private:
		GLid m_VBOID = 0;
		size_t m_BufferIntSize = 0;
		BufferUsage m_Usage;
	};
}
