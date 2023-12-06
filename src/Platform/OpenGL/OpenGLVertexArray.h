#pragma once
#include<pch.h>

#include<unordered_map>

#include"3D/VertexArray.h"
#include"OpenGLBuffers.h"

namespace SCESS
{
	class OpenGLVertexArray : public VertexArray
	{
	public:
		OpenGLVertexArray();
		virtual ~OpenGLVertexArray();

		virtual void Bind() override;
		virtual void Unbind() override;

		virtual void SubmitArrays(const GLenum drawMode = GL_TRIANGLES) override;

		virtual Index AddVertexBuffer(const sRef<VertexBuffer>& buffer, const VertexSpecification spec = { 3,0,0,0 }) override;
		virtual void AddIndexBuffer(const sRef<IndexBuffer>& buffer) override;

		virtual Index QueryBufferIndex(const sRef<VertexBuffer>& buffer) override;

		virtual sRef<VertexBuffer> GetVBuffer(Index index) override;
		virtual sRef<IndexBuffer> GetIBuffer() override;
		virtual bool IsIndexed() const override;

		virtual size_t GetNoOfVBuffers() const override
		{
			return m_Buffers.size();
		}

		const inline GLid GetGLID() const { return m_VAOId; }

	private:
		GLid m_VAOId = 0;
		std::unordered_map<Index, sRef<OpenGLVertexBuffer>> m_Buffers;
		sRef<OpenGLIndexBuffer> m_IndexBuffer;

		
	};
}