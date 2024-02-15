#pragma once
#include"Buffers.h"
#include"Utilities/Memory.h"
#include"Maths/Geometry.hpp"

namespace SCESS
{
	struct VertexSpecification
	{
		unsigned char vertexSize = 3;
		//In floats
		unsigned int vertexStride = 0;
		//In floats
		unsigned int vertexOffset = 0;
		//In instances. Specify the number of instances 
		// that will pass between updates of the generic attribute at slot index.
		unsigned int divisor = 0;
	};

	class VertexArray abstract
	{
	public:
		virtual ~VertexArray() {}

		static sRef<VertexArray> Create();
		static sRef<VertexArray> Create(const sRef<Mesh>& meshData);

		virtual void Bind() = 0;
		virtual void Unbind() = 0;

		virtual void SubmitArrays(const GLenum drawMode = GL_TRIANGLES, size_t instances = 1) = 0;

		virtual Index AddVertexBuffer(const sRef<VertexBuffer>& buffer, const VertexSpecification spec = { 3,0,0 }) = 0;
		virtual void AddIndexBuffer(const sRef<IndexBuffer>& buffer) = 0;

		virtual void StoreMesh(const sRef<Mesh>& meshData);

		virtual Index QueryBufferIndex(const sRef<VertexBuffer>& buffer) = 0;

		virtual sRef<VertexBuffer> GetVBuffer(Index index) = 0;
		virtual sRef<IndexBuffer> GetIBuffer() = 0;
		virtual bool IsIndexed() const = 0;

		virtual size_t GetNoOfVBuffers() const = 0;

		void UpdateElementLength(const size_t newEleLength)
		{
			m_elementLength = newEleLength;
		}
		size_t GetArrayElementLength() const { return m_elementLength; }

	protected:
		size_t m_elementLength = 0;
	};
}
