#include<pch.h>
#include "OpenGLVertexArray.h"
#include"Utilities/Assert.h"
//
//

namespace SCESS
{
	OpenGLVertexArray::OpenGLVertexArray() : m_Buffers(), m_IndexBuffer(nullptr)
	{
		glGenVertexArrays(1, &m_VAOId);
	}

	OpenGLVertexArray::~OpenGLVertexArray()
	{
		glDeleteVertexArrays(1, &m_VAOId);
		m_Buffers.clear();
	}

	void OpenGLVertexArray::Bind()
	{
		glBindVertexArray(m_VAOId);
		if (m_IndexBuffer != nullptr) m_IndexBuffer->Bind();
	}

	void OpenGLVertexArray::Unbind()
	{
		if (m_IndexBuffer != nullptr) m_IndexBuffer->Unbind();
		glBindVertexArray(0);
	}

	void OpenGLVertexArray::SubmitArrays(const GLenum drawMode)
	{
		if (IsIndexed())
			glDrawElements(drawMode, GetArrayElementLength(), GL_UNSIGNED_INT, nullptr);
		else
			glDrawArrays(drawMode, 0, GetArrayElementLength());
	}

	Index OpenGLVertexArray::AddVertexBuffer(const sRef<VertexBuffer>& buffer, const VertexSpecification spec)
	{
#ifdef SM_DEBUG
		if (std::dynamic_pointer_cast<OpenGLVertexBuffer>(buffer) == nullptr)
		{
			UNI_FATAL(OpenGLVertexArray::AddVertexBuffer, "Vertex buffer must be created within the OpenGL API in order to be bound to this vertex array.");
			return IndexUndefined;
		}
#endif
		sRef<OpenGLVertexBuffer> glBuffer = std::static_pointer_cast<OpenGLVertexBuffer>(buffer);

		unsigned int attribArray = m_Buffers.size();

		m_Buffers.insert({ attribArray, glBuffer });

		if (m_elementLength == 0) m_elementLength = glBuffer->GetBufferSize();

		glBindVertexArray(m_VAOId);
		glBindBuffer(GL_ARRAY_BUFFER, glBuffer->GetGLId());
		
		glVertexAttribPointer(attribArray, (int)spec.vertexSize, GL_FLOAT, GL_FALSE, spec.vertexStride * sizeof(float), (void*)(spec.vertexOffset * sizeof(float)));
		glVertexAttribDivisor(attribArray, spec.divisor);
		glEnableVertexAttribArray(attribArray);

		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindVertexArray(0);

		return attribArray;
	}

	void OpenGLVertexArray::AddIndexBuffer(const sRef<IndexBuffer>& buffer)
	{
#ifdef SM_DEBUG
		if (buffer == nullptr)
		{
			UNI_FATAL(OpenGLVertexArray::AddIndexBuffer, "Buffer cannot be nullptr.");
			return;
		}
		else if (std::dynamic_pointer_cast<OpenGLIndexBuffer>(buffer) == nullptr)
		{
			UNI_FATAL(OpenGLVertexArray::AddIndexBuffer, "Index buffer must be created within the OpenGL API in order to be bound to this vertex array.");
			return;
		}
#endif

		this->m_IndexBuffer = std::static_pointer_cast<OpenGLIndexBuffer>(buffer);
		m_elementLength = buffer->GetBufferSize();
	}

	size_t OpenGLVertexArray::QueryBufferIndex(const sRef<VertexBuffer>& buffer)
	{
		for (auto& kvPair : m_Buffers)
		{
			if (kvPair.second == buffer) return (size_t)kvPair.first;
		}
		return IndexUndefined;
	}

	sRef<VertexBuffer> OpenGLVertexArray::GetVBuffer(Index index)
	{
#ifdef SM_DEBUG
		//if (index >= m_Buffers.size())
		if(m_Buffers.find(index) == m_Buffers.end())
		{
			UNI_FATAL(OpenGLVertexArray::GetVBuffer, "Index of buffer is not valid");
			return nullptr;
		}
#endif

		return m_Buffers.at(index);
	}

	sRef<IndexBuffer> OpenGLVertexArray::GetIBuffer()
	{
#ifdef SM_DEBUG
		if (m_IndexBuffer == nullptr)
		{
			UNI_FATAL(OpenGLVertexArray::GetIBuffer, "There is no index buffer bound to this vertex array.");
			return nullptr;
		}
#endif
		return m_IndexBuffer;
	}

	bool OpenGLVertexArray::IsIndexed() const
	{
		return m_IndexBuffer != nullptr;
	}
}