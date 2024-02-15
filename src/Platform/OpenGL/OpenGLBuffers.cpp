#include<pch.h>
#include "OpenGLBuffers.h"
//
//

namespace SCESS
{

	static GLenum UsageToGLenum(BufferUsage usage)
	{
		switch (usage)
		{
		case BufferUsage::StaticUpdate:
			return GL_STATIC_DRAW;
		case BufferUsage::DynamicUpdate:
			return GL_DYNAMIC_DRAW;
		default:
			return GL_STATIC_DRAW;
		}
	}

#pragma region OpenGLVertexBuffer Class

	OpenGLVertexBuffer::OpenGLVertexBuffer(BufferUsage usageIn) : m_Usage(usageIn)
	{
		glGenBuffers(1, &m_VBOID);
	}

	OpenGLVertexBuffer::OpenGLVertexBuffer(BufferUsage usageIn, const float* data, size_t floatCount) : m_Usage(usageIn)
	{
		glGenBuffers(1, &m_VBOID);
		this->StoreData(data, floatCount);
	}

	OpenGLVertexBuffer::~OpenGLVertexBuffer()
	{
		glDeleteBuffers(1, &m_VBOID);
		m_VBOID = 0;
	}

	void OpenGLVertexBuffer::StoreData(const float* data, size_t floatCount)
	{
		m_BufferFloatSize = floatCount;
		glBindBuffer(GL_ARRAY_BUFFER, m_VBOID);
		glBufferData(GL_ARRAY_BUFFER, floatCount * sizeof(float), (const void*)data, UsageToGLenum(m_Usage));
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	void OpenGLVertexBuffer::UpdateData(const float* data, const size_t offset, const size_t floatCount)
	{
		glBindBuffer(GL_ARRAY_BUFFER, m_VBOID);
		glBufferSubData(GL_ARRAY_BUFFER, offset, floatCount * sizeof(float), data);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	void OpenGLVertexBuffer::Bind()
	{
		glBindBuffer(GL_ARRAY_BUFFER, m_VBOID);
	}

	void OpenGLVertexBuffer::Unbind()
	{
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	size_t  OpenGLVertexBuffer::GetBufferSize() const
	{
		return m_BufferFloatSize;
	}

	size_t  OpenGLVertexBuffer::GetBufferByteSize() const
	{
		return m_BufferFloatSize * sizeof(float);
	}

#pragma endregion

#pragma region OpenGLIndexBuffer Class

	OpenGLIndexBuffer::OpenGLIndexBuffer(BufferUsage usageIn) : m_Usage(usageIn)
	{
		glGenBuffers(1, &m_VBOID);
	}

	OpenGLIndexBuffer::OpenGLIndexBuffer(BufferUsage usageIn, const unsigned int* data,
		size_t intCount) : m_Usage(usageIn)
	{
		glGenBuffers(1, &m_VBOID);
		StoreData(data, intCount);
	}

	OpenGLIndexBuffer::~OpenGLIndexBuffer()
	{
		glDeleteBuffers(1, &m_VBOID);
		m_VBOID = 0;
	}

	void OpenGLIndexBuffer::StoreData(const unsigned int* data, size_t intCount)
	{
		m_BufferIntSize = intCount;

		if (!m_DataStored)
		{
			m_DataStored = true;

			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_VBOID);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, intCount * sizeof(unsigned int), (const void*)data, UsageToGLenum(m_Usage));
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
		else
		{
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_VBOID);
			glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, intCount * sizeof(unsigned int), data);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		}
	}

	void OpenGLIndexBuffer::Bind()
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_VBOID);
	}

	void OpenGLIndexBuffer::Unbind()
	{
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}

	size_t OpenGLIndexBuffer::GetBufferSize() const
	{
		return m_BufferIntSize;
	}

	size_t OpenGLIndexBuffer::GetBufferByteSize() const
	{
		return m_BufferIntSize * sizeof(unsigned int);
	}


#pragma endregion

}