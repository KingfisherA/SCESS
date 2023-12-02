#include<pch.h>
#include"3D/Buffers.h"
#include"Platform/OpenGL/OpenGLBuffers.h"
#include"Platform/API Core.h"
#include"Utilities/Memory.h"
//
//

namespace SCESS
{
	sRef<VertexBuffer> VertexBuffer::Create(BufferUsage usage)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			return make <OpenGLVertexBuffer>(usage);
		}
	}

	sRef<VertexBuffer> VertexBuffer::Create(const float* data, size_t floatCount, BufferUsage usage)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			return make <OpenGLVertexBuffer>(usage, data, floatCount);
		}
	}
	//
	//
	sRef<IndexBuffer> IndexBuffer::Create(BufferUsage usage)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			return make <OpenGLIndexBuffer>(usage);
		}
	}

	sRef<IndexBuffer> IndexBuffer::Create(const unsigned int* data, size_t intCount, BufferUsage usage)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			return make <OpenGLIndexBuffer>(usage, data, intCount);
		}
	}
}