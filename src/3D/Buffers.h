#pragma once
#include"Utilities/Memory.h"

namespace SCESS
{
	struct VertexLayout
	{
		unsigned char eleSize = 0;
		size_t offset = 0;
		size_t stride = 0;
	};

	enum class BufferUsage : unsigned char
	{
		StaticUpdate, DynamicUpdate
	};

	class VertexBuffer
	{
	public:
		virtual ~VertexBuffer() {}
		//
		virtual void StoreData(const float* data, const size_t floatCount) = 0;
		virtual void UpdateData(const float* data, const size_t offset, const size_t floatCount) = 0;
		virtual void Bind() = 0;
		virtual void Unbind() = 0;
		virtual size_t GetBufferSize() const = 0;
		virtual size_t GetBufferByteSize() const = 0;
		//
		static sRef<VertexBuffer> Create(BufferUsage usage = BufferUsage::StaticUpdate);
		static sRef<VertexBuffer> Create(const float* data, size_t floatCount, BufferUsage usage = BufferUsage::StaticUpdate);
	};

	class IndexBuffer
	{
	public:
		virtual ~IndexBuffer() {}
		//
		virtual void StoreData(const unsigned int* data, size_t intCount) = 0;
		virtual void Bind() = 0;
		virtual void Unbind() = 0;
		virtual size_t GetBufferSize() const = 0;
		virtual size_t GetBufferByteSize() const = 0;
		//
		static sRef<IndexBuffer> Create(BufferUsage usage = BufferUsage::StaticUpdate);
		static sRef<IndexBuffer> Create(const unsigned int* data, size_t intCount, BufferUsage usage = BufferUsage::StaticUpdate);
	};
}