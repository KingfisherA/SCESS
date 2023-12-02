#include<pch.h>
#include"VertexArray.h"
#include"Platform/OpenGL/OpenGLVertexArray.h"
#include"Platform/API Core.h"
//
//

namespace SCESS
{
	sRef<VertexArray> VertexArray::Create()
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
			return make <OpenGLVertexArray>();
		}
	}

	sRef<VertexArray> VertexArray::Create(const sRef<Mesh>& meshData)
	{
		switch (PROGRAM_API)
		{
		case API::OpenGL:
		{
			sRef<VertexArray> vArray = make <OpenGLVertexArray>();
			vArray->StoreMesh(meshData);
			return vArray;
		}
		}
	}

	void VertexArray::StoreMesh(const sRef<Mesh>& meshData)
	{
		if (meshData->GetVertexArray() != nullptr) this->AddVertexBuffer(VertexBuffer::Create(meshData->GetVertexArray(), meshData->GetVertexFloatCount()),
			{ meshData->vertexElementSize, 0, 0 });
		if (meshData->GetNormalsArray() != nullptr) this->AddVertexBuffer(VertexBuffer::Create(meshData->GetNormalsArray(), meshData->GetNormalsFloatCount()));
		if (meshData->GetTexCoordsArray() != nullptr) this->AddVertexBuffer(VertexBuffer::Create(meshData->GetTexCoordsArray(), meshData->GetTexCoordFloatCount()),
			{ 2,0,0 });

		if (meshData->IsIndexed()) this->AddIndexBuffer(IndexBuffer::Create(meshData->GetIndiciesArray(), meshData->GetIndexIntCount()));
	}
}
