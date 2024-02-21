#pragma once
#include<array>

#include"Utilities/Memory.h"

namespace SCESS
{

#pragma region Circle Constants
	//
	constexpr float DefaultCircleRadius = 1.0f;
	constexpr unsigned int DefaultCircleSegments = 32;
	//
#pragma endregion

	//
	class Mesh
	{
	public:
		using byte = unsigned char;

		Mesh() : Mesh(0) {}
		Mesh(const size_t bufferByteSize);
		~Mesh();
		//
		float* GetVertexArray(const size_t floatCount = 0);
		float* GetNormalsArray(const size_t floatCount = 0);
		float* GetTexCoordsArray(const size_t floatCount = 0);
		unsigned int* GetIndiciesArray(const size_t intCount = 0);
		//
		inline size_t GetVertexFloatCount() const { return vArrayCount; }
		inline size_t GetNormalsFloatCount() const { return nArrayCount; }
		inline size_t GetTexCoordFloatCount() const { return texArrayCount; }
		inline size_t GetIndexIntCount() const { return iArrayCount; }
		//
		inline bool IsIndexed() const { return iArray != nullptr; }
		//
		unsigned char vertexElementSize = 3;
		//
	private:
		float* vArray = nullptr;
		size_t vArrayCount = 0;
		float* nArray = nullptr;
		size_t nArrayCount = 0;
		float* texArray = nullptr;
		size_t texArrayCount = 0;
		unsigned int* iArray = nullptr;
		size_t iArrayCount = 0;

		byte* byteBuffer = nullptr;
		size_t byteBufferSize = 0;
		size_t currentOffset = 0;
	};

	constexpr double ICO_X = 0.525731112119133606;
	constexpr double ICO_Z = 0.850650808352039932;
	constexpr double ICO_N = 0.0;

	constexpr std::array<double, 36> ICOSAHEDRON_VERTS = {
				ICO_X, -ICO_N, -ICO_Z,
				-ICO_X, -ICO_N, -ICO_Z,
				ICO_X, -ICO_N, ICO_Z,
				-ICO_X, -ICO_N, ICO_Z,
				-ICO_N, -ICO_Z, -ICO_X,
				-ICO_N, -ICO_Z, ICO_X,
				-ICO_N, ICO_Z, -ICO_X,
				-ICO_N, ICO_Z, ICO_X,
				-ICO_Z, -ICO_X, -ICO_N,
				ICO_Z, -ICO_X, -ICO_N,
				-ICO_Z, ICO_X, -ICO_N,
				ICO_Z, ICO_X, -ICO_N
	};
	constexpr std::array<unsigned int, 60> ICOSAHEDRON_TRIS = {
0,4,1 , 0,9,4 , 9,5,4 , 4,5,8 , 4,8,1 ,
 8,10,1 , 8,3,10 , 5,3,8 , 5,2,3 , 2,7,3 ,
 7,10,3 , 7,6,10 , 7,11,6 , 11,0,6 , 0,1,6 ,
 6,1,10 , 9,0,11 , 9,11,2 , 9,2,5 , 7,2,11
	};


	//
	class Geometry final abstract
	{
	public:
		using TriangleI = glm::vec<3, unsigned int>;
		using Vertex = glm::vec3;
		using TexCoord = glm::vec2;
		//
		static glm::vec3 CalculateNormal(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2);
		//
		static sRef<Mesh> CreateUVSphere(float radius, unsigned int lattitudes, unsigned int longitudes);
		static sRef<Mesh> CreateQuadsphere(float radius, unsigned int subdivisions);
		//
		static sRef<Mesh> CreateIcosahedron(float size = 1.0f);
		static sRef<Mesh> CreateIcosphere(float radius, unsigned int subdivisions);
		//
		static sRef<Mesh> CreateCylinder(float radius, float height, size_t segments);
	private:
		//
		static std::vector<unsigned int> DetectSphereWrappedUVCoords(const std::vector<TriangleI>& triangles, const std::vector<glm::vec2>& texCoords, unsigned int hintCount = 0);
		static void FixSphereWrappedUV(const std::vector<unsigned int>& wrapped, std::vector<glm::vec3>& verticies,
			std::vector<TriangleI>& triangles, std::vector<glm::vec2>& texCoords, unsigned int hintCount = 0);
		static void FixPoleVerticies(std::vector<glm::vec3>& verticies, std::vector<glm::vec2>& texCoords, std::vector<TriangleI>& triangles);
		//
	public:
		//
		static  sRef<Mesh> CreateQuad();
		//
		static sRef<Mesh> CreateGrid(float scale, unsigned int slices, unsigned int stacks);
		//
		static sRef<Mesh> CreateCube(float size = 1.0f);
		//
		static sRef<Mesh> CreateCircle(float radius, unsigned int segments);
		static sRef<Mesh> CreateRings(float innerRadius, float outerRadius, unsigned int segments);
		static sRef<Mesh> CreateHalfRings(float innerRadius, float outerRadius, unsigned int segments);
		//
		static sRef<Mesh> CreateTorus(float innerRadius, float outerRadius, float depth, unsigned int sides, unsigned int rings, unsigned short twistyness = 1);
		//
	};
}