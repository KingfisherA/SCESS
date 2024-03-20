#include<pch.h>

#include<vector>
#include<map>
#include<set>

#include"Geometry.hpp"
#include"Utilities/Assert.h"
#include"Maths/Maths.hpp"
#include<numbers>

namespace SCESS
{
	//
	struct hash_pair
	{
		template <class T1, class T2>
		size_t operator()(const std::pair<T1, T2>& p) const
		{
			auto hash1 = std::hash<T1>{}(p.first);
			auto hash2 = std::hash<T2>{}(p.second);
			return hash1 ^ hash2;
		}
	};
	//
#pragma region Mesh
	//
	Mesh::Mesh(const size_t bufferByteSize) : byteBufferSize(bufferByteSize)
	{
		byteBuffer = (byte*)malloc(bufferByteSize);
		UNI_ASSERT(byteBuffer == nullptr, SCESS::Mesh::Mesh, "Could not dynamically allocate memory!");
	}
	Mesh::~Mesh()
	{
		free(byteBuffer);
		byteBuffer = nullptr;
	}
	//
	float* Mesh::GetVertexArray(const size_t floatCount)
	{
		if (floatCount == 0 || vArray != nullptr)
			return vArray;

		vArray = (float*)(byteBuffer + currentOffset);
		vArrayCount = floatCount;
		currentOffset += floatCount * sizeof(float);
		return vArray;
	}
	float* Mesh::GetNormalsArray(const size_t floatCount)
	{
		if (floatCount == 0 || nArray != nullptr)
			return nArray;

		nArray = (float*)(byteBuffer + currentOffset);
		nArrayCount = floatCount;
		currentOffset += floatCount * sizeof(float);
		return nArray;
	}
	float* Mesh::GetTexCoordsArray(const size_t floatCount)
	{
		if (floatCount == 0 || texArray != nullptr)
			return texArray;

		texArray = (float*)(byteBuffer + currentOffset);
		texArrayCount = floatCount;
		currentOffset += floatCount * sizeof(float);
		return texArray;
	}
	unsigned int* Mesh::GetIndiciesArray(const size_t intCount)
	{
		if (intCount == 0 || iArray != nullptr)
			return iArray;

		iArray = (unsigned int*)(byteBuffer + currentOffset);
		iArrayCount = intCount;
		currentOffset += intCount * sizeof(unsigned int);
		return iArray;
	}
	//
#pragma endregion
// 
//
#pragma region Geometry
	//
	glm::vec3 Geometry::CalculateNormal(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
	{
		//	glm::vec3 A = v1 - v0;
		//	glm::vec3 B = v2 - v0;
		//	glm::vec3 N = glm::cross(A, B);
		//	float sinAlpha = glm::length(N) / (glm::length(A) * glm::length(B));
		//	return glm::normalize(N) * sinAlpha;
		return glm::normalize(glm::cross(v1 - v0, v2 - v0));
	}
	//
	sRef<Mesh> Geometry::CreateUVSphere(float radius, unsigned int lattitudes, unsigned int longitudes)
	{
		UNI_ASSERT(lattitudes < 3, Geometry::CreateSphere, "The minimum number of slices is 3");
		UNI_ASSERT(longitudes < 3, Geometry::CreateSphere, "The minimum number of stacks is 3");

		unsigned int numberOfVerts = ((size_t)lattitudes) * (longitudes - 1) + 2;
		std::vector<Vertex> verticies;
		verticies.reserve(numberOfVerts);
		std::vector<glm::vec2> texCoord;
		texCoord.reserve(numberOfVerts);
		std::vector<TriangleI> triangles;
		triangles.reserve(((size_t)longitudes - 2) * lattitudes * 2 + 2 * lattitudes);

		unsigned int i = 0;

		{
			unsigned int j = 0;
			glm::vec3 currentVertex;

			//Write data for top most vertex
			verticies.emplace_back(0.0f, 1.0f, 0.0f);
			texCoord.emplace_back(0.5f, 0.0f);

			for (i = 1; i < longitudes; i++)
			{
				//Compute the V texture coordinate
				float V = float(i) / float(longitudes);
				//Compute the stack's angle
				float Phi = V * SCESS::Maths::PI;

				for (j = 0; j < lattitudes; j++)
				{
					//Compute the U texture coordinate
					float U = float(j) / float(lattitudes);
					//Compute the slice's angle
					float Theta = U * SCESS::Maths::Tau;

					currentVertex.x = cos(Theta) * sin(Phi);
					currentVertex.y = cos(Phi);
					currentVertex.z = sin(Theta) * sin(Phi);

					verticies.push_back(currentVertex);

					texCoord.emplace_back(
						((atan2(currentVertex.x, currentVertex.z) + Maths::PI) / Maths::PI) * 0.5f,
						((acos(currentVertex.y) + Maths::PI) / Maths::PI) - 1.0f);
				}
			}

			//Write data for bottom most vertex
			verticies.emplace_back(0.0f, -1.0f, 0.0f);
			texCoord.emplace_back(0.5f, 1.0f);
		}

		{
			const unsigned int topTriangleCount = lattitudes;
			const unsigned int mainFaceCount = (longitudes - 2) * lattitudes;

			for (i = 0; i < topTriangleCount; i++)
			{
				triangles.emplace_back(
					i == topTriangleCount - 1 ? 1 : i + 2,
					i + 1,
					0);
			}

			//Face count instead of triangle count as two triangles are written per iteration
			for (i = topTriangleCount; i < mainFaceCount + topTriangleCount; i++)
			{
				if ((i + 1) % (lattitudes) == 0)
				{
					//First triangle's indices
					triangles.emplace_back(
						i - lattitudes + 2,
						i + 1,
						i - lattitudes + 1);
					//Second triangle's indices
					triangles.emplace_back(
						i - lattitudes + 2,
						i - lattitudes + 1,
						i - lattitudes - (lattitudes - 1) + 1);

					continue;
				}

				//First triangle's indices
				triangles.emplace_back(
					i + 2,
					i + 1,
					i - lattitudes + 1);
				//Second triangle's indices
				triangles.emplace_back(
					i - lattitudes + 1,
					i - lattitudes + 2,
					i + 2);
			}

			for (i = 0; i < topTriangleCount; i++)
			{
				unsigned int triIndex = i + verticies.size() - 1 - lattitudes;
				triangles.emplace_back(
					i == topTriangleCount - 1 ? triIndex - (lattitudes - 1) : triIndex + 1,
					verticies.size() - 1,
					triIndex);
			}

			std::vector<unsigned int> wrappedCoords = move(DetectSphereWrappedUVCoords(triangles, texCoord, longitudes));
			FixSphereWrappedUV(wrappedCoords, verticies, triangles, texCoord);
		}

		sRef<Mesh>sphereMesh = make <Mesh>(
			verticies.size() * sizeof(glm::vec3) * 2 + //*2 for normals
			texCoord.size() * sizeof(glm::vec2) +
			triangles.size() * sizeof(TriangleI)
		);

		float* verticiesSpan = sphereMesh->GetVertexArray(verticies.size() * 3);
		float* normalsSpan = sphereMesh->GetNormalsArray(verticies.size() * 3);
		float* texCoordsSpan = sphereMesh->GetTexCoordsArray(texCoord.size() * 2);
		unsigned int* indiciesSpan = sphereMesh->GetIndiciesArray(triangles.size() * 3);

		memcpy(normalsSpan, &verticies[0], verticies.size() * sizeof(glm::vec3));

		for (glm::vec3& vertex : verticies)
		{
			vertex *= radius;
		}

		memcpy(verticiesSpan, &verticies[0], verticies.size() * sizeof(glm::vec3));
		memcpy(texCoordsSpan, &texCoord[0], texCoord.size() * sizeof(glm::vec2));
		memcpy(indiciesSpan, &triangles[0], triangles.size() * sizeof(TriangleI));

		return sphereMesh;
	}
	//
	sRef<Mesh>Geometry::CreateQuadsphere(float radius, unsigned int subdivisions)
	{
		unsigned int noOfVerts = subdivisions * subdivisions * 6;
		std::vector<glm::vec3> verticiesArray(noOfVerts);
		std::vector<TriangleI> trianglesArray((subdivisions - 1) * (subdivisions - 1) * 2 * 6);
		std::vector<glm::vec2> texCoordsArray(noOfVerts);

		constexpr glm::vec3 localUps[] = {
			glm::vec3{0.0f, 1.0f, 0.0f},
			glm::vec3{0.0f, -1.0f, 0.0f},

			glm::vec3{1.0f, 0.0f, 0.0f},
			glm::vec3{-1.0f, 0.0f, 0.0f},

			glm::vec3{0.0f, 0.0f, 1.0f},
			glm::vec3{0.0f, 0.0f, -1.0f}
		};

		glm::vec3 axisA;
		glm::vec3 axisB;

		unsigned int triIndex = 0;

		for (unsigned int faces = 0; faces < 6; faces++)
		{
			axisA.x = localUps[faces].y;
			axisA.y = localUps[faces].z;
			axisA.z = localUps[faces].x;

			axisB = glm::cross(localUps[faces], axisA);

			for (unsigned int y = 0; y < subdivisions; y++)
			{
				for (unsigned int x = 0; x < subdivisions; x++)
				{
					glm::vec2 percent(x, y);
					percent /= (float)(subdivisions - 1);

					unsigned int index = x + y * subdivisions + faces * subdivisions * subdivisions;
					verticiesArray[index] = glm::normalize(
						localUps[faces] +
						(percent.x - 0.5f) * 2.0f * axisA +
						(percent.y - 0.5f) * 2.0f * axisB);

					if (x != subdivisions - 1 && y != subdivisions - 1)
					{
						//trianglesArray[triIndex] = index;
						//trianglesArray[triIndex + 1] = index + subdivisions + 1;
						//trianglesArray[triIndex + 2] = index + subdivisions;

						//trianglesArray[triIndex + 3] = index;
						//trianglesArray[triIndex + 4] = index + 1;
						//trianglesArray[triIndex + 5] = index + subdivisions + 1;

						trianglesArray[triIndex].x = index;
						trianglesArray[triIndex].y = index + subdivisions + 1;
						trianglesArray[triIndex].z = index + subdivisions;

						trianglesArray[triIndex + 1].x = index;
						trianglesArray[triIndex + 1].y = index + 1;
						trianglesArray[triIndex + 1].z = index + subdivisions + 1;

						triIndex += 2;
					}
				}
			}
		}

		for (unsigned int i = 0; i < verticiesArray.size(); i++)
		{
			const glm::vec3& normal = verticiesArray[i];

			texCoordsArray[i].x = ((atan2(normal.x, normal.z) + Maths::PI) / Maths::PI) * 0.5f;
			texCoordsArray[i].y = ((acos(normal.y) + Maths::PI) / Maths::PI) - 1.0;
		}

		//std::vector<unsigned int> wrappedCoords = move(DetectSphereWrappedUVCoords(trianglesArray, texCoordsArray, subdivisions * 4));
		//FixSphereWrappedUV(wrappedCoords, verticiesArray, trianglesArray, texCoordsArray, unsigned int(subdivisions * 4.5f));

		sRef<Mesh> sphereMesh = make <Mesh>(
			verticiesArray.size() * sizeof(glm::vec3) * 2 /*For normals array */ +
			trianglesArray.size() * sizeof(TriangleI) +
			texCoordsArray.size() * sizeof(glm::vec2)
		);

		float* vertexArrayS = sphereMesh->GetVertexArray(verticiesArray.size() * 3);
		unsigned int* indiciesArrayS = sphereMesh->GetIndiciesArray(trianglesArray.size() * 3);
		float* normalsArrayS = sphereMesh->GetNormalsArray(verticiesArray.size() * 3);
		float* texCoordArrayS = sphereMesh->GetTexCoordsArray(texCoordsArray.size() * 2);

		memcpy(normalsArrayS, &verticiesArray[0], verticiesArray.size() * sizeof(glm::vec3));

		for (glm::vec3& vertex : verticiesArray)
		{
			vertex *= radius;
		}

		memcpy(vertexArrayS, &verticiesArray[0], verticiesArray.size() * sizeof(glm::vec3));
		memcpy(indiciesArrayS, &trianglesArray[0].x, trianglesArray.size() * sizeof(TriangleI));
		memcpy(texCoordArrayS, &texCoordsArray[0], texCoordsArray.size() * sizeof(glm::vec2));

		return sphereMesh;
	}
	//
	sRef<Mesh>Geometry::CreateIcosahedron(float size)
	{
		sRef<Mesh>icosahedronMesh = make <Mesh>(
			ICOSAHEDRON_VERTS.size() * sizeof(float) * 2 +
			ICOSAHEDRON_TRIS.size() * sizeof(unsigned int)
		);

		float* Verticies = icosahedronMesh->GetVertexArray(ICOSAHEDRON_VERTS.size());
		float* Normals = icosahedronMesh->GetNormalsArray(ICOSAHEDRON_VERTS.size());
		unsigned int* Indicies = icosahedronMesh->GetIndiciesArray(ICOSAHEDRON_TRIS.size());

		unsigned int i = 0;

		for (; i < ICOSAHEDRON_VERTS.size(); i++)
		{
			Normals[i] = (float)ICOSAHEDRON_VERTS[i];
			Verticies[i] = (float)ICOSAHEDRON_VERTS[i] * size;
		}
		for (i = 0; i < ICOSAHEDRON_TRIS.size(); i++)
		{
			Indicies[i] = (unsigned int)ICOSAHEDRON_TRIS[i];
		}

		return icosahedronMesh;
	}
	//
	sRef<Mesh>Geometry::CreateIcosphere(float radius, unsigned int subdivisions)
	{
		auto PrecomputeNoTriangles = [](unsigned int subdivisions) -> unsigned int
			{ return 20 * Maths::IntegerPow<unsigned int>(4, subdivisions); };

		auto PrecomputeVertexAmount = [](unsigned int subdivisions)->unsigned int
			{
				unsigned int result = 12;

				for (unsigned int i = 1; i <= subdivisions; i++)
				{
					result += (30 * Maths::IntegerPow((unsigned int)4, i - 1));
				}
				return result;
			};

		std::vector<Vertex> icoVerticies;
		icoVerticies.reserve(PrecomputeVertexAmount(subdivisions));
		std::vector<TriangleI> icoTriangles;
		icoTriangles.reserve(PrecomputeNoTriangles(subdivisions));
		std::vector<TriangleI> icoTrianglesBuffer;
		icoTrianglesBuffer.reserve(subdivisions > 0 ? PrecomputeNoTriangles(subdivisions - 1) : 0);

		unsigned int i = 0;

		//Write Icosahedron Data
		for (i = 0; i < ICOSAHEDRON_VERTS.size() / 3; i++)
		{
			icoVerticies.emplace_back(
				ICOSAHEDRON_VERTS[i * 3 + 0],
				ICOSAHEDRON_VERTS[i * 3 + 1],
				ICOSAHEDRON_VERTS[i * 3 + 2]);
		}
		for (i = 0; i < ICOSAHEDRON_TRIS.size() / 3; i++)
		{
			icoTriangles.emplace_back(
				ICOSAHEDRON_TRIS[i * 3 + 0],
				ICOSAHEDRON_TRIS[i * 3 + 1],
				ICOSAHEDRON_TRIS[i * 3 + 2]);
		}

		std::unordered_map<std::pair<unsigned int, unsigned int>, unsigned int, hash_pair> lookupTable;

		//Lamda functions for subdividing icosahedron
		auto CalculateVertex = [&lookupTable, &icoVerticies](unsigned int first, unsigned int second) -> unsigned int
			{
				std::pair<unsigned int, unsigned int> key(first, second);
				if (key.first > key.second)
				{
					std::swap(key.first, key.second);
				}

				auto inserted = lookupTable.insert({ key, icoVerticies.size() });
				if (inserted.second)
				{
					const glm::vec3& edge0 = icoVerticies[first];
					const glm::vec3& edge1 = icoVerticies[second];
					icoVerticies.push_back(glm::normalize(edge0 + edge1));
				}

				return inserted.first->second;
			};
		auto Subdivide = [&CalculateVertex, &icoTriangles, &icoTrianglesBuffer, &icoVerticies, &lookupTable]() ->void
			{
				TriangleI midPointTri;
				lookupTable.clear();
				icoTriangles.clear();

				for (TriangleI& currentTri : icoTrianglesBuffer)
				{
					for (int edge = 0; edge < 3; ++edge)
					{
						midPointTri[edge] = CalculateVertex(currentTri[edge], currentTri[(edge + 1) % 3]);
					}

					icoTriangles.emplace_back(currentTri[0], midPointTri[0], midPointTri[2]);
					icoTriangles.emplace_back(currentTri[1], midPointTri[1], midPointTri[0]);
					icoTriangles.emplace_back(currentTri[2], midPointTri[2], midPointTri[1]);
					icoTriangles.emplace_back(midPointTri[0], midPointTri[1], midPointTri[2]);
				}
			};

		if (subdivisions > 0)
		{
			icoTrianglesBuffer = icoTriangles;
		}

		//Perform subdivisions
		for (i = 0; i < subdivisions; ++i)
		{
			Subdivide();

			if (i != subdivisions - 1)
			{
				icoTrianglesBuffer = icoTriangles;
			}
		}

		icoTrianglesBuffer.clear();

		//Create an array to hold texture coordinates
		std::vector<glm::vec2> texCoords;
		texCoords.reserve(icoVerticies.size());
		//Populate texture coordinate array
		for (i = 0; i < icoVerticies.size(); i++)
		{
			const glm::vec3& normal = icoVerticies[i];

			texCoords[i].x = ((atan2(normal.x, normal.z) + Maths::PI) / Maths::PI) * 0.5f;
			texCoords[i].y = ((acos(normal.y) + Maths::PI) / Maths::PI) - 1.0;
		}

		//Fix texture seam issue by duplicating verticies and correcting texture coordinates
		std::vector<unsigned int> wrappedCoords = move(DetectSphereWrappedUVCoords(icoTriangles, texCoords));
		FixSphereWrappedUV(wrappedCoords, icoVerticies, icoTriangles, texCoords);

		//Create a mesh object to store data
		sRef<Mesh> icoSphereMesh = make <Mesh>(
			icoVerticies.size() * sizeof(glm::vec3) * 2 +
			icoTriangles.size() * sizeof(TriangleI) +
			texCoords.size() * 2 * sizeof(float)
		);

		float* icosphereVerts = icoSphereMesh->GetVertexArray(icoVerticies.size() * 3);
		float* icosphereNorms = icoSphereMesh->GetNormalsArray(icoVerticies.size() * 3);
		unsigned int* icosphereIndicies = icoSphereMesh->GetIndiciesArray(icoTriangles.size() * 3);
		float* icosphereUV = icoSphereMesh->GetTexCoordsArray(icoVerticies.size() * 2);

		//Copy normalised verticies first
		memcpy(icosphereNorms, &icoVerticies[0], icoVerticies.size() * sizeof(glm::vec3));

		//Scale all verticies by radius
		for (glm::vec3& vert : icoVerticies)
		{
			vert *= radius;
		}

		//Copy verticies, UV coordinates and triangles (indicies) into data buffer
		memcpy(icosphereVerts, &icoVerticies[0], icoVerticies.size() * sizeof(glm::vec3));
		memcpy(icosphereIndicies, &icoTriangles[0], icoTriangles.size() * sizeof(TriangleI));
		memcpy(icosphereUV, &texCoords[0], texCoords.size() * sizeof(glm::vec2));

		//Return icosphere mesh
		return icoSphereMesh;
	}
	//
	sRef<Mesh> Geometry::CreateCylinder(float radius, float height, size_t segments)
	{
		std::vector<glm::vec3> vertices(2 * segments + 2);
		//std::vector<glm::vec3> normals(2 * segments + 2);
		size_t noTriangles = 4 * segments;
		std::vector<unsigned int> indicies;
		indicies.reserve(noTriangles * 3);

		vertices[0] = { 0.0f, 0.0f, height * 0.5f };
		for (size_t i = 0; i < segments; i++)
		{
			float theta = 2.0f * std::numbers::pi * float(i) / float(segments);
			glm::vec3 v1(radius * cos(theta), radius * sin(theta), height * 0.5f);
			glm::vec3 v2(v1.x, v1.y, -height * 0.5f);
			vertices[i + 1] = v1;
			//normals[i + 1] = glm::vec3(v1.x, v2.y, 0.0f);
			vertices[i + 1 + segments] = v2;
			//normals[i + 1 + segments] = glm::vec3(v1.x, v2.y, 0.0f);
		}
		vertices[2 * segments + 1] = { 0.0f, 0.0f, -height * 0.5f };

		for (size_t i = 0; i < segments; i++)
		{
			indicies.push_back(i + 1);
			indicies.push_back(0);
			indicies.push_back(i == segments - 1 ? 1 : i + 2);
		}

		for (size_t i = 0; i < segments; i++)
		{
			indicies.push_back(i == segments - 1 ? 1 : i + 2);
			indicies.push_back(i + 1);
			indicies.push_back(segments + i + 1);

			indicies.push_back(segments + i + 1);
			indicies.push_back(i == segments - 1 ? segments + 1 : segments + i + 2);
			indicies.push_back(i == segments - 1 ? 1 : i + 2);
		}

		for (size_t i = 0; i < segments; i++)
		{
			indicies.push_back(i + segments + 1);
			indicies.push_back(2 * segments + 1);
			indicies.push_back(i == segments - 1 ? segments + 1 : i + segments + 2);
		}
		
		sRef<Mesh> cylinderMesh = std::make_shared<Mesh>(vertices.size() * sizeof(glm::vec3) +
			indicies.size() * sizeof(unsigned int));

		float* vertexArray = cylinderMesh->GetVertexArray(vertices.size() * 3);
		unsigned int* indicesArray = cylinderMesh->GetIndiciesArray(indicies.size());

		memcpy(vertexArray, vertices.data(), vertices.size() * sizeof(glm::vec3));
		memcpy(indicesArray, indicies.data(), indicies.size() * sizeof(unsigned int));

		return cylinderMesh;
	}
	//
	sRef<Mesh>Geometry::CreateQuad()
	{
		constexpr float quadVerticies[] = { -1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 1.0f, -1.0f, 1.0f };
		constexpr unsigned int quadIndicies[] = { 0, 1, 3, 3, 1, 2 };

		constexpr unsigned int vertexCount = sizeof(quadVerticies) / sizeof(float);
		constexpr unsigned int  indiciesCount = sizeof(quadIndicies) / sizeof(unsigned int);

		sRef<Mesh> quadMesh = make <Mesh>(vertexCount * sizeof(float) + indiciesCount * sizeof(unsigned int));
		memcpy(quadMesh->GetVertexArray(vertexCount), &quadVerticies, vertexCount * sizeof(float));
		memcpy(quadMesh->GetIndiciesArray(indiciesCount), &quadIndicies, indiciesCount * sizeof(unsigned int));
		quadMesh->vertexElementSize = 2;

		return quadMesh;
	}
	//
	sRef<Mesh>Geometry::CreateGrid(float scale, unsigned int slices, unsigned int stacks)
	{
		//unsigned int vertexCount = (slices + 1) * (stacks + 1) * 3;

		std::vector<float> vertices;
		std::vector<unsigned int> indices;

		unsigned int i = 0, j = 0;

		for (j = 0; j <= slices; ++j)
		{
			for (i = 0; i <= slices; ++i)
			{
				float x = ((float)i / (float)slices) * scale;
				float y = 0;
				float z = ((float)j / (float)slices) * scale;
				vertices.push_back(x);
				vertices.push_back(y);
				vertices.push_back(z);
			}
		}

		for (j = 0; j < slices; ++j)
		{
			for (i = 0; i < slices; ++i)
			{
				unsigned int  row1 = j * (slices + 1);
				unsigned int row2 = (j + 1) * (slices + 1);


				indices.push_back(row1 + i);
				indices.push_back(row1 + i + 1);
				indices.push_back(row1 + i + 1);
				indices.push_back(row2 + i + 1);

				indices.push_back(row2 + i + 1);
				indices.push_back(row2 + i);
				indices.push_back(row2 + i);
				indices.push_back(row1 + i);
			}
		}

		sRef<Mesh>gridMesh = make <Mesh>(vertices.size() * sizeof(float) + indices.size() * sizeof(unsigned int));
		float* vArraySpan = gridMesh->GetVertexArray(vertices.size());
		unsigned int* iArraySpan = gridMesh->GetIndiciesArray(indices.size());

		memcpy(vArraySpan, &vertices[0], vertices.size() * sizeof(float));
		memcpy(iArraySpan, &indices[0], indices.size() * sizeof(unsigned int));

		return gridMesh;
	}
	//
	sRef<Mesh>Geometry::CreateCube(float size)
	{
		constexpr glm::vec3 cubeVerticies[] = {
			//Front
{-1,  -1, 1},
{1,  -1, 1},
{1,  1, 1},
{1,  1, 1},
{-1,  1, 1},
{-1,  -1, 1},

//Right
{1,  -1, 1},
{1,  -1, -1},
{1,  1, -1},
{1,  1, -1},
{1,  1, 1},
{1,  -1, 1},

//Back
{-1,  1, -1},
{1,  1, -1},
{1,  -1, -1},
{1,  -1, -1},
{-1,  -1, -1},
{-1,  1, -1},

//Left
{-1,  -1, -1},
{-1,  -1, 1},
{-1,  1, 1},
{-1,  1, 1},
{-1,  1, -1},
{-1,  -1, -1},

//Bottom
{-1,  -1, -1},
{1,  -1, -1},
{1,  -1, 1},
{1,  -1, 1},
{-1,  -1, 1},
{-1,  -1, -1},

//Top
{-1,  1, 1},
{1,  1, 1},
{1,  1, -1},
{1,  1, -1},
{-1,  1, -1},
{-1,  1, 1}
		};
		constexpr glm::vec3 cubeNormals[] = {
	{0, 0, 1},
{0, 0, 1},
{0, 0, 1},
{0, 0, 1},
{0, 0, 1},
{0, 0, 1},

{1, 0, 0},
{1, 0, 0},
{1, 0, 0},
{1, 0, 0},
{1, 0, 0},
{1, 0, 0},

{0, 0, -1},
{0, 0, -1},
{0, 0, -1},
{0, 0, -1},
{0, 0, -1},
{0, 0, -1},

{-1, 0, 0},
{-1, 0, -0},
{-1, 0, 0},
{-1, 0, -0},
{-1, 0, 0},
{-1, 0, -0},

{0, -1, 0},
{-0, -1, 0},
{0, -1, 0},
{-0, -1, 0},
{0, -1, 0},
{-0, -1, 0},

{-0, 1, 0},
{0, 1, 0},
{-0, 1, 0},
{0, 1, 0},
{-0, 1, 0},
{0, 1, 0}
		};
		constexpr glm::vec2 cubeUV[] = {
			//Front
			{0.0f, 0.0f},
			{1.0f, 0.0f},
			{1.0f, -1.0f},
			{1.0f, -1.0f},
			{0.0f, -1.0f},
			{0.0f, 0.0f},

			//Right
			{0.0f, 0.0f},
			{1.0f, 0.0f},
			{1.0f, -1.0f},
			{1.0f, -1.0f},
			{0.0f, -1.0f},
			{0.0f, 0.0f},

			//Back
			{1.0f, -1.0f},
			{0.0f, -1.0f},
			{0.0f, 0.0f},
			{0.0f, 0.0f},
			{1.0f, 0.0f},
			{1.0f, -1.0f},

			//Left
			{0.0f, 0.0f},
			{1.0f, 0.0f},
			{1.0f, -1.0f},
			{1.0f, -1.0f},
			{0.0f, -1.0f},
			{0.0f, 0.0f},

			//Bottom
			{0.0f, 0.0f},
			{-1.0f, 0.0f},
			{-1.0f, 1.0f},
			{-1.0f, 1.0f},
			{0.0f, 1.0f},
			{0.0f, 0.0f},

			//Top
			{0.0f, 0.0f},
			{-1.0f, 0.0f},
			{-1.0f, 1.0f},
			{-1.0f, 1.0f},
			{0.0f, 1.0f},
			{0.0f, 0.0f}
		};

		sRef<Mesh> cubeMesh = make <Mesh>(
			sizeof(cubeVerticies) +
			sizeof(cubeNormals) +
			sizeof(cubeUV)
		);

		float* vertSpan = cubeMesh->GetVertexArray(sizeof(cubeVerticies) / sizeof(float));
		float* normalSpan = cubeMesh->GetNormalsArray(sizeof(cubeNormals) / sizeof(float));
		float* texCoordSpan = cubeMesh->GetTexCoordsArray(sizeof(cubeUV) / sizeof(float));

		memcpy(vertSpan, &cubeVerticies, sizeof(cubeVerticies));
		memcpy(normalSpan, &cubeNormals, sizeof(cubeNormals));
		memcpy(texCoordSpan, &cubeUV, sizeof(cubeUV));

		for (unsigned int i = 0; i < sizeof(cubeVerticies) / sizeof(float); i++)
		{
			vertSpan[i] *= size;
		}

		return cubeMesh;
	}
	//
	sRef<Mesh>Geometry::CreateCircle(float radius, unsigned int segments)
	{
		sRef<Mesh> circleMesh = make<Mesh>((size_t)segments * 2 * sizeof(float));
		float* Verticies = circleMesh->GetVertexArray((size_t)segments * 2);

		float theta = 0.0f;
		for (unsigned int i = 0; i < segments; i++)
		{
			theta = 2.0f * Maths::PI * float(i) / float(segments);

			Verticies[(size_t)i * 2 + 0] = radius * cos(theta);
			Verticies[(size_t)i * 2 + 1] = radius * sin(theta);
		}

		return circleMesh;
	}
	//
	sRef<Mesh>Geometry::CreateRings(float innerRadius, float outerRadius, unsigned int segments)
	{
		unsigned int& faceCount = segments;
		size_t vertexCount = (size_t)faceCount * 4;
		size_t indexCount = (size_t)segments * 6;

		sRef<Mesh> ringMesh = make <Mesh>(vertexCount * sizeof(float) + indexCount * sizeof(unsigned int));

		float* Verticies = ringMesh->GetVertexArray(vertexCount);
		unsigned int* Indicies = ringMesh->GetIndiciesArray(indexCount);

		float theta;
		unsigned int i = 0;
		for (i = 0; i < faceCount; i++)
		{
			theta = 2.0f * Maths::PI * float(i) / float(segments);

			Verticies[(size_t)i * 4 + 0] = innerRadius * cos(theta);
			Verticies[(size_t)i * 4 + 1] = innerRadius * sin(theta);

			Verticies[(size_t)i * 4 + 2] = outerRadius * cos(theta);
			Verticies[(size_t)i * 4 + 3] = outerRadius * sin(theta);
		}

		for (i = 0; i < segments; i++)
		{
			if (i + 1 == segments)
			{
				Indicies[(size_t)i * 6 + 0] = i * 2 + 0;
				Indicies[(size_t)i * 6 + 1] = 1;
				Indicies[(size_t)i * 6 + 2] = i * 2 + 1;

				Indicies[(size_t)i * 6 + 3] = i * 2 + 0;
				Indicies[(size_t)i * 6 + 4] = 0;
				Indicies[(size_t)i * 6 + 5] = 1;
			}
			else
			{
				Indicies[(size_t)i * 6 + 0] = i * 2 + 0;
				Indicies[(size_t)i * 6 + 1] = i * 2 + 3;
				Indicies[(size_t)i * 6 + 2] = i * 2 + 1;

				Indicies[(size_t)i * 6 + 3] = i * 2 + 0;
				Indicies[(size_t)i * 6 + 4] = i * 2 + 2;
				Indicies[(size_t)i * 6 + 5] = i * 2 + 3;
			}
		}

		return ringMesh;
	}
	//
	sRef<Mesh> Geometry::CreateHalfRings(float innerRadius, float outerRadius, unsigned int segments)
	{
		unsigned int faceCount = segments / 2 + 1;
		size_t vertexCount = (size_t)faceCount * 4;
		size_t indexCount = (size_t)(faceCount - 1) * 6;

		sRef<Mesh> ringMesh = make <Mesh>(vertexCount * sizeof(float) + indexCount * sizeof(unsigned int));

		float* Verticies = ringMesh->GetVertexArray(vertexCount);
		unsigned int* Indicies = ringMesh->GetIndiciesArray(indexCount);

		float theta;
		unsigned int i = 0;
		for (i = 0; i < faceCount; i++)
		{
			theta = 2.0f * Maths::PI * float(i) / float(segments);

			Verticies[(size_t)i * 4 + 0] = innerRadius * cos(theta);
			Verticies[(size_t)i * 4 + 1] = innerRadius * sin(theta);

			Verticies[(size_t)i * 4 + 2] = outerRadius * cos(theta);
			Verticies[(size_t)i * 4 + 3] = outerRadius * sin(theta);
		}

		for (i = 0; i < faceCount - 1; i++)
		{
			Indicies[(size_t)i * 6 + 0] = i * 2 + 0;
			Indicies[(size_t)i * 6 + 1] = i * 2 + 3;
			Indicies[(size_t)i * 6 + 2] = i * 2 + 1;

			Indicies[(size_t)i * 6 + 3] = i * 2 + 0;
			Indicies[(size_t)i * 6 + 4] = i * 2 + 2;
			Indicies[(size_t)i * 6 + 5] = i * 2 + 3;
		}

		return ringMesh;
	}
	//
	sRef<Mesh> Geometry::CreateTorus(float innerRadius, float outerRadius, float depth, unsigned int sides, unsigned int rings, unsigned short twistyness)
	{
		//Array for final verticies and normal vectors
		unsigned int elementCount = (sides + 1) * ((rings + 1) * (twistyness + 1));
		size_t vertexCount = ((size_t)elementCount * 3);
		size_t texCoordsCount = ((size_t)elementCount * 2);
		size_t indiciesCount = ((size_t)elementCount * 3) - 6;

		sRef<Mesh>torusMesh = make <Mesh>(
			vertexCount * sizeof(float) * 2 + //For verticies and normals
			texCoordsCount * sizeof(float) +
			indiciesCount * sizeof(unsigned int)
		);

		float* Vertices = torusMesh->GetVertexArray(vertexCount);
		float* Normals = torusMesh->GetNormalsArray(vertexCount);
		float* TextureCoordinates = torusMesh->GetTexCoordsArray(texCoordsCount);
		unsigned int* Indicies = torusMesh->GetIndiciesArray(indiciesCount);

		//Incremental varibles
		unsigned int i = 0, j = 0;
		unsigned int a = 0;

		//Tesselation "twistyness". Varible used for correct quad tesselation
		short k = twistyness;
		//Angles for torus calculation
		float theta = 0, phi = 0;

		for (i = 0; i <= sides; i++)
		{
			for (j = 0; j <= rings; j++)
			{
				for (k = twistyness; k >= 0; k--)
				{
					//	theta = (i + k) % rings+0.5f;
					//phi = j % rings;
					theta = (float)fmod(i + k, float(rings) + 0.5f);
					phi = (float)fmod(j, rings);

					Vertices[(size_t)a * 3 + 0] = (outerRadius + innerRadius * cos(theta * Maths::Tau / rings)) * cos(phi * Maths::Tau / rings);
					Vertices[(size_t)a * 3 + 1] = (outerRadius + innerRadius * cos(theta * Maths::Tau / rings)) * sin(phi * Maths::Tau / rings);
					Vertices[(size_t)a * 3 + 2] = depth * sin(theta * Maths::Tau / rings);

					Normals[(size_t)a * 3 + 0] = cos(theta * Maths::Tau / rings) * cos(phi * Maths::Tau / rings);
					Normals[(size_t)a * 3 + 1] = cos(theta * Maths::Tau / rings) * sin(phi * Maths::Tau / rings);
					Normals[(size_t)a * 3 + 2] = sin(theta * Maths::Tau / rings);

					//	TextureCoordinates.SetObjectValue(acos(currentVertex.y / innerRadius) / (2*phi), a * 2);
					//	TextureCoordinates.SetObjectValue(acos(currentVertex.x/(innerRadius + outerRadius*cos(theta*TwoPI / rings)))*(2*phi), a * 2 + 1);
					TextureCoordinates[(size_t)a * 2 + 0] = (i + k) / (float)sides;
					TextureCoordinates[(size_t)a * 2 + 1] = phi / (float)rings;

					a++;
				}
			}
		}

		j = (unsigned int)(indiciesCount / 6);
		for (i = 0; i < j; i++)
		{
			Indicies[(size_t)i * 6 + 0] = i * 2;
			Indicies[(size_t)i * 6 + 1] = i * 2 + 1;
			Indicies[(size_t)i * 6 + 2] = i * 2 + 2;

			Indicies[(size_t)i * 6 + 3] = i * 2 + 2;
			Indicies[(size_t)i * 6 + 4] = i * 2 + 1;
			Indicies[(size_t)i * 6 + 5] = i * 2 + 3;
		}

		return torusMesh;
	}
	//
	std::vector<unsigned int> Geometry::DetectSphereWrappedUVCoords(const std::vector<TriangleI>& triangles, const std::vector<glm::vec2>& texCoords, unsigned int hintCount)
	{
		std::vector<unsigned int> returnIndicies;
		returnIndicies.reserve(hintCount > 0 ? hintCount : 100);

		for (unsigned int i = 0; i < triangles.size(); ++i)
		{
			const TriangleI& currentTri = triangles[i];
			glm::vec3 texA(texCoords[currentTri[0]], 0.0f);
			glm::vec3 texB(texCoords[currentTri[1]], 0.0f);
			glm::vec3 texC(texCoords[currentTri[2]], 0.0f);

			glm::vec3 normal = glm::cross(texB - texA, texC - texA);
			if (normal.z > 0.0f)
			{
				returnIndicies.push_back(i);
			}
		}

		returnIndicies.shrink_to_fit();

		return move(returnIndicies);
	}
	//
	void Geometry::FixSphereWrappedUV(const std::vector<unsigned int>& wrapped, std::vector<glm::vec3>& verticies, std::vector<TriangleI>& triangles, std::vector<glm::vec2>& texCoords, unsigned int hintCount)
	{
		unsigned int verticieIndex = verticies.size() - 1;
		std::unordered_map<unsigned int, unsigned int> visited;

		if (hintCount > 0)
		{
			verticies.reserve(hintCount);
			texCoords.reserve(hintCount);
		}

		for (const unsigned int i : wrapped)
		{
			unsigned int a = triangles[i][0];
			unsigned int b = triangles[i][1];
			unsigned int c = triangles[i][2];

			if (texCoords[a].x < 0.25f)
			{
				unsigned int tempA = a;
				if (visited.find(a) == visited.end())
				{
					glm::vec3 dVertex = verticies[a];
					verticies.push_back(dVertex);
					glm::vec2 dTexCoord = texCoords[a];
					texCoords.emplace_back(dTexCoord.x + 1.0f, dTexCoord.y);
					verticieIndex++;
					visited.insert({ a, verticieIndex });
					tempA = verticieIndex;
				}
				else
				{
					tempA = visited[a];
				}
				a = tempA;
			}

			if (texCoords[b].x < 0.25f)
			{
				unsigned int tempB = b;
				if (visited.find(b) == visited.end())
				{
					glm::vec3 dVertex = verticies[b];
					verticies.push_back(dVertex);
					glm::vec2 dTexCoord = texCoords[b];
					texCoords.emplace_back(dTexCoord.x + 1.0f, dTexCoord.y);
					verticieIndex++;
					visited.insert({ b, verticieIndex });
					tempB = verticieIndex;
				}
				else
				{
					tempB = visited[b];
				}
				b = tempB;
			}

			if (texCoords[c].x < 0.25f)
			{
				unsigned int tempC = c;
				if (visited.find(b) == visited.end())
				{
					glm::vec3 dVertex = verticies[c];
					verticies.push_back(dVertex);
					glm::vec2 dTexCoord = texCoords[c];
					texCoords.emplace_back(dTexCoord.x + 1.0f, dTexCoord.y);
					verticieIndex++;
					visited.insert({ c, verticieIndex });
					tempC = verticieIndex;
				}
				else
				{
					tempC = visited[c];
				}
				c = tempC;
			}

			triangles[i].x = a;
			triangles[i].y = b;
			triangles[i].z = c;
		}

		verticies.shrink_to_fit();
		triangles.shrink_to_fit();
		texCoords.shrink_to_fit();
	}
	//
	void Geometry::FixPoleVerticies(std::vector<glm::vec3>& verticies, std::vector<glm::vec2>& texCoords, std::vector<TriangleI>& triangles)
	{
		unsigned int northPole = 0, southPole = 0;
		unsigned int i = 0;

		for (; i < verticies.size(); i++)
		{
			if (verticies[i].y >= 0.9995f) northPole = i;
			else if (verticies[i].y <= -0.9995f) southPole = i;
		}

		unsigned int verticeIndex = (unsigned int)verticies.size() - 1;

		auto FixVertex = [&triangles, &texCoords, &verticies, &verticeIndex](unsigned int i, unsigned int poleIndex, unsigned int b, unsigned int c) -> void
			{
				const TriangleI& currentTri = triangles[i];
				const glm::vec2& B = texCoords[currentTri[b]];
				const glm::vec2& C = texCoords[currentTri[c]];
				verticies.push_back(verticies[currentTri[poleIndex]]);
				//texCoords.emplace_back((B.x + C.x) * 0.5f, texCoords[currentTri[poleIndex]].y);
				texCoords.emplace_back(1, 1);
				verticeIndex++;
				triangles[i][poleIndex] = verticeIndex;
			};

		for (i = 0; i < triangles.size(); i++)
		{
			if (triangles[i].x == northPole)
			{
				FixVertex(i, 0, 1, 2);
			}
			else if (triangles[i].y == northPole)
			{
				FixVertex(i, 1, 0, 2);
			}
			else if (triangles[i].z == northPole)
			{
				FixVertex(i, 2, 0, 1);
			}

			if (triangles[i].x == southPole)
			{
				FixVertex(i, 0, 1, 2);
			}
			else if (triangles[i].y == southPole)
			{
				FixVertex(i, 1, 0, 2);
			}
			else if (triangles[i].z == southPole)
			{
				FixVertex(i, 2, 0, 1);
			}
		}

	}
	/*
	private void FixSharedPoleVertices(Mesh sphereMesh)
{
	Vertex north = sphereMesh.Vertices.First((v) => v.Position.Y == 1);
	Vertex south = sphereMesh.Vertices.First((v) => v.Position.Y == -1);
	uint northIndex = (uint)Array.IndexOf(sphereMesh.Vertices, north);
	uint southIndex = (uint)Array.IndexOf(sphereMesh.Vertices, south);
	List<Vertex> vertices = new List<Vertex>(sphereMesh.Vertices);
	uint verticeIndex = (uint)sphereMesh.Vertices.Length - 1;
	for (int i = 0; i < sphereMesh.Triangles.Length; ++i)
	{
		if (sphereMesh.Triangles[i].A == northIndex)
		{
			Vertex A = sphereMesh.Vertices[sphereMesh.Triangles[i].A];
			Vertex B = sphereMesh.Vertices[sphereMesh.Triangles[i].B];
			Vertex C = sphereMesh.Vertices[sphereMesh.Triangles[i].C];
			Vertex newNorth = north;
			newNorth.TexCoord.X = (B.TexCoord.X + C.TexCoord.X) / 2;
			verticeIndex++;
			vertices.Add(newNorth);
			sphereMesh.Triangles[i].A = verticeIndex;
		}
		else if (sphereMesh.Triangles[i].A == southIndex)
		{
			Vertex A = sphereMesh.Vertices[sphereMesh.Triangles[i].A];
			Vertex B = sphereMesh.Vertices[sphereMesh.Triangles[i].B];
			Vertex C = sphereMesh.Vertices[sphereMesh.Triangles[i].C];
			Vertex newSouth = south;
			newSouth.TexCoord.X = (B.TexCoord.X + C.TexCoord.X) / 2;
			verticeIndex++;
			vertices.Add(newSouth);
			sphereMesh.Triangles[i].A = verticeIndex;
		}
	}
	sphereMesh.Vertices = vertices.ToArray();
}
	*/
	//
#pragma endregion
//
}
