#pragma once

#include "GraphicsEngine/GEng_API.h"
#include "GraphicsEngine/Utilities/Memory.h"

#include <array>
#include <vector>
#include<functional>
#include <glm/glm.hpp>

namespace GEng {
    template<typename FloatType, typename IndexType>
    GENG_API class Mesh {
    public:
        explicit Mesh(const uint8_t vertexSize): mVertexSize(vertexSize) {
        }
        Mesh(const uint8_t vertexSize, const size_t vArraySize, const size_t nArraySize, const size_t tArraySize,
             const size_t iArraySize)
            : mVertexSize(vertexSize), mVertexArray(vArraySize), mNormalsArray(nArraySize), mTexCoordsArray(tArraySize),
              mIndicesArray(iArraySize) {
        }
        Mesh(const Mesh &mesh): mVertexSize(mesh.mVertexSize), mVertexArray(mesh.mVertexArray),
                                mNormalsArray(mesh.mNormalsArray), mTexCoordsArray(mesh.mTexCoordsArray),
                                mIndicesArray(mesh.mIndicesArray) {
        }
        Mesh(Mesh &&mesh) noexcept: mVertexSize(mesh.mVertexSize), mVertexArray(std::move(mesh.mVertexArray)),
                                    mNormalsArray(std::move(mesh.mNormalsArray)),
                                    mTexCoordsArray(std::move(mesh.mTexCoordsArray)),
                                    mIndicesArray(std::move(mesh.mIndicesArray)) {
        }
        ~Mesh() = default;

        static Mesh RemoveIndexing(const Mesh &mesh) {
            if (!mesh.HasIndices())
                return mesh;

            Mesh meshCopy(mesh.mVertexSize, mesh.GetIndexIntCount() * mesh.mVertexSize,
                          mesh.GetNormalsFloatCount() == 0 ? 0 : mesh.GetIndexIntCount() * mesh.mVertexSize,
                          mesh.GetTexCoordFloatCount() == 0 ? 0 : mesh.GetIndexIntCount() * 2,
                          0);

            std::vector<FloatType> &vertices = meshCopy.GetVertexArray();
            std::vector<FloatType> &normals = meshCopy.GetNormalsArray();
            std::vector<FloatType> &texCoords = meshCopy.GetTexCoordsArray();

            size_t j = 0;
            for (size_t i = 0; i < mesh.GetIndexIntCount(); i++) {
                IndexType index = mesh.GetIndicesArray()[i];
                if (mesh.HasVertex()) {
                    for (uint8_t k = 0; k < mesh.mVertexSize; k++) {
                        vertices[j++] = mesh.GetVertexArray()[index * (size_t) mesh.mVertexSize + k];
                    }
                }
                if (mesh.HasNormals()) {
                    for (uint8_t k = 0; k < mesh.mVertexSize; k++) {
                        normals[j++] = mesh.GetNormalsArray()[index * (size_t) mesh.mVertexSize + k];
                    }
                }
                if (mesh.HasTexCoords()) {
                    texCoords[j++] = mesh.GetTexCoordsArray()[index * 2 + 0];
                    texCoords[j++] = mesh.GetTexCoordsArray()[index * 2 + 1];
                }
            }

            return meshCopy;
        }
        //
        std::vector<FloatType> InterleaveVertexData() { return {}; }
        //
        std::vector<FloatType> &GetVertexArray() { return mVertexArray; }
        const std::vector<FloatType> &GetVertexArray() const { return mVertexArray; }
        std::vector<FloatType> &GetNormalsArray() { return mNormalsArray; }
        const std::vector<FloatType> &GetNormalsArray() const { return mNormalsArray; }
        std::vector<FloatType> &GetTexCoordsArray() { return mTexCoordsArray; }
        const std::vector<FloatType> &GetTexCoordsArray() const { return mTexCoordsArray; }
        std::vector<IndexType> &GetIndicesArray() { return mIndicesArray; }
        const std::vector<IndexType> &GetIndicesArray() const { return mIndicesArray; }
        //
        [[nodiscard]] size_t GetVertexFloatCount() const { return mVertexArray.size(); }
        [[nodiscard]] size_t GetNormalsFloatCount() const { return mNormalsArray.size(); }
        [[nodiscard]] size_t GetTexCoordFloatCount() const { return mTexCoordsArray.size(); }
        [[nodiscard]] size_t GetIndexIntCount() const { return mIndicesArray.size(); }
        //
        [[nodiscard]] bool HasVertex() const { return !mVertexArray.empty(); }
        [[nodiscard]] bool HasNormals() const { return !mNormalsArray.empty(); }
        [[nodiscard]] bool HasTexCoords() const { return !mTexCoordsArray.empty(); }
        [[nodiscard]] bool HasIndices() const { return !mIndicesArray.empty(); }
        //
        [[nodiscard]] uint8_t GetVertexSize() const { return mVertexSize; }
        //
    private:
        std::vector<FloatType> mVertexArray, mNormalsArray, mTexCoordsArray;
        std::vector<IndexType> mIndicesArray;

        uint8_t mVertexSize;
    };


    using Mesh_F16 = Mesh<float, uint16_t>;
    using Mesh_D16 = Mesh<double, uint16_t>;
    using Mesh_F32 = Mesh<float, uint32_t>;
    using Mesh_D32 = Mesh<double, uint32_t>;

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
    constexpr std::array<uint32_t, 60> ICOSAHEDRON_TRIS = {
        0, 4, 1, 0, 9, 4, 9, 5, 4, 4, 5, 8, 4, 8, 1,
        8, 10, 1, 8, 3, 10, 5, 3, 8, 5, 2, 3, 2, 7, 3,
        7, 10, 3, 7, 6, 10, 7, 11, 6, 11, 0, 6, 0, 1, 6,
        6, 1, 10, 9, 0, 11, 9, 11, 2, 9, 2, 5, 7, 2, 11
    };

    //
    GENG_API class Geometry final {
    public:
        using TriangleI = glm::vec<3, uint32_t>;
        using Vertex = glm::vec3;
        using TexCoord = glm::vec2;
        //
        static glm::vec3 CalculateNormalCounterclockwise(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2);
        static glm::vec3 CalculateNormalClockwise(const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2);
        static void CalculateNormals(const std::vector<float> &vertices, const std::vector<uint32_t> &triangles,
                                     std::vector<float> &normals);
        //
        static Mesh_F32 CreateUVSphere(float radius, unsigned int lattitudes, unsigned int longitudes);
        static Mesh_F32 CreateQuadsphere(float radius, unsigned int subdivisions);
        //
        static Mesh_F32 CreateIcosahedron(float size);
        static Mesh_F32 CreateIcosphere(float radius, unsigned int subdivisions);
        //
        static Mesh_F32 CreateCylinder(float radius, float height, size_t segments);
        //
        static Mesh_F32 CreateCone(float height, float radius, size_t segments);
        static Mesh_F32 CreatePyramid(float size);
        static Mesh_F32 CreateOctahedron(float height, float size);
        static Mesh_F32 CreateTetrahedron(float size);
        static Mesh_F32 CreateDodecahedron(float size);
        //
    private:
        //
        static std::vector<uint32_t> DetectSphereWrappedUVCoords(const std::vector<uint32_t> &triangles,
                                                                 const std::vector<float> &texCoords,
                                                                 unsigned int hintCount = 0);
        static void FixSphereWrappedUV(const std::vector<uint32_t> &wrapped, std::vector<float> &verticies,
                                       std::vector<uint32_t> &triangles, std::vector<float> &texCoords,
                                       unsigned int hintCount = 0);
        //static void FixPoleVerticies(std::vector<float> &verticies, std::vector<float> &texCoords,
        //                             std::vector<uint32_t> &triangles);
        //
    public:
        //
        static Mesh_F32 CreateQuad();
        //
        static Mesh_F32 CreateSquircle(float radius, uint32_t roundness, uint32_t segments);
        //
        static Mesh_F32 CreateZSurface(const glm::vec3 scaling, const glm::vec2 axisStart, const glm::vec2 axisEnd,
                                       uint32_t resolution, const std::function<float(float, float)> &func);
        //
        static Mesh_F32 CreateGrid(float scale, unsigned int slices, unsigned int stacks);
        //
        static Mesh_F32 CreateCube(float size = 1.0f);
        //
        static Mesh_F32 CreateCircle(float radius, size_t segments);
        static Mesh_F32 CreateRings(float innerRadius, float outerRadius, unsigned int segments);
        static Mesh_F32 CreateHalfRings(float innerRadius, float outerRadius, unsigned int segments);
        //
        static Mesh_F32 CreateTorus(float innerRadius,
                                    float outerRadius,
                                    float depth,
                                    uint32_t sides,
                                    uint32_t rings,
                                    uint32_t twistyness = 1);
        //
    };
}
