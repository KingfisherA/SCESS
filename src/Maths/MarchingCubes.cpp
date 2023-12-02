//Adapted from https://github.com/nihaljn/marching-cubes

#include<pch.h>
#include "MarchingCubes.h"

//
//
namespace SCESS
{
    int MarchingCubes::calculate_cube_index(GridCell& cell, float isovalue)
    {
        int cubeIndex = 0;
        for (int i = 0; i < 8; i++)
            if (cell.value[i] < isovalue) cubeIndex |= (1 << i);
        return cubeIndex;
    }

    std::vector<glm::vec3> MarchingCubes::get_intersection_coordinates(GridCell& cell, float isovalue)
    {
        std::vector<glm::vec3> intersections(12);

        int cubeIndex = calculate_cube_index(cell, isovalue);
        int intersectionsKey = edgeTable[cubeIndex];

        int idx = 0;
        while (intersectionsKey)
        {
            if (intersectionsKey & 1)
            {
                int v1 = edgeToVertices[idx].first, v2 = edgeToVertices[idx].second;
                glm::vec3 intersectionPoint = interpolate(cell.vertex[v1], cell.value[v1],
                    cell.vertex[v2], cell.value[v2], isovalue);
                intersections[idx] = intersectionPoint;
            }
            idx++;
            intersectionsKey >>= 1;
        }


        return intersections;
    }

    glm::vec3 MarchingCubes::interpolate(glm::vec3& v1, float val1, glm::vec3& v2, float val2, float isovalue)
    {
        glm::vec3 interpolated;
        float mu = (isovalue - val1) / (val2 - val1);

        interpolated.x = mu * (v2.x - v1.x) + v1.x;
        interpolated.y = mu * (v2.y - v1.y) + v1.y;
        interpolated.z = mu * (v2.z - v1.z) + v1.z;

        return interpolated;
    }

    std::vector<Triangle> MarchingCubes::GetTriangles(std::vector<glm::vec3>& intersections, int cubeIndex)
    {
        std::vector<Triangle> triangles;
        for (int i = 0; triangleTable[cubeIndex][i] != -1; i += 3)
        {
            Triangle triangle;
            for (int j = 0; j < 3; j++)
                triangle[j] = intersections[triangleTable[cubeIndex][i + j]];
            triangles.push_back(triangle);
        }

        return triangles;
    }

    void MarchingCubes::PrintTriangles(const std::vector<Triangle>& triangles)
    {
        for (int i = 0; i < (int)triangles.size(); i++)
        {
            for (int j = 0; j < 3; j++)
                std::cout << triangles[i][j].x << ",\t" << triangles[i][j].y << ",\t" << triangles[i][j].z << "\n";
            std::cout << "\n";
        }
    }

    std::vector<Triangle> MarchingCubes::TriangulateCell(GridCell& cell, float isovalue)
    {
        int cubeIndex = calculate_cube_index(cell, isovalue);
        std::vector<glm::vec3> intersections = get_intersection_coordinates(cell, isovalue);
        std::vector<Triangle> triangles = GetTriangles(intersections, cubeIndex);

        // int flag = 0;
        // for (auto& triangle: triangles)
        //     for (int i = 0; i < 3; i++)
        //         for (int j = i + 1; j < 3; j++)
        //             if (triangle[i].x == triangle[j].x && triangle[i].y == triangle[j].y && triangle[i].z == triangle[j].z)
        //                 flag++;

        // if (flag)
        // {
        //     print_triangles(triangles);
        //     for (int i = 0; i < 8; i++)
        //         std::cout << cell.vertex[i].x << ", " << cell.vertex[i].y << ", " << cell.vertex[i].z << ": " << cell.value[i] << "\n";
        //     std::cout << "\n";
        // }

        return triangles;
    }

    void MarchingCubes::TriangulateField(
        const std::vector<std::vector<std::vector<float>>>& scalarFunction,
        float isovalue, std::vector<Triangle>* trianglesOut)
    {

        if (trianglesOut == nullptr)
        {
            throw std::exception("trianglesOut cannot be nullptr");
        }

        size_t max = scalarFunction.size();
        
        for (size_t i = 0; i + 1 < max; i++)
        {
            for (size_t j = 0; j + 1 < max; j++)
            {
                for (size_t k = 0; k + 1 < max; k++)
                {
                    float x = i, y = j, z = k;
                    // cell ordered according to convention in referenced website
                    GridCell cell =
                    {
                        {
                            {x, y, z}, {x + 1.0f, y, z},
                            {x + 1.0f, y, z + 1.0f}, {x, y, z + 1.0f},
                            {x, y + 1.0f, z}, {x + 1.0f, y + 1.0f, z},
                            {x + 1.0f, y + 1.0f, z + 1.0f}, {x, y + 1.0f, z + 1.0f}
                        },
                        {
                            scalarFunction[i][j][k], scalarFunction[i + 1][j][k],
                            scalarFunction[i + 1][j][k + 1], scalarFunction[i][j][k + 1],
                            scalarFunction[i][j + 1][k], scalarFunction[i + 1][j + 1][k],
                            scalarFunction[i + 1][j + 1][k + 1], scalarFunction[i][j + 1][k + 1]
                        }
                    };

                    std::vector<Triangle> cellTriangles = TriangulateCell(cell, isovalue);

                    trianglesOut->reserve(cellTriangles.size());
                    for (int i = 0; i < (int)cellTriangles.size(); i++)
                    {
                        trianglesOut->push_back(cellTriangles[i]);
                    }
                }
            }
        }
    }
}