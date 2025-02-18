#include "AtomDraw.h"
#include<glm/glm.hpp>

//
//
namespace SCESS::Viewer {
    //
    //
    std::unordered_map<std::string_view, AtomStyle> AtomStyle::sAtomStyle =
    {
        {"H", {{1.0f, 1.0f, 1.0f}, 1.0f}},
        {"C", {{0.2f, 0.2f, 0.2f}, 1.5f}},
        {"Li", {{0.8f, 0.9f, 1.0f}, 3.15f}},
        {"N", {{0.0f, 0.5f, 0.5f}, 1.8f}},
        {"O", {{1.0f, 0.0f, 0.0f}, 1.1f}},
        {"B", {{0.100f, 0.714f, 0.757f}, 2.4f}}
    };
    //
    glm::mat3 CalculateBondRotMatrix(const glm::vec3 &bondDirNorm, const glm::vec3 &startingDirNorm) {
        glm::vec3 v = glm::cross(bondDirNorm, startingDirNorm);
        float c = glm::dot(bondDirNorm, startingDirNorm);

        glm::mat3 vSkew = {
            0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0
        };

        return glm::mat3(1.0f) + vSkew + vSkew * vSkew / (1 + c);
    }
    //
    //
}
