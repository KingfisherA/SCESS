#pragma once

#include<string_view>
#include<unordered_map>
#include<glm/fwd.hpp>
#include<glm/vec3.hpp>

namespace SCESS {
    struct AtomStyle {
        //RGB colour
        glm::vec3 mColour;
        //In Bohr radii
        float mRadius;

        static std::unordered_map<std::string_view, AtomStyle> sAtomStyle;

        inline static AtomStyle GetStyle(const std::string_view atom) {
            auto it = sAtomStyle.find(atom);
            if (it != sAtomStyle.end())
                return it->second;
            return {{0.1, 0.1, 0.1}, 1.0f};
        }
    };

    glm::mat3 CalculateBondRotMatrix(const glm::vec3 &bondDirNorm, const glm::vec3 &startingDirNorm);
}
