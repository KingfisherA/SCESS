#pragma once

#include <cstdint>
#include <span>

//
//
namespace SCESS::Res::ShaderSources {
    //
    // SPIR-V for the viewer's shaders, compiled at build time from src/Resources/Shaders
    // (see dahlia_compile_shaders in cmake/ViewerSCESS.cmake)
    std::span<const uint32_t> GetMoleculeDiffuseVert();
    std::span<const uint32_t> GetMoleculeDiffuseFrag();
    std::span<const uint32_t> GetIsosurfaceDiffuseVert();
    std::span<const uint32_t> GetIsosurfaceDiffuseFrag();
    std::span<const uint32_t> GetSolidColourVert();
    std::span<const uint32_t> GetSolidColourFrag();
    //
    //
}
