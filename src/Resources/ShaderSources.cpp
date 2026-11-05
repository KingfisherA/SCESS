#include "ShaderSources.h"

//
//
namespace SCESS::Res::ShaderSources {
    //
    //
    namespace {
        const uint32_t MOLECULE_DIFFUSE_VERT[] =
#include "MoleculeDiffuse.vert.inc"
        ;
        const uint32_t MOLECULE_DIFFUSE_FRAG[] =
#include "MoleculeDiffuse.frag.inc"
        ;
        const uint32_t ISOSURFACE_DIFFUSE_VERT[] =
#include "IsosurfaceDiffuse.vert.inc"
        ;
        const uint32_t ISOSURFACE_DIFFUSE_FRAG[] =
#include "IsosurfaceDiffuse.frag.inc"
        ;
        const uint32_t SOLID_COLOUR_VERT[] =
#include "SolidColour.vert.inc"
        ;
        const uint32_t SOLID_COLOUR_FRAG[] =
#include "SolidColour.frag.inc"
        ;
    }
    //
    std::span<const uint32_t> GetMoleculeDiffuseVert() { return MOLECULE_DIFFUSE_VERT; }
    std::span<const uint32_t> GetMoleculeDiffuseFrag() { return MOLECULE_DIFFUSE_FRAG; }
    std::span<const uint32_t> GetIsosurfaceDiffuseVert() { return ISOSURFACE_DIFFUSE_VERT; }
    std::span<const uint32_t> GetIsosurfaceDiffuseFrag() { return ISOSURFACE_DIFFUSE_FRAG; }
    std::span<const uint32_t> GetSolidColourVert() { return SOLID_COLOUR_VERT; }
    std::span<const uint32_t> GetSolidColourFrag() { return SOLID_COLOUR_FRAG; }
    //
    //
}
