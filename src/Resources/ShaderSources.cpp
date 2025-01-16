#include<pch.h>
#include "AssetFolder.h"
#include "ShaderSources.h"

namespace SCESS::Res {
    static const char *MOLECULE_COLOUR_DIFFUSE_GLSL = R"(#version 410

#pragma remove
#include "GLSL Language Specification.h.glsl"
#pragma endremove

#begin(vShader)

layout(location = 0) in vec3 aVertex;
layout(location = 1) in float aAtomScale;
layout(location = 2) in vec3 aAtomTranslation;
layout(location = 3) in vec3 aAtomColour;
//
uniform mat4 uMPVMatrix;
//
out vec3 vNormal;
out vec3 vAtomColour;

void main()
{
    vNormal = normalize(aVertex);
    vAtomColour = aAtomColour;
    gl_Position = uMPVMatrix * vec4(aVertex * aAtomScale + aAtomTranslation, 1.0f);
}

#end(vShader)
#begin(fShader)

in vec3 vNormal;
in vec3 vAtomColour;
//
out vec4 fragColour;

uniform vec3 uLightPos = vec3(10.0f, 10.0f, 10.0f);
uniform float uOpacity = 1.0f;

void main()
{
    vec3 unitLightVector = normalize(uLightPos);
    float diffuseAmount = clamp(dot(unitLightVector, vNormal), 0.1f, 1.0f);

    fragColour = vec4(vAtomColour * diffuseAmount, uOpacity);
}
#end(fShader))";

    static const char *SOLID_COLOUR_GLSL = R"(#version 410

#pragma remove
#include "GLSL Language Specification.h.glsl"
#pragma endremove

#begin(vShader)

layout(location = 0) in vec3 aVertex;
//
uniform mat4 uTransformMatrix;
uniform mat4 uPVMatrix;
//

void main()
{
    gl_Position = uPVMatrix * uTransformMatrix * vec4(aVertex, 1.0f);
}

#end(vShader)
#begin(fShader)

//
out vec4 fragColour;

uniform vec3 uColour = vec3(1.0f, 1.0f, 1.0f);
uniform float uOpacity = 1.0f;

void main()
{
    fragColour = vec4(uColour, uOpacity);
}
#end(fShader))";

    const char *ShaderSources::GetMoleculeColourDiffuseGLSL() {
        return MOLECULE_COLOUR_DIFFUSE_GLSL;
    }

    const char *ShaderSources::GetSolidColourGLSL() {
        return SOLID_COLOUR_GLSL;
    }
}
