#version 410

#pragma remove
#include "GLSL Language Specification.h.glsl"
#pragma endremove

#begin(vShader)

layout(location = 0) in vec3 aVertex;
//
uniform mat4 uMPVMatrix;
//
out vec3 vNormal;

void main()
{
    vNormal = normalize(aVertex);
    gl_Position = uMPVMatrix * vec4(aVertex, 1.0f);
}

#end(vShader)
#begin(fShader)

in vec3 vNormal;
in float vSign;
//
out vec4 fragColour;

uniform vec3 uLightPos = vec3(10.0f, 10.0f, 10.0f);
uniform vec3 uColour = vec3(1.0f, 1.0f, 1.0f);
uniform float uOpacity = 1.0f;

void main()
{
    vec3 unitLightVector = normalize(uLightPos);
    float diffuseAmount = clamp(dot(unitLightVector, vNormal), 0.1f, 1.0f);

    fragColour = vec4(uColour * diffuseAmount, uOpacity);
}
#end(fShader)