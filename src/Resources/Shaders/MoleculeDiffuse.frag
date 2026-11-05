#version 450

layout(location = 0) in vec3 vNormal;
layout(location = 1) in vec3 vAtomColour;

layout(push_constant) uniform PushConstants
{
    mat4 uPVMatrix;
    vec4 uLightPos; // xyz used
    vec4 uParams;   // x = opacity
} pc;

layout(location = 0) out vec4 fragColour;

void main()
{
    vec3 unitLightVector = normalize(pc.uLightPos.xyz);
    float diffuseAmount = clamp(dot(unitLightVector, vNormal), 0.1, 1.0);

    fragColour = vec4(vAtomColour * diffuseAmount, pc.uParams.x);
}
