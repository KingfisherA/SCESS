#version 450

// Isosurface triangles are generated in world space, so there is no model transform
layout(location = 0) in vec3 aVertex;

layout(push_constant) uniform PushConstants
{
    mat4 uPVMatrix;
    vec4 uLightPos; // xyz used
    vec4 uColour;   // rgb colour, a = opacity
} pc;

layout(location = 0) out vec3 vNormal;

void main()
{
    vNormal = normalize(aVertex);
    gl_Position = pc.uPVMatrix * vec4(aVertex, 1.0);
}
