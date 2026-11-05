#version 450

// Instanced atom spheres. Binding 0 holds the unit sphere's positions; binding 3 holds one
// element per atom: (scale; float) (translation; vec3) (colour; vec3)
layout(location = 0) in vec3 aVertex;
layout(location = 1) in float aAtomScale;
layout(location = 2) in vec3 aAtomTranslation;
layout(location = 3) in vec3 aAtomColour;

layout(push_constant) uniform PushConstants
{
    mat4 uPVMatrix;
    vec4 uLightPos; // xyz used
    vec4 uParams;   // x = opacity
} pc;

layout(location = 0) out vec3 vNormal;
layout(location = 1) out vec3 vAtomColour;

void main()
{
    vNormal = normalize(aVertex);
    vAtomColour = aAtomColour;
    gl_Position = pc.uPVMatrix * vec4(aVertex * aAtomScale + aAtomTranslation, 1.0);
}
