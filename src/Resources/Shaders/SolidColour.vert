#version 450

layout(location = 0) in vec3 aVertex;

layout(push_constant) uniform PushConstants
{
    mat4 uMVPMatrix; // projection * view * model
    vec4 uColour;    // rgb colour, a = opacity
} pc;

void main()
{
    gl_Position = pc.uMVPMatrix * vec4(aVertex, 1.0);
}
