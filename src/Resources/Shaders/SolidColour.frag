#version 450

layout(push_constant) uniform PushConstants
{
    mat4 uMVPMatrix; // projection * view * model
    vec4 uColour;    // rgb colour, a = opacity
} pc;

layout(location = 0) out vec4 fragColour;

void main()
{
    fragColour = pc.uColour;
}
