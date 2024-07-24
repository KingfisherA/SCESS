#version 410

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
#end(fShader)