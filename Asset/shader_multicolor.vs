#version 330 core

layout (location = 0) in vec3 vertex_pos;
layout (location = 1) in vec3 vertex_color;

out vec4 mid_color;

uniform mat4 mv_mat;
uniform mat4 proj_mat;

void main()
{
    gl_Position = proj_mat * mv_mat * vec4(vertex_pos, 1.0f);
    mid_color = vec4(vertex_color, 1.0f);
}
