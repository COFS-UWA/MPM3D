#version 330 core

layout (location = 0) in vec3 vertex_pos;
layout (location = 1) in vec3 vertex_color;

out vec3 mid_color;

uniform mat4 view_mat;
uniform mat4 proj_mat;

void main()
{
    mid_color = vertex_color;
    gl_Position = proj_mat * view_mat * vec4(vertex_pos, 1.0f);
}
