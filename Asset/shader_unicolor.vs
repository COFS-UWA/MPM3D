#version 330 core

layout(location = 0) in vec3 in_pos;

out vec3 mid_color;

uniform vec3 color;
uniform mat4 view_mat;
uniform mat4 proj_mat;

void main()
{
	mid_color = color;
    gl_Position = proj_mat * view_mat * vec4(in_pos, 1.0f);
}