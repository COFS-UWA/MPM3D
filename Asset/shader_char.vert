#version 330 core

layout (location = 0) in vec2 pt_pos;
layout (location = 1) in vec2 pt_tex;

out vec2 obj_tex;

uniform mat4 view_mat;

void main()
{
    obj_tex = pt_tex;
    gl_Position = view_mat * vec4(pt_pos, 0.0f, 1.0f);
}
