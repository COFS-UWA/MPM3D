#version 330 core

layout (location = 0) in vec3 v_pos;
layout (location = 1) in vec3 v_color;
layout (location = 2) in vec3 v_normal;

out vec3 obj_pos;
out vec3 obj_color;
out vec3 obj_normal;

uniform mat4 view_mat;
uniform mat4 proj_mat;

void main()
{
    obj_pos = v_pos;
    obj_color = v_color;
    obj_normal = v_normal;
    
    gl_Position = proj_mat * view_mat * vec4(v_pos, 1.0f);
}
