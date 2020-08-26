#version 330 core

layout (location = 0) in uint v_type;
layout (location = 1) in vec3 v_pos;
layout (location = 2) in vec3 v_normal;
layout (location = 3) in vec3 v_color;

flat out uint obj_type;
out vec3 obj_pos;
out vec3 obj_normal;
out vec3 obj_color;

uniform vec3 g_color;
uniform mat4 view_mat;
uniform mat4 proj_mat;

void main()
{
    obj_type = v_type;
    obj_pos = v_pos;
    obj_normal = v_normal;

    if (v_type == 0u)
    {
        obj_color = g_color;
    }
    else if (v_type == 1u)
    {
        obj_color = v_color;
    }
    else
    {
        obj_color = vec3(0.0f, 0.0f, 0.0f);
    }

    gl_Position = proj_mat * view_mat * vec4(v_pos, 1.0f);
}
