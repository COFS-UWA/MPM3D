#version 330 core

layout (location = 0) in int v_type;
layout (location = 1) in vec2 v_pos;
layout (location = 2) in vec3 v_color;

flat out int obj_type;
out vec3 obj_color;

uniform vec3 g_color;
uniform mat4 model_mat;
uniform mat4 view_mat;

void main()
{    
    obj_type = v_type;

    if (v_type == 0)
    {
        obj_color = g_color;
    }
    else if (v_type == 1)
    {
        obj_color = v_color;
    }
    else
    {
        obj_color = vec3(0.0f, 0.0f, 0.0f);
    }

    gl_Position = view_mat * model_mat * vec4(v_pos, 0.0f, 1.0f);
}
