#version 330 core

layout (location = 0) in uint v_type;
layout (location = 1) in vec3 v_pos;
layout (location = 2) in vec3 v_normal;

flat out uint obj_type;
out vec3 obj_pos;
out vec3 obj_normal;
out vec3 obj_color;

uniform vec3 g_color;
uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

void main()
{
    obj_type = v_type;
    vec4 pos_tmp = model_mat * vec4(v_pos, 1.0);
    obj_pos = pos_tmp.xyz;
    vec4 normal_tmp = model_mat * vec4(v_normal, 0.0);
    obj_normal = normal_tmp.xyz;
    obj_color = g_color;
    
    gl_Position = proj_mat * view_mat * model_mat * vec4(v_pos, 1.0f);
}
