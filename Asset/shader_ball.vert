#version 330 core

layout (location = 0) in vec3 ball_pt_pos; // ball points
layout (location = 1) in vec3 pcl_pos; // particle position
layout (location = 2) in float pcl_r;
layout (location = 3) in vec3 pcl_color;

out vec3 obj_pos;
out vec3 obj_color;
out vec3 obj_normal;

uniform mat4 view_mat;
uniform mat4 proj_mat;

void main()
{
    obj_pos = ball_pt_pos * pcl_r + pcl_pos;
    obj_color = pcl_color;
    obj_normal = ball_pt_pos;
    
    gl_Position = proj_mat * view_mat * vec4(obj_pos, 1.0f);
}
