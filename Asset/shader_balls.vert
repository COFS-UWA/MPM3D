#version 330 core

layout (location = 0) in vec3 ball_pt_pos; // ball points position
layout (location = 1) in uint pt_type;    // point rendering method type
layout (location = 2) in vec3 pt_pos;     // point position
layout (location = 3) in float pt_radius; // point radius
layout (location = 4) in float pt_value;  // point value

out vec3 obj_color;
out vec3 obj_normal;
flat out uint obj_type;
out vec3 obj_pos;

uniform vec3 g_color;

uniform float color_value_lower;
uniform float color_value_range;
uniform sampler1D color_map;

uniform mat4 view_mat;
uniform mat4 proj_mat;

void main()
{
    if (pt_type == 0u)
    {
        obj_color = g_color;
    }
    else if (pt_type == 1u)
    {
        float tex_coord = (pt_value - color_value_lower) / color_value_range;
        obj_color = texture(color_map, tex_coord).xyz;
    }
    else // default
    {
        obj_color = vec3(0.0f, 0.0f, 0.0f);
    }
    
    obj_normal = ball_pt_pos;
    obj_type = pt_type;

    obj_pos = pt_pos + ball_pt_pos * pt_radius;
    gl_Position = proj_mat * view_mat * vec4(obj_pos, 1.0f);
}
