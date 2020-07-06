#version 330 core

layout (location = 0) in vec2 circle_pt_pos; // circle point position
layout (location = 1) in uint pt_type;    // point rendering method type
layout (location = 2) in vec2 pt_pos;     // point position
layout (location = 3) in float pt_radius; // point radius
layout (location = 4) in float pt_value;  // point value

out vec3 obj_color;
flat out uint obj_type;

uniform vec3 g_color;

uniform float color_value_lower;
uniform float color_value_range;
uniform sampler1D color_map;

uniform mat4 view_mat;

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
    
    obj_type = pt_type;
    
    vec2 obj_pos = pt_pos + circle_pt_pos * pt_radius;
    gl_Position = view_mat * vec4(obj_pos, 0.0f, 1.0f);
}
