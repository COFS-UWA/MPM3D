#version 330 core

out vec4 frag_color;

flat in uint obj_type;
in vec3 obj_pos;
in vec3 obj_normal;
in vec3 obj_color;

uniform vec3 view_pos;

uniform float fog_coef;
uniform vec3 fog_color;

// phong model parameters
uniform vec3 light_pos;
uniform vec3 light_color;

uniform float amb_coef;
uniform float diff_coef;
uniform float spec_coef;
uniform float spec_shininess;

vec3 phong_lighting()
{
    // ambient
    //float amb_coef = 0.1;
    vec3 ambient = amb_coef * light_color;

    // diffuse 
    vec3 norm = normalize(obj_normal);
    vec3 light_dir = normalize(light_pos - obj_pos);
    float diff_mid = max(dot(norm, light_dir), 0.0);
    vec3 diffuse = diff_coef * diff_mid * light_color;
    
    // specular
    //float spec_coef = 0.5f;
    //float spec_shininess = 32.0f;
    vec3 view_dir = normalize(view_pos - obj_pos);
    vec3 reflect_dir = reflect(-light_dir, norm);
    float spec_mid = pow(max(dot(view_dir, reflect_dir), 0.0), spec_shininess);
    vec3 specular = spec_coef * spec_mid * light_color;
    
    return (diffuse + ambient + specular);
}

void main()
{
    vec3 result = vec3(0.0, 0.0, 0.0);
    if (obj_type == 1u)
    {
        result = obj_color * phong_lighting();
    }
    else if (obj_type == 2u)
    {
        result = obj_color;
    }
    // frog
    float dist = distance(obj_pos, view_pos);
    float fog_mid = exp(-fog_coef * dist);
    result = fog_mid * result + (1.0f - fog_mid) * fog_color;
    frag_color = vec4(result, 1.0f);
} 
