#version 330 core

in vec2 obj_tex;

out vec4 frag_color;

uniform vec3 g_color;
uniform sampler2D char_texture;

void main()
{
    float frag_alpha = texture(char_texture, obj_tex).r;
    frag_color = vec4(g_color, frag_alpha);
}
