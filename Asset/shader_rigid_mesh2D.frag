#version 330 core

flat in int obj_type;
in vec3 obj_color;

out vec4 frag_color;

void main()
{
	frag_color = vec4(obj_color, 1.0f);
}
