#version 330 core

in vec3 obj_color;
flat in int obj_type;

out vec4 frag_color;

void main()
{
	frag_color = vec4(obj_color, 1.0f);
	//frag_color = vec4(1.0f, 1.0f, 1.0f, 1.0f);
}
