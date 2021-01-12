#include "ModelViewer_pcp.h"

#include "UniformColorMap_Mono.h"

UniformColorMap_Mono::UniformColorMap_Mono()
{
	unsigned char color_map[2][3];
	color_map[0][0] = unsigned char(255.0f);
	color_map[0][1] = unsigned char(0.8941f * 255.0f);
	color_map[0][2] = unsigned char(0.7098f * 255.0f);
	color_map[1][0] = color_map[0][0];
	color_map[1][1] = color_map[0][1];
	color_map[1][2] = color_map[0][2];
	UniformColorMap::set_color(color_map, 2);
}

UniformColorMap_Mono::~UniformColorMap_Mono() {}

void UniformColorMap_Mono::set_color(
	float r, float g, float b)
{
	unsigned char color_map[2][3];
	color_map[0][0] = unsigned char(r * 255.0f);
	color_map[0][1] = unsigned char(g * 255.0f);
	color_map[0][2] = unsigned char(b * 255.0f);
	color_map[1][0] = color_map[0][0];
	color_map[1][1] = color_map[0][1];
	color_map[1][2] = color_map[0][2];
	UniformColorMap::set_color(color_map, 2);
}
