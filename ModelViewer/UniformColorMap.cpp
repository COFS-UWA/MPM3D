#include "ModelViewer_pcp.h"

#include "UniformColorMap.h"

void UniformColorMap::clear()
{
	if (colors)
	{
		delete[] colors;
		colors = nullptr;
	}
	color_num = 0;
}

UniformColorMap::UniformColorMap() :
	colors(nullptr), color_num(0),
	lb_color(0.3f, 0.3f, 0.3f), // grey
	ub_color(1.0f, 1.0f, 1.0f), // white
	lb_value(0.0f), ub_value(0.0f)
{

}

UniformColorMap::~UniformColorMap() { clear(); }

int UniformColorMap::set_color(
	unsigned char color_datas[][3],
	size_t color_data_num,
	bool apply_ob_color
	)
{
	clear();

	if (!color_datas || !color_data_num)
		return -1;

	color_num = color_data_num;
	colors = new Color[color_num];
	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		Color& c = colors[c_id];
		c.r = float(color_datas[c_id][0]) / 255.0f;
		c.g = float(color_datas[c_id][1]) / 255.0f;
		c.b = float(color_datas[c_id][2]) / 255.0f;
	}

	if (apply_ob_color)
	{
		lb_color = colors[0];
		ub_color = colors[color_num-1];
	}

	inv_len = (lb_value - ub_value) / float(color_num-1);

	return 0;
}

unsigned char* UniformColorMap::gen_1Dtexture(
	size_t inv_resolution,
	size_t& texture_size
	)
{
	if (inv_resolution == 0)
	{
		texture_size = 0;
		return nullptr;
	}

	texture_size = (color_num - 1) * inv_resolution + 1;
	unsigned char* texture = new unsigned char[texture_size*3];
	// init texture content
	size_t tex_num = color_num - 1;
	struct color_ub
	{
		unsigned char r, g, b;
	} *cur_tex;
	cur_tex = reinterpret_cast<color_ub*>(texture);
	float w1, w2;
	float w_inv = 1.0f / float(inv_resolution);
	for (size_t t_id = 0; t_id < tex_num; ++t_id)
	{
		Color &c1 = colors[t_id];
		Color &c2 = colors[t_id+1];
		w2 = 0.0f;
		for (size_t r_id = 0; r_id < inv_resolution; ++r_id)
		{
			w1 = 1.0f - w2;
			cur_tex->r = unsigned char((w1 * c1.r + w2 * c2.r) * 255.0f);
			cur_tex->g = unsigned char((w1 * c1.g + w2 * c2.g) * 255.0f);
			cur_tex->b = unsigned char((w1 * c1.b + w2 * c2.b) * 255.0f);
			++cur_tex;
			w2 += w_inv;
		}
	}
	cur_tex->r = unsigned char(colors[tex_num].r * 255.0f);
	cur_tex->g = unsigned char(colors[tex_num].g * 255.0f);
	cur_tex->b = unsigned char(colors[tex_num].b * 255.0f);
	return texture;
}
