#include "ModelViewer_pcp.h"

#include "UniformColourMap.h"

void UniformColourMap::clear()
{
	if (colours)
	{
		delete[] colours;
		colours = nullptr;
	}
	colour_num = 0;
}

UniformColourMap::UniformColourMap() :
	colours(nullptr), colour_num(0),
	lb_colour(0.3f, 0.3f, 0.3f), // grey
	ub_colour(1.0f, 1.0f, 1.0f), // white
	lb_value(0.0f), ub_value(0.0f)
{

}

UniformColourMap::~UniformColourMap() { clear(); }

int UniformColourMap::set_colour(
	unsigned char colour_datas[][3],
	size_t colour_data_num,
	bool apply_ob_Colour
	)
{
	clear();

	if (!colour_datas || !colour_data_num)
		return -1;

	colour_num = colour_data_num;
	colours = new Colour[colour_num];
	for (size_t c_id = 0; c_id < colour_num; ++c_id)
	{
		Colour& c = colours[c_id];
		c.r = float(colour_datas[c_id][0]) / 255.0f;
		c.g = float(colour_datas[c_id][1]) / 255.0f;
		c.b = float(colour_datas[c_id][2]) / 255.0f;
	}

	if (apply_ob_Colour)
	{
		lb_colour = colours[0];
		ub_colour = colours[colour_num-1];
	}

	inv_len = (lb_value - ub_value) / float(colour_num-1);

	return 0;
}

unsigned char* UniformColourMap::gen_1Dtexture(
	size_t inv_resolution,
	size_t& texture_size
	)
{
	if (inv_resolution == 0)
	{
		texture_size = 0;
		return nullptr;
	}

	texture_size = (colour_num - 1) * inv_resolution + 1;
	unsigned char* texture = new unsigned char[texture_size*3];
	// init texture content
	size_t tex_num = colour_num - 1;
	struct Colour_ub
	{
		unsigned char r, g, b;
	} *cur_tex;
	cur_tex = reinterpret_cast<Colour_ub*>(texture);
	float w1, w2;
	float w_inv = 1.0f / float(inv_resolution);
	for (size_t t_id = 0; t_id < tex_num; ++t_id)
	{
		Colour& c1 = colours[t_id];
		Colour& c2 = colours[t_id+1];
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
	cur_tex->r = unsigned char(colours[tex_num].r * 255.0f);
	cur_tex->g = unsigned char(colours[tex_num].g * 255.0f);
	cur_tex->b = unsigned char(colours[tex_num].b * 255.0f);
	return texture;
}
