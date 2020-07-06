#ifndef __Uniform_Colour_Map_h__
#define __Uniform_Colour_Map_h__

#include <cmath>

class UniformColourMap
{
public:
	struct Colour
	{
		float r, g, b;
		Colour() {}
		Colour(float _r, float _g, float _b) :
			r(_r), g(_g), b(_b) {}
	};

protected:
	size_t colour_num;
	Colour *colours;
	Colour lb_colour, ub_colour;

	float lb_value, ub_value;
	float inv_len;

	void clear();

public:
	UniformColourMap();
	~UniformColourMap();

	int set_colour(
		unsigned char colour_datas[][3],
		size_t colour_data_num,
		bool apply_ob_colour = true
		);
	inline void set_range(float min, float max)
	{
		lb_value = min;
		ub_value = max;
		if (colour_num)
			inv_len = (max - min) / float(colour_num-1);
	}

	inline float get_lower_bound() { return lb_value; }
	inline float get_upper_bound() { return ub_value; }
	inline float get_range() { return (ub_value - lb_value); }

	inline Colour get_colour(float value)
	{
		if (value < lb_value)
			return lb_colour;
		else if (value >= ub_value)
			return ub_colour;
		// interpolate
		float w1 = (value - lb_value) / inv_len;
		size_t index = size_t(w1);
		w1 -= float(index);
		float w2 = 1.0f - w1;
		Colour& c1 = colours[index];
		Colour& c2 = colours[index + 1];
		Colour res;
		res.r = w1 * c1.r + w2 * c2.r;
		res.g = w1 * c1.g + w2 * c2.g;
		res.b = w1 * c1.b + w2 * c2.b;
		return res;
	}

	unsigned char *gen_1Dtexture(size_t inv_resolution, size_t &texture_size);
};

#endif