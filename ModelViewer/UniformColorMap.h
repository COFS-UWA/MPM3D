#ifndef __Uniform_Color_Map_h__
#define __Uniform_Color_Map_h__

#include <cmath>

class UniformColorMap
{
public:
	struct Color
	{
		float r, g, b;
		Color() {}
		Color(float _r, float _g, float _b) :
			r(_r), g(_g), b(_b) {}
	};

protected:
	size_t color_num;
	Color *colors;
	Color lb_color, ub_color;

	float lb_value, ub_value;
	float inv_len;

	void clear();

public:
	UniformColorMap();
	~UniformColorMap();

	int set_color(
		const unsigned char color_datas[][3],
		size_t color_data_num,
		bool apply_ob_color = true
		);

	inline void set_range(float min, float max)
	{
		lb_value = min;
		ub_value = max;
		if (color_num)
			inv_len = (max - min) / float(color_num-1);
	}

	inline float get_lower_bound() { return lb_value; }
	inline float get_upper_bound() { return ub_value; }
	inline float get_range() { return (ub_value - lb_value); }
	
	inline size_t get_color_num() { return color_num; }
	inline Color* get_colors() { return colors; }

	inline Color get_color(float value)
	{
		if (value < lb_value)
			return lb_color;
		else if (value >= ub_value)
			return ub_color;
		// interpolate
		float w1 = (value - lb_value) / inv_len;
		size_t index = size_t(w1);
		w1 -= float(index);
		float w2 = 1.0f - w1;
		Color& c1 = colors[index];
		Color& c2 = colors[index + 1];
		Color res;
		res.r = w1 * c1.r + w2 * c2.r;
		res.g = w1 * c1.g + w2 * c2.g;
		res.b = w1 * c1.b + w2 * c2.b;
		return res;
	}

	unsigned char *gen_1Dtexture(size_t inv_resolution, size_t &texture_size);
};

#endif