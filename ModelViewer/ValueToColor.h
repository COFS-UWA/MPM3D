#ifndef __Value_To_Color_h__
#define __Value_To_Color_h__

#include <cmath>

// Map value to color
class ValueToColor
{
public:
	struct Colori { unsigned char r, g, b; };

	struct Colorf
	{
		float r, g, b;
		Colorf() {}
		Colorf(float _r, float _g, float _b)
		{
			r = _r;
			g = _g;
			b = _b;
		}
		Colorf(Colori &c)
		{
			r = float(c.r) / 255.0f;
			g = float(c.g) / 255.0f;
			b = float(c.b) / 255.0f;
		}
	};

	struct ValueColorPair
	{
		double value;
		union
		{
			Colorf color;
			struct { float r, g, b; };
		};
		ValueColorPair() {}
	};

protected:
	size_t pair_num;
	ValueColorPair *pairs;
	Colorf lb_color; // lower bound color
	Colorf ub_color; // upper bound color

public:
	ValueToColor() : pair_num(0), pairs(nullptr)
	{
		// white
		ub_color.r = 1.0f;
		ub_color.g = 1.0f;
		ub_color.b = 1.0f;
		// grey
		lb_color.r = 0.3f;
		lb_color.g = 0.3f;
		lb_color.b = 0.3f;
	}
	~ValueToColor() { clear(); }
	void clear()
	{
		if (pairs)
		{
			delete[] pairs;
			pairs = nullptr;
		}
		pair_num = 0;
	}

	// support un-uniform color scale 
	int init(ValueColorPair *vcps, size_t vcp_num)
	{
		clear();
		if (!vcps || vcp_num <2)
			return -1;

		pairs = new ValueColorPair[vcp_num];
		pair_num = vcp_num;
		for (size_t p_id = 0; p_id < pair_num; ++p_id)
		{
			pairs[p_id].value = vcps[p_id].value;
			pairs[p_id].r = vcps[p_id].r;
			pairs[p_id].g = vcps[p_id].g;
			pairs[p_id].b = vcps[p_id].b;
		}
		return 0;
	}

	inline Colorf get_color(double value)
	{
		if (value < pairs[0].value)
			return lb_color;
		else if (value > pairs[pair_num-1].value)
			return ub_color;

		// find value is in which internval
		Colorf res;
		size_t l_id = 0;
		size_t u_id = pair_num - 1;
		size_t m_id = (l_id + u_id) / 2;
		do
		{
			if (pairs[m_id].value > value)
				u_id = m_id;
			else
				l_id = m_id;
			m_id = (l_id + u_id) / 2;
		} while (l_id != m_id);
		
		// interpolate color (can use nonlinear interpolation in the future)
		ValueColorPair &pair1 = pairs[l_id];
		ValueColorPair &pair2 = pairs[u_id];
		res.r = pair1.r + (pair2.r - pair1.r) / (pair2.value - pair1.value) * (value - pair1.value);
		res.g = pair1.g + (pair2.g - pair1.g) / (pair2.value - pair1.value) * (value - pair1.value);
		res.b = pair1.b + (pair2.b - pair1.b) / (pair2.value - pair1.value) * (value - pair1.value);
		return res;
	}

	// unform color scale
	int init(double lower, double upper,
		Colori *colors, size_t color_num,
		bool apply_out_of_bound_color = false)
	{
		clear();
		if (color_num < 2 || !colors)
			return -1;

		pair_num = color_num;
		pairs = new ValueColorPair[pair_num];
		double inv_len = (upper - lower) / double(pair_num - 1);
		for (size_t p_id = 0; p_id < pair_num; ++p_id)
		{
			pairs[p_id].value = lower + inv_len * double(p_id);
			pairs[p_id].color = colors[p_id];
		}

		if (apply_out_of_bound_color)
		{
			lb_color = pairs[0].color;
			ub_color = pairs[pair_num-1].color;
		}

		return 0;
	}
	
	//inline Colorf get_color(double value)
	// assume uniform color scale
};


#endif