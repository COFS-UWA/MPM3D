#ifndef __Color_Graph_h__
#define __Color_Graph_h__

// map value to color
class ColorGraph
{
public:
	struct Colori { int   r, g, b; };
	struct Colorf { float r, g, b; };
	struct ValueColorPair { double va; float r, g, b; };

	ColorGraph() : pairs(nullptr), color_num(0)
	{
		// white
		lower_bound_color.r = 1.0f;
		lower_bound_color.g = 1.0f;
		lower_bound_color.b = 1.0f;
		// grey
		upper_bound_color.r = 0.5f;
		upper_bound_color.g = 0.5f;
		upper_bound_color.b = 0.5f;
	}
	~ColorGraph() { clear(); }
	void clear()
	{
		if (pairs)
		{
			delete[] pairs;
			pairs = nullptr;
		}
		color_num = 0;
	}

	int init(ValueColorPair *vcps, size_t num)
	{
		clear();
		if (!vcps || num == 0)
			return -1;

		pairs = new ValueColorPair[num];
		color_num = num;
		for (size_t i = 0; i < num; ++i)
		{
			pairs[i].va = vcps[i].va;
			pairs[i].r = vcps[i].r;
			pairs[i].g = vcps[i].g;
			pairs[i].b = vcps[i].b;
		}
		return 0;
	}

	int init(double lower, double upper, Colori *colors, size_t num,
			 bool apply_out_of_bound_color = false)
	{
		clear();
		if (!colors || num == 0) return -1;
		pairs = new ValueColorPair[num];
		color_num = num;
		double intv = (upper - lower) / double(num - 1);
		for (size_t i = 0; i < num; ++i)
		{
			pairs[i].va = lower + intv * double(i);
			pairs[i].r = float(colors[i].r) / 255.0f;
			pairs[i].g = float(colors[i].g) / 255.0f;
			pairs[i].b = float(colors[i].b) / 255.0f;
		}
		if (apply_out_of_bound_color)
		{
			lower_bound_color.r = colors[0].r;
			lower_bound_color.g = colors[0].g;
			lower_bound_color.b = colors[0].b;
			upper_bound_color.r = colors[num - 1].r;
			upper_bound_color.g = colors[num - 1].g;
			upper_bound_color.b = colors[num - 1].b;
		}
		return 0;
	}

	inline Colorf get_color(double va)
	{
		if (va < pairs[0].va)
		{
			return lower_bound_color;
		}
		else if (va > pairs[color_num - 1].va)
		{
			return upper_bound_color;
		}

		// find va is in which internval
		Colorf res;
		size_t low_id = 0, up_id = color_num - 1, mid_id;
		mid_id = (low_id + up_id) / 2;
		do
		{
			if (pairs[mid_id].va > va)
				up_id  = mid_id;
			else
				low_id = mid_id;
			mid_id = (low_id + up_id) / 2;
		} while (low_id != mid_id);
		
		// interpolate color (can use nonlinear interpolation in the future)
		res.r = pairs[low_id].r + (pairs[up_id].r - pairs[low_id].r) / (pairs[up_id].va - pairs[low_id].va) * (va - pairs[low_id].va);
		res.g = pairs[low_id].g + (pairs[up_id].g - pairs[low_id].g) / (pairs[up_id].va - pairs[low_id].va) * (va - pairs[low_id].va);
		res.b = pairs[low_id].b + (pairs[up_id].b - pairs[low_id].b) / (pairs[up_id].va - pairs[low_id].va) * (va - pairs[low_id].va);
		
		return res;
	}

protected:
	size_t color_num;
	ValueColorPair *pairs;
	Colorf lower_bound_color, upper_bound_color;

protected:
};


#endif