#ifndef __Grid_2D_hpp__
#define __Grid_2D_hpp__

template <typename Grid = void>
struct Grid2D
{
	double xl, yl;
	double hx, hy;
	size_t x_num, y_num;
	double xu, yu;
	Grid *grids;

	size_t num;
	// max_x_id = grid_x_num - 1
	size_t max_x_id, max_y_id;
	// xu_min_hx = grid_xu - grid_hx
	double xu_min_hx, yu_min_hy;
	
	inline Grid2D() : grids(nullptr), x_num(0), y_num(0), num(0) {}
	~Grid2D() { clear(); }

	inline size_t offset_from_xy_id(size_t x_id, size_t y_id) const noexcept
	{ return y_id * x_num + x_id;}

	inline Grid& grid_by_xy_id(size_t x_id, size_t y_id) noexcept
	{ return grids[y_id * x_num + x_id]; }
	
	inline size_t get_x_id(double x) const noexcept
	{ return x < xl ? 0 : (x > xu_min_hx ? max_x_id : size_t((x - xl) / hx)); }
	inline size_t get_y_id(double y) const noexcept
	{ return y < yl ? 0 : (y > yu_min_hy ? max_y_id : size_t((y - yl) / hy)); }
	
	inline void clear() noexcept
	{
		if (grids)
		{
			delete[] grids;
			grids = nullptr;
		}
		x_num = 0;
		y_num = 0;
		num = 0;
	}

	int alloc_grid(
		double _xl, double _yl,
		double _hx, double _hy,
		size_t _x_num, size_t _y_num)
	{
		clear();
		if (_hx == 0.0 || _hy == 0.0 || _x_num == 0 || _y_num == 0 )
			return -1;
		xl = _xl; yl = _yl;
		xu = _xl + double(_x_num) * _hx;
		yu = _yl + double(_y_num) * _hy;
		hx = _hx; hy = _hy;
		x_num = _x_num;
		y_num = _y_num;
		init_cal_var();
		grids = new Grid[num];
		return 0;
	}

	int alloc_grid(double _xl, double _yl,
				   double _xu, double _yu,
				   double _hx, double _hy,
				   double exp_ratio = 0.0)
	{
		if (_xl >= _xu || _yl >= _yu || _hx == 0.0 || _hy == 0.0)
			return -1;
		clear();
		double x_len = _xu - _xl;
		double y_len = _yu - _yl;
		double x_padding = x_len * exp_ratio;
		double y_padding = y_len * exp_ratio;
		x_len += x_padding;
		y_len += y_padding;
		x_padding *= 0.5;
		y_padding *= 0.5;
		xl = _xl - x_padding; xu = _xu + x_padding;
		yl = _yl - y_padding; yu = _yu + y_padding;
		x_num = size_t(ceil(x_len / _hx));
		y_num = size_t(ceil(y_len / _hy));
		hx = x_len / double(x_num);
		hy = y_len / double(y_num);
		init_cal_var();
		grids = new Grid[num];
		return 0;
	}

	template <typename Grid2>
	int alloc_grid(const Grid2 &other)
	{
		clear();
		if (other.hx == 0.0 || other.hy == 0.0 || other.x_num == 0 || other.y_num == 0)
			return -1;
		xl = other.xl; yl = other.yl;
		xu = other.xl + double(other.x_num) * other.hx;
		yu = other.yl + double(other.y_num) * other.hy;
		hx = other.hx; hy = other.hy;
		x_num = other.x_num;
		y_num = other.y_num;
		init_cal_var();
		grids = new Grid[num];
		return 0;
	}

protected:
	inline void init_cal_var() noexcept
	{
		num = x_num * y_num;
		max_x_id = x_num - 1;
		max_y_id = y_num - 1;
		xu_min_hx = xu - hx;
		yu_min_hy = yu - hy;
	}
};

template <>
struct Grid2D<void>
{
	double xl, yl;
	double xu, yu;
	double hx, hy;
	size_t x_num, y_num, num;
	// max_x_id = grid_x_num - 1
	size_t max_x_id, max_y_id;
	// xu_min_hx = grid_xu - grid_hx
	double xu_min_hx, yu_min_hy;

	inline Grid2D() : x_num(0), y_num(0), num(0) {}
	~Grid2D() { clear(); }

	inline size_t offset_from_xy_id(size_t x_id, size_t y_id) const noexcept
	{ return y_id * x_num + x_id; }

	inline size_t get_x_id(double x) const noexcept
	{ return x < xl ? 0 : (x > xu_min_hx ? max_x_id : size_t((x - xl) / hx)); }
	inline size_t get_y_id(double y) const noexcept
	{ return y < yl ? 0 : (y > yu_min_hy ? max_y_id : size_t((y - yl) / hy)); }

	inline void clear() noexcept
	{ x_num = 0; y_num = 0; num = 0; }

	int alloc_grid(double _xl, double _yl,
				   double _hx, double _hy,
				   size_t _x_num, size_t _y_num)
	{
		if (_hx == 0.0 || _hy == 0.0 || _x_num == 0 || _y_num == 0)
			return -1;
		xl = _xl; yl = _yl;
		xu = _xl + double(_x_num) * _hx;
		yu = _yl + double(_y_num) * _hy;
		hx = _hx; hy = _hy;
		x_num = _x_num;
		y_num = _y_num;
		init_cal_var();
		return 0;
	}

	int alloc_grid(
		double _xl, double _yl,
		double _xu, double _yu,
		double _hx, double _hy,
		double exp_ratio = 0.0)
	{
		if (_xl >= _xu || _yl >= _yu || _hx == 0.0 || _hy == 0.0)
			return -1;
		double x_len = _xu - _xl;
		double y_len = _yu - _yl;
		double x_padding = x_len * exp_ratio;
		double y_padding = y_len * exp_ratio;
		x_len += x_padding;
		y_len += y_padding;
		x_padding *= 0.5;
		y_padding *= 0.5;
		xl = _xl - x_padding; xu = _xu + x_padding;
		yl = _yl - y_padding; yu = _yu + y_padding;
		x_num = size_t(ceil(x_len / _hx));
		y_num = size_t(ceil(y_len / _hy));
		hx = x_len / double(x_num);
		hy = y_len / double(y_num);
		init_cal_var();
		return 0;
	}

	template <typename Grid2>
	int alloc_grid(const Grid2& other)
	{
		if (other.hx == 0.0 || other.hy == 0.0 || other.x_num == 0 || other.y_num == 0)
			return -1;
		xl = other.xl; yl = other.yl;
		xu = other.xl + double(other.x_num) * other.hx;
		yu = other.yl + double(other.y_num) * other.hy;
		hx = other.hx; hy = other.hy;
		x_num = other.x_num;
		y_num = other.y_num;
		init_cal_var();
		return 0;
	}

protected:
	inline void init_cal_var() noexcept
	{
		num = x_num * y_num;
		max_x_id = x_num - 1;
		max_y_id = y_num - 1;
		xu_min_hx = xu - hx;
		yu_min_hy = yu - hy;
	}
};

#endif