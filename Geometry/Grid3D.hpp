#ifndef __Grid_3D_hpp__
#define __Grid_3D_hpp__

template <typename Grid = void>
struct Grid3D
{
	double xl, yl, zl;
	double hx, hy, hz;
	size_t x_num, y_num, z_num;
	double xu, yu, zu;
	Grid *grids;

	size_t xy_num, num;
	// max_x_id = grid_x_num - 1
	size_t max_x_id, max_y_id, max_z_id;
	// xu_min_hx = grid_xu - grid_hx
	double xu_min_hx, yu_min_hy, zu_min_hz;
	
	inline Grid3D() : grids(nullptr),
		x_num(0), y_num(0), z_num(0), num(0) {}
	~Grid3D() { clear(); }

	inline size_t offset_from_xyz_id(size_t x_id,
		size_t y_id, size_t z_id) const noexcept
	{ return z_id * xy_num + y_id * x_num + x_id;}

	inline Grid& grid_by_xyz_id(size_t x_id,
		size_t y_id, size_t z_id) noexcept
	{ return grids[z_id * xy_num + y_id * x_num + x_id]; }
	
	inline size_t get_x_id(double x) const noexcept
	{ return x < xl ? 0 : (x > xu_min_hx ? max_x_id : size_t((x - xl) / hx)); }
	inline size_t get_y_id(double y) const noexcept
	{ return y < yl ? 0 : (y > yu_min_hy ? max_y_id : size_t((y - yl) / hy)); }
	inline size_t get_z_id(double z) const noexcept
	{ return z < zl ? 0 : (z > zu_min_hz ? max_z_id : size_t((z - zl) / hz)); }

	inline void clear() noexcept
	{
		if (grids)
		{
			delete[] grids;
			grids = nullptr;
		}
		x_num = 0;
		y_num = 0;
		z_num = 0;
		num = 0;
	}

	int alloc_grid(double _xl, double _yl, double _zl,
				   double _hx, double _hy, double _hz,
				   size_t _x_num, size_t _y_num, size_t _z_num)
	{
		clear();
		if (_hx == 0.0 || _hy == 0.0 || _hz == 0.0 ||
			_x_num == 0 || _y_num == 0 || _z_num == 0)
			return -1;
		xl = _xl; yl = _yl; zl = _zl;
		xu = _xl + double(_x_num) * _hx;
		yu = _yl + double(_y_num) * _hy;
		zu = _zl + double(_z_num) * _hz;
		hx = _hx; hy = _hy; hz = _hz;
		x_num = _x_num;
		y_num = _y_num;
		z_num = _z_num;
		init_cal_var();
		grids = new Grid[num];
		return 0;
	}

	int alloc_grid(double _xl, double _yl, double _zl,
				   double _xu, double _yu, double _zu,
				   double _hx, double _hy, double _hz,
				   double exp_ratio = 0.0)
	{
		if (_xl >= _xu || _yl >= _yu || _zl >= _zu ||
			_hx == 0.0 || _hy == 0.0 || _hz == 0.0)
			return -1;
		clear();
		double x_len = _xu - _xl;
		double y_len = _yu - _yl;
		double z_len = _zu - _zl;
		double x_padding = x_len * exp_ratio;
		double y_padding = y_len * exp_ratio;
		double z_padding = z_len * exp_ratio;
		x_len += x_padding;
		y_len += y_padding;
		z_len += z_padding;
		x_padding *= 0.5;
		y_padding *= 0.5;
		z_padding *= 0.5;
		xl = _xl - x_padding; xu = _xu + x_padding;
		yl = _yl - y_padding; yu = _yu + y_padding;
		zl = _zl - z_padding; zu = _zu + z_padding;
		x_num = size_t(ceil(x_len / _hx));
		y_num = size_t(ceil(y_len / _hy));
		z_num = size_t(ceil(z_len / _hz));
		hx = x_len / double(x_num);
		hy = y_len / double(y_num);
		hz = z_len / double(z_num);
		init_cal_var();
		grids = new Grid[num];
		return 0;
	}

	template <typename Grid2>
	int alloc_grid(const Grid2 &other)
	{
		clear();
		if (other.hx == 0.0 || other.hy == 0.0 || other.hz == 0.0 ||
			other.x_num == 0 || other.y_num == 0 || other.z_num == 0)
			return -1;
		xl = other.xl; yl = other.yl; zl = other.zl;
		xu = other.xl + double(other.x_num) * other.hx;
		yu = other.yl + double(other.y_num) * other.hy;
		zu = other.zl + double(other.z_num) * other.hz;
		hx = other.hx; hy = other.hy; hz = other.hz;
		x_num = other.x_num;
		y_num = other.y_num;
		z_num = other.z_num;
		init_cal_var();
		grids = new Grid[num];
		return 0;
	}

protected:
	inline void init_cal_var() noexcept
	{
		xy_num = x_num * y_num;
		num = xy_num * z_num;
		max_x_id = x_num - 1;
		max_y_id = y_num - 1;
		max_z_id = z_num - 1;
		xu_min_hx = xu - hx;
		yu_min_hy = yu - hy;
		zu_min_hz = zu - hz;
	}
};

template <>
struct Grid3D<void>
{
	double xl, yl, zl;
	double hx, hy, hz;
	size_t x_num, y_num, z_num;
	double xu, yu, zu;

	size_t xy_num, num;
	// max_x_id = grid_x_num - 1
	size_t max_x_id, max_y_id, max_z_id;
	// xu_min_hx = grid_xu - grid_hx
	double xu_min_hx, yu_min_hy, zu_min_hz;

	inline Grid3D() : x_num(0), y_num(0), z_num(0), num(0) {}
	~Grid3D() { clear(); }

	inline size_t offset_from_xyz_id(size_t x_id,
		size_t y_id, size_t z_id) const noexcept
	{ return z_id * xy_num + y_id * x_num + x_id; }

	inline size_t get_x_id(double x) const noexcept
	{ return x < xl ? 0 : (x > xu_min_hx ? max_x_id : size_t((x - xl) / hx)); }
	inline size_t get_y_id(double y) const noexcept
	{ return y < yl ? 0 : (y > yu_min_hy ? max_y_id : size_t((y - yl) / hy)); }
	inline size_t get_z_id(double z) const noexcept
	{ return z < zl ? 0 : (z > zu_min_hz ? max_z_id : size_t((z - zl) / hz)); }

	inline void clear() noexcept
	{ x_num = 0; y_num = 0; z_num = 0; num = 0; }

	int alloc_grid(double _xl, double _yl, double _zl,
		double _hx, double _hy, double _hz,
		size_t _x_num, size_t _y_num, size_t _z_num)
	{
		if (_hx == 0.0 || _hy == 0.0 || _hz == 0.0 ||
			_x_num == 0 || _y_num == 0 || _z_num == 0)
			return -1;
		xl = _xl; yl = _yl; zl = _zl;
		xu = _xl + double(_x_num) * _hx;
		yu = _yl + double(_y_num) * _hy;
		zu = _zl + double(_z_num) * _hz;
		hx = _hx; hy = _hy; hz = _hz;
		x_num = _x_num;
		y_num = _y_num;
		z_num = _z_num;
		init_cal_var();
		return 0;
	}

	int alloc_grid(double _xl, double _yl, double _zl,
		double _xu, double _yu, double _zu,
		double _hx, double _hy, double _hz,
		double exp_ratio = 0.0)
	{
		if (_xl >= _xu || _yl >= _yu || _zl >= _zu ||
			_hx == 0.0 || _hy == 0.0 || _hz == 0.0)
			return -1;
		double x_len = _xu - _xl;
		double y_len = _yu - _yl;
		double z_len = _zu - _zl;
		double x_padding = x_len * exp_ratio;
		double y_padding = y_len * exp_ratio;
		double z_padding = z_len * exp_ratio;
		x_len += x_padding;
		y_len += y_padding;
		z_len += z_padding;
		x_padding *= 0.5;
		y_padding *= 0.5;
		z_padding *= 0.5;
		xl = _xl - x_padding; xu = _xu + x_padding;
		yl = _yl - y_padding; yu = _yu + y_padding;
		zl = _zl - z_padding; zu = _zu + z_padding;
		x_num = size_t(ceil(x_len / _hx));
		y_num = size_t(ceil(y_len / _hy));
		z_num = size_t(ceil(z_len / _hz));
		hx = x_len / double(x_num);
		hy = y_len / double(y_num);
		hz = z_len / double(z_num);
		init_cal_var();
		return 0;
	}

	template <typename Grid2>
	int alloc_grid(const Grid2& other)
	{
		if (other.hx == 0.0 || other.hy == 0.0 || other.hz == 0.0 ||
			other.x_num == 0 || other.y_num == 0 || other._z_num == 0)
			return -1;
		xl = other.xl; yl = other.yl; zl = other.zl;
		xu = other.xl + double(other.x_num) * other.hx;
		yu = other.yl + double(other.y_num) * other.hy;
		zu = other.zl + double(other.z_num) * other.hz;
		hx = other.hx; hy = other.hy; hz = other.hz;
		x_num = other.x_num;
		y_num = other.y_num;
		z_num = other.z_num;
		init_cal_var();
		return 0;
	}

protected:
	inline void init_cal_var() noexcept
	{
		xy_num = x_num * y_num;
		num = xy_num * z_num;
		max_x_id = x_num - 1;
		max_y_id = y_num - 1;
		max_z_id = z_num - 1;
		xu_min_hx = xu - hx;
		yu_min_hy = yu - hy;
		zu_min_hz = zu - hz;
	}
};

#endif