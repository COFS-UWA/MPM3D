#ifndef __Grid_3D_hpp__
#define __Grid_3D_hpp__

template <typename Grid = void>
struct Grid3D
{
	double xl, yl, zl;
	double hx, hy, hz;
	size_t x_num, y_num, z_num;
	double xu, yu, zu;
	size_t xy_num, num;
	Grid *grids;

	inline Grid3D() : grids(nullptr),
		x_num(0), y_num(0), z_num(0), num(0) {}
	~Grid3D() { clear(); }

	inline size_t offset_from_xyz_id(size_t x_id,
		size_t y_id, size_t z_id) const noexcept
	{ return z_id * xy_num + y_id * x_num + x_id;}

	inline Grid& grid_by_xyz_id(size_t x_id,
		size_t y_id, size_t z_id) noexcept
	{ return grids[z_id * xy_num + y_id * x_num + x_id]; }
	
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
		xy_num = x_num * y_num;
		num = xy_num * z_num;
		grids = new Grid[num];
		return 0;
	}

	int alloc_grid(double _xl, double _yl, double _zl,
				   double _xu, double _yu, double _zu,
				   double _hx, double _hy, double _hz)
	{
		if (_xl >= _xu || _yl >= _yu || _zl >= _zu ||
			_hx == 0.0 || _hy == 0.0 || _hz == 0.0)
			return -1;
		clear();
		xl = _xl; yl = _yl; zl = _zl;
		xu = _xu; yu = _yu; zu = _zu;
		x_num = size_t(ceil((_xu - _xl) / _hx));
		y_num = size_t(ceil((_yu - _yl) / _hy));
		z_num = size_t(ceil((_zu - _zl) / _hz));
		hx = (_xu - _xl) / double(x_num);
		hy = (_yu - _yl) / double(y_num);
		hz = (_zu - _zl) / double(z_num);
		xy_num = x_num * y_num;
		num = xy_num * z_num;
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
		xy_num = x_num * y_num;
		num = xy_num * z_num;
		grids = new Grid[num];
		return 0;
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

	inline Grid3D() : x_num(0), y_num(0), z_num(0), num(0) {}
	~Grid3D() { clear(); }

	inline size_t offset_from_xyz_id(size_t x_id,
		size_t y_id, size_t z_id) const noexcept
	{ return z_id * xy_num + y_id * x_num + x_id; }

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
		xy_num = x_num * y_num;
		num = xy_num * z_num;
		return 0;
	}

	int alloc_grid(double _xl, double _yl, double _zl,
		double _xu, double _yu, double _zu,
		double _hx, double _hy, double _hz)
	{
		if (_xl >= _xu || _yl >= _yu || _zl >= _zu ||
			_hx == 0.0 || _hy == 0.0 || _hz == 0.0)
			return -1;
		xl = _xl; yl = _yl; zl = _zl;
		xu = _xu; yu = _yu; zu = _zu;
		x_num = size_t(ceil((_xu - _xl) / _hx));
		y_num = size_t(ceil((_yu - _yl) / _hy));
		z_num = size_t(ceil((_zu - _zl) / _hz));
		hx = (_xu - _xl) / double(x_num);
		hy = (_yu - _yl) / double(y_num);
		hz = (_zu - _zl) / double(z_num);
		xy_num = x_num * y_num;
		num = xy_num * z_num;
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
		xy_num = x_num * y_num;
		num = xy_num * z_num;
		return 0;
	}
};

#endif