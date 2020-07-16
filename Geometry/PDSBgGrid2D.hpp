#ifndef __PDS_Bg_Grid_2D_hpp__
#define __PDS_Bg_Grid_2D_hpp__

#include <cmath>

// next_offset means the offset of next pointer in Point2D type
template <typename Point2D, size_t next_offset>
class PDSBgGrid2D
{
protected:
	struct Cell
	{
		size_t x_id;
		size_t y_id;
		Point2D *top;
	};

	double xl, yl;
	double xu, yu;
	double hx, hy;
	size_t x_num, y_num;
	Cell* cells;

public:
	PDSBgGrid2D() : x_num(0), y_num(0), cells(nullptr) {}
	~PDSBgGrid2D() { clear(); }

	void clear()
	{
		if (cells)
			delete[] cells;
		cells = nullptr;
		x_num = 0;
		y_num = 0;
	}

	inline Cell& get_cell(size_t x_id, size_t y_id)
	{ return cells[y_id * x_num + x_id]; }

	int init(double _xl, double _xu,
			 double _yl, double _yu,
			 size_t _x_num, size_t _y_num)
	{
		clear();
		xl = _xl;
		xu = _xu;
		x_num = _x_num;
		hx = (xu - xl) / double(x_num);
		yl = _yl;
		yu = _yu;
		y_num = _y_num;
		hy = (yu - yl) / double(y_num);
		cells = new Cell[x_num * y_num];
		size_t cell_id = 0;
		for (size_t y_id = 0; y_id < y_num; ++y_id)
			for (size_t x_id = 0; x_id < x_num; ++x_id)
			{
				Cell& c = cells[cell_id];
				c.x_id = x_id;
				c.y_id = y_id;
				c.top = nullptr;
				++cell_id;
			}
		return 0;
	}

	inline bool is_in_grid(Point2D& p)
	{
		return p.x >= xl && p.x < xu &&
			   p.y >= yl && p.y < yu;
	}

	inline bool add_point(Point2D& p)
	{
		if (p.x < xl || p.x >= xu ||
			p.y < yl || p.y >= yu)
			return false;
		size_t x_id = size_t(floor((p.x - xl) / hx));
		size_t y_id = size_t(floor((p.y - yl) / hy));
		Cell& c = get_cell(x_id, y_id);
		// p.next = c.top;
		*(Point2D **)((char*)(&p) + next_offset) = c.top;
		c.top = &p;
		return true;
	}

	bool has_point_nearby(Point2D& p, double dist)
	{
		double pxl = p.x - dist;
		if (pxl < xl)
			pxl = xl;
		size_t xl_id = size_t(floor(pxl - xl) / hx);
		if (xl_id > x_num - 1)
			xl_id = x_num - 1;

		double pxu = p.x + dist;
		if (pxu > xu)
			pxu = xu;
		size_t xu_id = size_t(ceil(pxu - xl) / hx);
		if (xu_id > x_num - 1)
			xu_id = x_num - 1;

		double pyl = p.y - dist;
		if (pyl < yl)
			pyl = yl;
		size_t yl_id = size_t(floor(pyl - yl) / hy);
		if (yl_id > y_num - 1)
			yl_id = y_num - 1;

		double pyu = p.y + dist;
		if (pyu > yu)
			pyu = yu;
		size_t yu_id = size_t(ceil(pyu - yl) / hy);
		if (yu_id > y_num - 1)
			yu_id = y_num - 1;

		double dist2 = dist * dist;
		double dx, dy, dd2;
		for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
			for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
			{
				Cell& c = get_cell(x_id, y_id);
				for (Point2D* pt = c.top; pt;
					 pt = *(Point2D **)((char *)pt + next_offset))
				{
					dx = pt->x - p.x;
					dy = pt->y - p.y;
					dd2 = dx * dx + dy * dy;
					if (dd2 < dist2)
						return true;
				}
			}
		return false;
	}
};

#endif