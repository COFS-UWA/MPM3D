#include <cstdlib>

#include "BgGrid.h"

int BgGrid::init(double _xl, double _xu,
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

void BgGrid::clear()
{
	if (cells)
		delete[] cells;
	cells = nullptr;
	x_num = 0;
	y_num = 0;
}

bool BgGrid::add_point(Point2D& p)
{
	if (p.x < xl || p.x > xu ||
		p.y < yl || p.y > yu)
		return false;

	size_t x_id = size_t(floor((p.x - xl) / hx));
	size_t y_id = size_t(floor((p.y - yl) / hy));
	Cell& c = get_cell(x_id, y_id);
	p.next = c.top;
	c.top = &p;
	return true;
}

bool BgGrid::is_in_grid(Point2D& p)
{
	return p.x >= xl && p.x <= xu && p.y >= yl && p.y <= yu;
}

bool BgGrid::has_point_nearby(Point2D& p, double dist)
{
	double pxl = p.x - dist;
	if (pxl < xl)
		pxl = xl;
	size_t xl_id = size_t(floor(pxl - xl) / hx);

	double pxu = p.x + dist;
	if (pxu > xu)
		pxu = xu;
	size_t xu_id = size_t(ceil(pxu - xl) / hx);

	double pyl = p.y - dist;
	if (pyl < yl)
		pyl = yl;
	size_t yl_id = size_t(floor(pyl - yl) / hy);

	double pyu = p.y + dist;
	if (pyu > yu)
		pyu = yu;
	size_t yu_id = size_t(ceil(pyu - yl) / hy);

	double dist2 = dist * dist;
	double dx, dy, dd2;
	for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
		for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
		{
			Cell& c = get_cell(x_id, y_id);
			for (Point2D *pt = c.top; pt; pt = pt->next)
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