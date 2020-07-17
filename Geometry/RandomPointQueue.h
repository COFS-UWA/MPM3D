#ifndef __Random_Point_Queue_h__
#define __Random_Point_Queue_h__

#include <unordered_map>

#include "Geometry.h"

class RandomPointQueue
{
public:
	struct Point
	{
		double x, y;
		double l_pt_dist;
	};

protected:
	typedef std::unordered_map<size_t, Point> PointBuffer;
	PointBuffer point_buf;
	size_t point_num;

public:
	RandomPointQueue();
	~RandomPointQueue();

	template <typename Point>
	void add_point(Point &p);

	void add_point(Point2D &p);
	bool get_point(Point2D &p);
	bool is_empty();
};

template <typename Point>
inline void RandomPointQueue::add_point(Point& p)
{
	Point2D p2d;
	p2d.x = p.x;
	p2d.y = p.y;
	add_point(p2d);
}

#endif