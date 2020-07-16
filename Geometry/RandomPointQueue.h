#ifndef __Random_Point_Queue_h__
#define __Random_Point_Queue_h__

#include <unordered_map>

#include "Geometry.h"

class RandomPointQueue
{
protected:
	typedef std::unordered_map<size_t, Point2D> PointBuffer;
	//typedef std::pair<size_t, Point2D> PointItem;
	
	PointBuffer point_buf;
	size_t point_num;

public:
	RandomPointQueue();
	~RandomPointQueue();

	void add_point(Point2D &p);
	bool get_point(Point2D &p);
	bool is_empty();
};

#endif