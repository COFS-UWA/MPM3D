#ifndef __Geometry_Utils_h__
#define __Geometry_Utils_h__

// limit theta to [-pi, pi]
#define PI 3.14159265359
inline void trim_to_pi(double& theta)
{
	if (theta > PI)
		theta -= (2.0 * PI) * long long((theta + PI) / (2.0 * PI));
	else if (theta < -PI)
		theta -= (2.0 * PI) * long long((theta - PI) / (2.0 * PI));
}
#undef PI

template <typename Item>
inline void sort_array_3_acc(Item ids[3])
{
	Item tmp;
	size_t min_id = 0;
	if (ids[1] < ids[0])
		min_id = 1;
	if (ids[2] < ids[min_id])
		min_id = 2;
	if (min_id != 0)
	{
		tmp = ids[0];
		ids[0] = ids[min_id];
		ids[min_id] = tmp;
	}
	if (ids[2] < ids[1])
	{
		tmp = ids[1];
		ids[1] = ids[2];
		ids[2] = tmp;
	}
}

template <typename Item>
inline void sort_array_4_acc(Item ids[4])
{
	Item tmp;
	size_t min_id;
	for (size_t i = 0; i < 3; ++i)
	{
		min_id = i;
		for (size_t j = i + 1; j < 4; ++j)
		{
			if (ids[j] < ids[min_id])
				min_id = j;
		}
		if (min_id != i)
		{
			tmp = ids[min_id];
			ids[min_id] = ids[i];
			ids[i] = tmp;
		}
	}
}

#endif