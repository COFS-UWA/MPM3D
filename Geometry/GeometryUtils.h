#ifndef __Geometry_Utils_h__
#define __Geometry_Utils_h__

#include <unordered_map>

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

namespace GeometryUtils_Internal
{
	class LineMap
	{
	public:
		inline LineMap() {}
		inline ~LineMap() {}

		inline void add_line(size_t n1, size_t n2)
		{
			if (n1 > n2) // swap
			{
				size_t tmp = n1;
				n1 = n2;
				n2 = tmp;
			}
			char nid_str[50];
			snprintf(nid_str, sizeof(nid_str), "%zu_%zu", n1, n2);
			auto res = map.emplace(std::string(nid_str), LineMapItem(n1, n2));
		}

		inline size_t get_line_num() const noexcept { return map.size(); }

		template <typename Line>
		void output_lines(Line* lines)
		{
			Line* cur_line = lines;
			for (Map::iterator iter = map.begin(); iter != map.end(); ++iter)
			{
				LineMapItem& li = iter->second;
				cur_line->n1 = li.n1;
				cur_line->n2 = li.n2;
				++cur_line;
			}
		}

	protected:
		struct LineMapItem
		{
			size_t n1, n2;
			inline LineMapItem(size_t _n1, size_t _n2) : n1(_n1), n2(_n2) {}
		};
		typedef std::unordered_map<std::string, LineMapItem> Map;
		Map map;
	};
}

template <typename Line, typename Triangle>
Line* extract_edges_from_triangles(
	const Triangle *elems,
	size_t elem_num,
	size_t &line_num
	)
{
	GeometryUtils_Internal::LineMap line_map;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		const Triangle& e = elems[e_id];
		line_map.add_line(e.n1, e.n2);
		line_map.add_line(e.n1, e.n3);
		line_map.add_line(e.n2, e.n3);
	}

	line_num = line_map.get_line_num();
	if (line_num == 0)
		return nullptr;
	Line* res = new Line[line_num];
	line_map.output_lines(res);
	return res;
}


namespace GeometryUtils_Internal
{
	class BoundaryFaceMap
	{
	public:
		inline BoundaryFaceMap() {}
		inline ~BoundaryFaceMap() {}

		inline void add_face(size_t n1, size_t n2, size_t n3)
		{
			size_t ns[3];
			ns[0] = n1;
			ns[1] = n2;
			ns[2] = n3;
			sort_array_3_acc(ns);
			char nid_str[50];
			snprintf(nid_str, sizeof(nid_str), "%zu_%zu_%zu", ns[0], ns[1], ns[2]);
			auto res = map.emplace(std::string(nid_str), FaceMapItem(n1, n2, n3));
			if (!res.second)
				++(res.first->second.elem_num);
		}

		inline size_t get_bface_num() const noexcept
		{
			size_t bface_num = 0;
			for (Map::const_iterator iter = map.begin(); iter != map.end(); ++iter)
			{
				const FaceMapItem& fi = iter->second;
				if (fi.elem_num == 0)
					++bface_num;
			}
			return bface_num;
		}

		template <typename Face>
		void output_bfaces(Face *bfaces)
		{
			Face* cur_face = bfaces;
			for (Map::iterator iter = map.begin(); iter != map.end(); ++iter)
			{
				FaceMapItem& fi = iter->second;
				if (fi.elem_num == 0)
				{
					cur_face->n1 = fi.n1;
					cur_face->n2 = fi.n2;
					cur_face->n3 = fi.n3;
					++cur_face;
				}
			}
		}

	protected:
		struct FaceMapItem
		{
			size_t n1, n2, n3, elem_num;
			inline FaceMapItem(size_t _n1, size_t _n2, size_t _n3) :
				n1(_n1), n2(_n2), n3(_n3), elem_num(0) {}
		};
		typedef std::unordered_map<std::string, FaceMapItem> Map;
		Map map;
	};
}

template <typename Triangle, typename Tetrahedron>
Triangle* extract_boundary_triangles_from_tetrahedrons(
	const Tetrahedron *elems,
	size_t elem_num,
	size_t &triangle_num
	)
{
	GeometryUtils_Internal::BoundaryFaceMap bface_map;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		const Tetrahedron& e = elems[e_id];
		bface_map.add_face(e.n1, e.n2, e.n3);
		bface_map.add_face(e.n1, e.n4, e.n2);
		bface_map.add_face(e.n2, e.n4, e.n3);
		bface_map.add_face(e.n1, e.n3, e.n4);
	}

	triangle_num = bface_map.get_bface_num();
	if (triangle_num == 0)
		return nullptr;
	
	Triangle *res = new Triangle[triangle_num];
	bface_map.output_bfaces(res);
	return res;
}

#endif