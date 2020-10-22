#ifndef __Triangle_Utils_f_h__
#define __Triangle_Utils_f_h__

template <typename Node2D, typename Point2D>
inline float cal_triangle_area_f(
	const Node2D& p1,
	const Node2D& p2,
	const Point2D& p3
	) noexcept
{
	return 0.5f * ((p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y));
}

struct PointInTriangle_f
{
protected:
	float a1, b1, coef1;
	float a2, b2, coef2;
	float a3, b3, coef3;

public:
	template <typename Node2D>
	inline void init_triangle(
		const Node2D& n1,
		const Node2D& n2,
		const Node2D& n3,
		float area
	) noexcept
	{
		float area2 = area * 2.0f;
		a1 = (n2.y - n3.y) / area2;
		b1 = (n3.x - n2.x) / area2;
		coef1 = (n2.x * n3.y - n3.x * n2.y) / area2;
		a2 = (n3.y - n1.y) / area2;
		b2 = (n1.x - n3.x) / area2;
		coef2 = (n3.x * n1.y - n1.x * n3.y) / area2;
		a3 = (n1.y - n2.y) / area2;
		b3 = (n2.x - n1.x) / area2;
		coef3 = (n1.x * n2.y - n2.x * n1.y) / area2;
	}

	template <typename Node2D>
	inline void init_triangle(
		const Node2D& n1,
		const Node2D& n2,
		const Node2D& n3
	) noexcept
	{
		init_triangle(n1, n2, n3, cal_triangle_area_f<Node2D, Node2D>(n1, n2, n3));
	}

	inline bool is_in_triangle(float x, float y) const noexcept
	{
		float N1v = N1(x, y);
		float N2v = N2(x, y);
		float N3v = N3(x, y);
		return !(N1v < 0.0 || N1v > 1.0 ||
				 N2v < 0.0 || N2v > 1.0 ||
				 N3v < 0.0 || N3v > 1.0);
	}
	template <typename Point2D>
	inline bool is_in_triangle(const Point2D& p) const noexcept
	{
		return is_in_triangle(p.x, p.y);
	}

	// shape functions
	inline float N1(float x, float y) const noexcept
	{
		return a1 * x + b1 * y - coef1;
	}
	inline float N2(float x, float y) const noexcept
	{
		return a2 * x + b2 * y - coef2;
	}
	inline float N3(float x, float y) const noexcept
	{
		return a3 * x + b3 * y - coef3;
	}
	inline void cal_N(float x, float y,
		float& N1v, float& N2v, float& N3v) const noexcept
	{
		N1v = N1(x, y); N2v = N2(x, y); N3v = 1.0 - N1v - N2v;
	}

	// shape function derivatives
	inline float dN1_dx() const noexcept { return a1; }
	inline float dN1_dy() const noexcept { return b1; }
	inline float dN2_dx() const noexcept { return a2; }
	inline float dN2_dy() const noexcept { return b2; }
	inline float dN3_dx() const noexcept { return a3; }
	inline float dN3_dy() const noexcept { return b3; }

	inline float get_a1() const noexcept { return a1; }
	inline float get_b1() const noexcept { return b1; }
	inline float get_coef1() const noexcept { return coef1; }
	inline float get_a2() const noexcept { return a2; }
	inline float get_b2() const noexcept { return b2; }
	inline float get_coef2() const noexcept { return coef2; }
	inline float get_a3() const noexcept { return a3; }
	inline float get_b3() const noexcept { return b3; }
	inline float get_coef3() const noexcept { return coef3; }
};

#endif