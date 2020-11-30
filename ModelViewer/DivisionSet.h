#ifndef __Division_Set_h__
#define __Division_Set_h__

class EmptyDivisionSet
{
public:
	EmptyDivisionSet() {}
	~EmptyDivisionSet() {}
	inline bool inside(double x, double y, double z) const noexcept
	{ return false; }
};

class PlaneDivisionSet
{
protected:
	double a, b, c, d;

public:
	PlaneDivisionSet() : a(0.0), b(0.0), c(0.0), d(-1.0) {}
	~PlaneDivisionSet() {}
	inline void set_param(double _a, double _b, double _c, double _d) noexcept
	{ _a = a; _b = b; c = _c; d = _d; }
	inline void set_by_normal_and_point(
		double nx, double ny, double nz,
		double px, double py, double pz) noexcept
	{
		a = nx; b = ny; c = nz;
		d = -(a * px + b * py + c * pz);
	}
	inline bool inside(double x, double y, double z) const noexcept
	{ return (a * x + b * y + c * z + d) >= 0.0; }
};

class BoxDivisionSet
{
protected:
	double xl, xu, yl, yu, zl, zu;

public:
	BoxDivisionSet() : xl(0.0), xu(0.0),
		yl(0.0), yu(0.0), zl(0.0), zu(0.0) {}
	~BoxDivisionSet() {}
	inline void set_param(double _xl, double _xu,
		double _yl, double _yu, double _zl, double _zu) noexcept
	{
		xl = _xl; xu = _xu;
		yl = _yl; yu = _yu;
		zl = _zl; zu = _zu;
	}
	inline bool inside(double x, double y, double z) const noexcept
	{ return !(x < xl || x > xu || y < yl || y > yu || z < zl || z > zu); }
};

template <class SetA, class SetB>
class IntersectionDivisionSet
{
protected:
	SetA a;
	SetB b;

public:
	IntersectionDivisionSet() {}
	~IntersectionDivisionSet() {}
	inline SetA& seta() noexcept { return a; }
	inline SetB& setb() noexcept { return b; }
	inline bool inside(double x, double y, double z) const noexcept
	{ return a.inside(x, y, z) && b.inside(x, y, z); }
};

typedef IntersectionDivisionSet<PlaneDivisionSet, PlaneDivisionSet> TwoPlaneDivisionSet;

#endif