#ifndef __Ratio_H__
#define __Ratio_H__

class Ratio
{
public:
	virtual double operator() (double x) const noexcept { return 0.0; }
};

class OneRatio : public Ratio
{
public:
	double operator() (double x) const noexcept override { return 1.0; }
};

class QuadraticRampUpRatio : public Ratio
{
private:
	double ramp_up_range;
public:
	QuadraticRampUpRatio() : ramp_up_range(1.0) {}
	inline void set_ramp_up_range(double ru_range)
	{
		if (ru_range != 0.0)
			ramp_up_range = ru_range;
	}
	double operator() (double x) const noexcept override
	{
		const double t = x / ramp_up_range;
		if (t < 0.0)
			return 0.0;
		
		if (t > 1.0)
			return 1.0;

		if (t < 0.5)
			return 2.0 * t * t;

		return -2.0 * t * t + 4.0 * t - 1.0;
	}
};

#endif