#ifndef __BC_H__
#define __BC_H__

struct BodyForce
{
	double bf;
	union
	{
		size_t elem_id;
		size_t pcl_id;
	};
};

// Surface traction on soil skeleton (effective stress)
// sign: tensile - positive
struct TractionBC
{
	double t;
	size_t node_id;
};

struct TractionBC_2DFEM
{
	double t0, t1;
	size_t elem_id;
	double xi0, eta0;
	double xi1, eta1;
};

struct TractionBC_MPM
{
	size_t pcl_id;
	double t;
};

// Pore pressure boundary conditions
// sign: compressive - positive
struct PressureBC
{
	double p;
	size_t node_id;
};

struct PressureBC_2D_FEM
{
	double p;
	size_t elem_id;
	size_t edge_id;
};

// Acceleration boundary conditions
struct AccelerationBC
{
	double a;
	size_t node_id;
};

// boundary condition of velocity
struct VelocityBC
{
	double v;
	size_t node_id;
};

// boundary condition of increment of displacement
struct DisplacementBC
{
	double u;
	size_t node_id;
};

#endif