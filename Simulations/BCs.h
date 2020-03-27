#ifndef __BCs_h__
#define __BCs_h__

// ============ Body force ============
struct BodyForceAtPcl
{
	double bf;
	size_t pcl_id;
};

// ======== Surface traction =========
// sign: tensile - positive
struct TractionBCAtPcl
{
	size_t pcl_id;
	double t;
};

struct TractionBCAtNode
{
	double t;
	size_t node_id;
};

// =========== Acceleration =========== 
struct AccelerationBC
{
	double a;
	size_t node_id;
};

// ============= Velocity ============= 
struct VelocityBC
{
	double v;
	size_t node_id;
};

// ============ Displacement ========== 
struct DisplacementBC
{
	double u;
	size_t node_id;
};

#endif