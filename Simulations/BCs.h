#ifndef __BCs_h__
#define __BCs_h__

// ============ Body force ============
struct BodyForceAtPcl
{
	size_t pcl_id;
	double bf;
};

struct BodyForceAtElem
{
	size_t elem_id;
	double bf;
};

// ======== Surface traction =========
// sign: tensile - positive
struct TractionBCAtPcl
{
	size_t pcl_id;
	double t;
};

struct TractionBCAtFace
{
	size_t elem_id;
	size_t face_id;
	double t;
};

// =========== Acceleration =========== 
struct AccelerationBC
{
	size_t node_id;
	double a;
};

// ============= Velocity ============= 
struct VelocityBC
{
	size_t node_id;
	double v;
};

// ============ Displacement ========== 
struct DisplacementBC
{
	size_t node_id;
	double u;
};

#endif