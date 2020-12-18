#ifndef __Rigid_Object_By_T3D_Mesh_h__
#define __Rigid_Object_By_T3D_Mesh_h__

#include "RigidObjectMotion3D.h"
#include "RigidMeshT3D.h"

class RigidObjectByT3DMesh :
	public RigidObjectMotion3D,
	protected RigidMeshT3D
{
protected:
	double density;

public:
	explicit RigidObjectByT3DMesh();
	~RigidObjectByT3DMesh();

	int init(double _density, const char* file_name,
		double dx, double dy, double dz,
		double dx_ang, double dy_ang, double dz_ang);

	using RigidObjectMotion3D::get_local_point;
	using RigidObjectMotion3D::get_global_point;
	using RigidObjectMotion3D::get_local_vector;
	using RigidObjectMotion3D::get_global_vector;

	bool detect_collision_with_point(
		double p_x, double p_y, double p_z, double p_r,
		double& dist, Vector3D& lnorm, Point3D& lcontpos
		) const noexcept;
};

#endif