#include "Tests_pcp.h"

#include "ItemArray.hpp"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_FEM_T3D_ME_s.h"
#include "Step_FEM_T3D_ME_s.h"

#include "ModelData_FEM_T3D_ME_s.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "TimeHistory_FEM_T3D_ME_s_complete.h"

#include "test_simulations.h"

#include "utils.h"
#include "test_model_view.h"

namespace
{

template <typename Point3D>
inline bool test_face_on_z_plane(
	Point3D &p1,
	Point3D &p2,
	Point3D &p3,
	double z_coord
	)
{
	static double tol = 1.0e-5;
	double zl = z_coord - tol;
	double zu = z_coord + tol;
	return p1.z > zl && p1.z < zu
		&& p2.z > zl && p2.z < zu
		&& p3.z > zl && p3.z < zu;
}

struct ElemFace
{
	size_t elem_id;
	size_t face_id;
};

void find_elem_face_on_z_plane(
	Model_FEM_T3D_ME_s &md,
	double z_coord,
	MemoryUtils::ItemArray<ElemFace> &faces)
{
	size_t elem_num = md.get_elem_num();
	Model_FEM_T3D_ME_s::Element *elems = md.get_elems();
	Model_FEM_T3D_ME_s::Node* nodes = md.get_nodes();
	ElemFace *ef;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Model_FEM_T3D_ME_s::Element& e = elems[e_id];
		Model_FEM_T3D_ME_s::Node& n1 = nodes[e.n1];
		Model_FEM_T3D_ME_s::Node& n2 = nodes[e.n2];
		Model_FEM_T3D_ME_s::Node& n3 = nodes[e.n3];
		Model_FEM_T3D_ME_s::Node& n4 = nodes[e.n4];
		if (test_face_on_z_plane(n1, n2, n3, z_coord))
		{
			ef = faces.alloc();
			ef->elem_id = e_id;
			ef->face_id = 0;
		}
		if (test_face_on_z_plane(n1, n4, n2, z_coord))
		{
			ef = faces.alloc();
			ef->elem_id = e_id;
			ef->face_id = 1;
		}
		if (test_face_on_z_plane(n1, n3, n4, z_coord))
		{
			ef = faces.alloc();
			ef->elem_id = e_id;
			ef->face_id = 2;
		}
		if (test_face_on_z_plane(n2, n4, n3, z_coord))
		{
			ef = faces.alloc();
			ef->elem_id = e_id;
			ef->face_id = 3;
		}
	}
}

}

void test_fem_t3d_me_s_1d_compression(int argc, char **argv)
{
	double mh_length = 1.0;
	size_t inv_x = 2;
	size_t inv_y = 2;
	size_t inv_z = 2;

	char mesh_filename[50];
	snprintf(mesh_filename, 50,
		"..\\..\\Asset\\brick_mesh_%.2f_%zdx%zdx%zd.h5",
		mh_length, inv_x, inv_y, inv_z);
	std::cout << mesh_filename << "\n";

	Model_FEM_T3D_ME_s model;
	model.load_mesh_from_hdf5(mesh_filename, 10.0);
	std::cout << "node num: " << model.get_node_num() << "\n"
		<< "elem num: " << model.get_elem_num() << "\n"
		<< "pcl num: " << model.get_pcl_num() << "\n";

	size_t pcl_num = model.get_pcl_num();
	Model_FEM_T3D_ME_s::Particle* pcls = model.get_pcls();
	MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_FEM_T3D_ME_s::Particle& pcl = pcls[pcl_id];
		MatModel::LinearElasticity& mm = mms[pcl_id];
		mm.set_param(100.0, 0.0);
		pcl.set_mat_model(mm);
		// stress field
		pcl.s33 = -1.0;
	}

	MemoryUtils::ItemArray<ElemFace> faces;
	faces.reserve(30);
	find_elem_face_on_z_plane(model, mh_length, faces);
	ElemFace* efs = faces.get_mem();
	model.init_tzs(faces.get_num());
	for (size_t t_id = 0; t_id < model.tz_num; ++t_id)
	{
		TractionBCAtFace& tbc = model.tzs[t_id];
		tbc.elem_id = efs[t_id].elem_id;
		tbc.face_id = efs[t_id].face_id;
		tbc.t = -1.0;
	}
	std::cout << "tz_num: " << model.tz_num << "\n";

	IndexArray pt_array(100);

	find_3d_nodes_on_x_plane(model, pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, pt_array, mh_length / inv_z * inv_x, false);
	size_t* vx_bc_n_id = pt_array.get_mem();
	model.init_vxs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC& vbc = model.vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_3d_nodes_on_y_plane(model, pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, pt_array, mh_length / inv_z * inv_x, false);
	size_t* vy_bc_n_id = pt_array.get_mem();
	model.init_vys(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC& vbc = model.vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_3d_nodes_on_z_plane(model, pt_array, 0.0);
	size_t* vz_bc_n_id = pt_array.get_mem();
	model.init_vzs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vz_num; ++v_id)
	{
		VelocityBC& vbc = model.vzs[v_id];
		vbc.node_id = vz_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	MemoryUtils::ItemArray<Point3D> ptlist(50);
	//init_tz_face_bcs_display(model, ptlist);
	//init_vx_bcs_display(model, ptlist);
	//init_vy_bcs_display(model, ptlist);
	//init_vz_bcs_display(model, ptlist);
	//display_model(argc, argv, -90.0, 20.0, 45.0, 35.0, model, ptlist, 1.0e-4, true, false);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("fem_t3d_me_s_1d_compression.h5");

	ModelData_FEM_T3D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_ConsoleProgressBar out_cpb;
	TimeHistory_FEM_T3D_ME_s_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);

	Step_FEM_T3D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0e-5);
	step.set_dtime(1.0e-5);
	//step.add_time_history(out_cpb);
	//step.add_time_history(out1);
	step.solve();
}
