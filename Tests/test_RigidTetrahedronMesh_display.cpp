#include "Tests_pcp.h"

#include <Eigen/Dense>

#include "Model_T3D_ME_s.h"
#include "QtApp_Prep_T3D_ME_s.h"

#include "test_model_view.h"

void test_RigidTetrahedronMesh_display(int argc, char** argv)
{
	Model_T3D_ME_s model;

	// load mesh
	model.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_1x1x1.h5");
	model.init_search_grid(0.05, 0.05, 0.05);

	// init rigid tetrahedron mesh
	model.init_rb(1.0, "../../Asset/brick_mesh_1.00_1x1x1.h5", 0.0, 0.0, 1.0, 0.7854);
	model.set_contact_params(100.0, 0.0, 0.0);
	
	RigidTetrahedronMesh& rb = model.get_rb();
	
	Point3D cen = rb.get_centre();
	std::cout << "centre: " << cen.x << ", " << cen.y << ", " << cen.z << "\n";
	
	RigidTetrahedronMesh::Face* bfaces = rb.get_bfaces();
	size_t bface_num = rb.get_bface_num();
	std::cout << "faces: " << bface_num << "\n";
	for (size_t f_id = 0; f_id < bface_num; ++f_id)
	{
		RigidTetrahedronMesh::Face& f = bfaces[f_id];
		std::cout << f.n1 << ", " << f.n2 << ", " << f.n3 << "\n";
	}

	rb.init_bg_grids(0.5, 0.3);

	Cube rb_bbox;
	rb.get_cur_bbox(rb_bbox);
	std::cout << "bbox: " << rb_bbox.xl << ", " << rb_bbox.xu
				  << ", " << rb_bbox.yl << ", " << rb_bbox.yu
				  << ", " << rb_bbox.zl << ", " << rb_bbox.zu << "\n";

	// generate particles
	ParticleGenerator3D<Model_T3D_ME_s> pcl_gen;
	pcl_gen.generate_pcls_second_order_gauss(model);
	model.init_pcls(pcl_gen, 10.0);

	QtApp_Prep_T3D_ME_s view_app(argc, argv);
	view_app.set_win_size(950, 950);
	view_app.set_view_dir(0.0, 0.0);
	view_app.set_view_dist_scale(2.0f);
	view_app.set_light_dir(20.0, 20.0);
	view_app.set_rb_display_mode(QtRigidTetrahedronMeshGLObject::LineFrame);
	view_app.set_model(model);
	view_app.start();
}

void test_RigidTetrahedronMesh_moi(int argc, char** argv)
{
	double rb_m;
	typedef Eigen::Matrix<double, 3, 3> Matrix3x3;

	//RigidTetrahedronMesh rb1;
	//rb1.init_mesh("../../Asset/teh_mesh.h5", 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0);
	//rb_m = rb1.get_m();
	//const Matrix3x3& rb_mat1 = rb1.get_moi();
	//std::cout << "m: " << rb_m << "\nmoi:\n" << rb_mat1 << "\n\n";

	//RigidTetrahedronMesh rb2;
	//rb2.init_mesh("../../Asset/brick_mesh_1.00_1x1x1.h5", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//rb_m = rb2.get_m();
	//const Matrix3x3& rb_mat2 = rb2.get_moi();
	//std::cout << "m: " << rb_m << "\nmoi:\n" << rb_mat2 << "\n\n";

	RigidTetrahedronMesh rb3;
	rb3.init_mesh("../../Asset/brick_mesh_1.00_2x2x10.h5", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	rb_m = rb3.get_m();
	const Matrix3x3& rb_mat3 = rb3.get_moi();
	std::cout << "m: " << rb_m << "\nmoi:\n" << rb_mat3 << "\n\n";

	RigidTetrahedronMesh rb4;
	rb4.init_mesh("../../Asset/ball_mesh_r1.00.h5", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	rb_m = rb4.get_m();
	const Matrix3x3& rb_mat4 = rb4.get_moi();
	std::cout << "m: " << rb_m << "\nmoi:\n" << rb_mat4 << "\n\n";
}