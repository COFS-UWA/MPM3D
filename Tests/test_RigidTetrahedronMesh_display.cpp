#include "Tests_pcp.h"

#include "Model_T3D_ME_s.h"
#include "QtApp_Prep_T3D_ME_s.h"

#include "test_model_view.h"

void test_RigidTetrahedronMesh_display(int argc, char** argv)
{
	Model_T3D_ME_s model;

	// load mesh
	model.load_mesh_from_hdf5("..\\..\\Asset\\brick_mesh_1.00_1x1x1.h5");
	model.init_search_grid(0.05, 0.05, 0.05);

	// init rigid tetrahedron mesh
	model.init_rb("..\\..\\Asset\\brick_mesh_1.00_1x1x1.h5", 100.0, 0.0, 0.0, 1.0);
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

	// generate particles
	ParticleGenerator3D<Model_T3D_ME_s> pcl_gen;
	pcl_gen.generate_pcls_second_order_gauss(model);
	model.init_pcls(pcl_gen, 10.0);

	QtApp_Prep_T3D_ME_s view_app(argc, argv);
	view_app.set_win_size(900, 900);
	view_app.set_view_dir(70.0, 35.0);
	view_app.set_light_dir(20.0, 20.0);
	view_app.set_model(model);
	view_app.start();
}
