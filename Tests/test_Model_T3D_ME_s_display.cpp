#include "Tests_pcp.h"

#include "Model_T3D_ME_s.h"
#include "PrepMPM3DApp.h"

#include "test_simulations.h"

void test_Model_T3D_ME_s_display(int argc, char **argv)
{
	Model_T3D_ME_s model;

	// load mesh
	//model.load_mesh_from_hdf5("..\\..\\Asset\\teh_mesh.h5");
	//model.load_mesh_from_hdf5("..\\..\\Asset\\brick_mesh.h5");
	model.load_mesh_from_hdf5("..\\..\\Asset\\column_mesh.h5");
	model.init_edges();

	size_t elem_num = model.get_elem_num();
	Model_T3D_ME_s::Element *elems = model.get_elems();
	std::cout << "elem: " << elem_num << "\n";
	//for (size_t e_id = 0; e_id < elem_num; ++e_id)
	//{
	//	Model_T3D_ME_s::Element &e = elems[e_id];
	//	std::cout << e.id << " " << e.n1 << " " << e.n2 << " "
	//		<< e.n3 << " " << e.n4 << "\n";
	//}

	size_t edge_num = model.get_edge_num();
	Model_T3D_ME_s::Edge *edges = model.get_edges();
	std::cout << "edge: " << edge_num << "\n";
	//for (size_t e_id = 0; e_id < edge_num; ++e_id)
	//{
	//	Model_T3D_ME_s::Edge &e = edges[e_id];
	//	std::cout << e.n1 << " " << e.n2 << "\n";
	//}

	Point3D cen = model.get_centre();
	std::cout << "cen: " << cen.x << ", " << cen.y << ", " << cen.z << "\n";

	Cube bbox = model.get_bounding_box();
	std::cout << "bbox: "
		<< bbox.xl << " " << bbox.xu << " "
		<< bbox.yl << " " << bbox.yu << " "
		<< bbox.zl << " " << bbox.zu << "\n";

	// generate particles
	ParticleGenerator3D<Model_T3D_ME_s> pcl_gen;
	pcl_gen.generate_pcls_first_order_gauss(model);
	model.init_pcls(pcl_gen, 10.0);
	std::cout << model.get_pcl_num() << "\n";

	PrepMPM3DApp view_app(argc, argv);
	view_app.set_view_dir(-1.0f, -2.0f, -1.0f);
	view_app.set_model(model);
	view_app.start();
}
