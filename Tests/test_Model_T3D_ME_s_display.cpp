#include "Tests_pcp.h"

#include <iostream>
#include <fstream>

#include "Model_T3D_ME_s.h"
#include "QtApp_Prep_T3D_ME_s.h"

#include "test_simulations.h"

namespace
{

typedef Model_T3D_ME_s::Element Element;
typedef Model_T3D_ME_s::Node Node;
typedef Model_T3D_ME_s::Particle Particle;

void print_pcl_shape_func(
	Model_T3D_ME_s &md,
	Particle& pcl,
	std::ostream &os
	)
{
	Node* nodes = md.get_nodes();
	pcl.pe = md.find_in_which_element(pcl);
	if (pcl.pe)
	{
		Element &e = *pcl.pe;
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		Node& n4 = nodes[e.n4];
		os << "p  = Point(" << pcl.x << ", " << pcl.y << ", " << pcl.z << ")\n"
		   << "p1 = Point(" << n1.x << ", " << n1.y << ", " << n1.z << ")\n"
		   << "p2 = Point(" << n2.x << ", " << n2.y << ", " << n2.z << ")\n"
		   << "p3 = Point(" << n3.x << ", " << n3.y << ", " << n3.z << ")\n"
		   << "p4 = Point(" << n4.x << ", " << n4.y << ", " << n4.z << ")\n"
			<< "vol: " << e.vol << "\n"
			<< "N1: " << pcl.N1 << "\n"
			<< "N2: " << pcl.N2 << "\n"
			<< "N3: " << pcl.N3 << "\n"
			<< "N4: " << pcl.N4 << "\n"
			<< "dN1_dx: " << e.dN1_dx << "\n"
			<< "dN1_dy: " << e.dN1_dy << "\n"
			<< "dN1_dz: " << e.dN1_dz << "\n"
			<< "dN2_dx: " << e.dN2_dx << "\n"
			<< "dN2_dy: " << e.dN2_dy << "\n"
			<< "dN2_dz: " << e.dN2_dz << "\n"
			<< "dN3_dx: " << e.dN3_dx << "\n"
			<< "dN3_dy: " << e.dN3_dy << "\n"
			<< "dN3_dz: " << e.dN3_dz << "\n"
			<< "dN4_dx: " << e.dN4_dx << "\n"
			<< "dN4_dy: " << e.dN4_dy << "\n"
			<< "dN4_dz: " << e.dN4_dz << "\n";
	}
}

}

void test_Model_T3D_ME_s_display(int argc, char **argv)
{
	Model_T3D_ME_s model;

	// load mesh
	//model.load_mesh_from_hdf5("../..\\Asset/teh_mesh.h5");
	//model.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_1x1x1.h5");
	//model.load_mesh_from_hdf5("../../Asset/brick_mesh_plus.h5");
	model.load_mesh_from_hdf5("../../Asset/ball_mesh_r1.00.h5");

	size_t elem_num = model.get_elem_num();
	Model_T3D_ME_s::Element *elems = model.get_elems();
	std::cout << "elem num: " << elem_num << "\n";
	//for (size_t e_id = 0; e_id < elem_num; ++e_id)
	//{
	//	Model_T3D_ME_s::Element &e = elems[e_id];
	//	std::cout << e.id << " " << e.n1 << " " << e.n2 << " "
	//		<< e.n3 << " " << e.n4 << "\n";
	//}

	size_t edge_num = model.get_edge_num();
	Model_T3D_ME_s::Edge *edges = model.get_edges();
	std::cout << "edge num: " << edge_num << "\n";
	//for (size_t e_id = 0; e_id < edge_num; ++e_id)
	//{
	//	Model_T3D_ME_s::Edge &e = edges[e_id];
	//	std::cout << e.n1 << " " << e.n2 << "\n";
	//}

	Point3D cen = model.get_centre();
	std::cout << "cen: " << cen.x << ", " << cen.y << ", " << cen.z << "\n";

	Cube bbox = model.get_bounding_box();
	std::cout << "bbox: "
		<< bbox.xl << ", " << bbox.xu << ", "
		<< bbox.yl << ", " << bbox.yu << ", "
		<< bbox.zl << ", " << bbox.zu << "\n";

	model.init_search_grid(0.05, 0.05, 0.05);

	// generate particles
	ParticleGenerator3D<Model_T3D_ME_s> pcl_gen;
	pcl_gen.generate_pcls_second_order_gauss(model);
	//pcl_gen.generate_pcls_first_order_gauss(model);
	model.init_pcls(pcl_gen, 10.0);
	std::cout << "pcl num: " << model.get_pcl_num() << "\n";

	QtApp_Prep_T3D_ME_s view_app(argc, argv);
	view_app.set_win_size(900, 900);
	view_app.set_view_dir(70.0, 35.0);
	view_app.set_light_dir(20.0, 20.0);
	view_app.set_model(model);
	//view_app.set_display_pcls(false);
	view_app.start();;
}
