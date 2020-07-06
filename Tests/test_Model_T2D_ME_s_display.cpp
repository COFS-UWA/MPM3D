#include "Tests_pcp.h"

#include <iostream>
#include <fstream>

#include "Model_T2D_ME_s.h"
#include "QtApp_Prep_2DMPM.h"

#include "ParticleGenerator2D.hpp"

#include "test_simulations.h"

namespace
{
	typedef Model_T2D_ME_s::Element Element;
	typedef Model_T2D_ME_s::Node Node;
	typedef Model_T2D_ME_s::Particle Particle;

	void print_pcl_shape_func(
		Model_T2D_ME_s& md,
		Particle& pcl,
		std::ostream& os
		)
	{
		Node* nodes = md.get_nodes();
		pcl.pe = md.find_in_which_element(pcl);
		if (pcl.pe)
		{
			Element& e = *pcl.pe;
			Node& n1 = nodes[e.n1];
			Node& n2 = nodes[e.n2];
			Node& n3 = nodes[e.n3];
			os << "p  = Point(" << pcl.x << ", " << pcl.y << ")\n"
				<< "p1 = Point(" << n1.x << ", " << n1.y << ")\n"
				<< "p2 = Point(" << n2.x << ", " << n2.y << ")\n"
				<< "p3 = Point(" << n3.x << ", " << n3.y << ")\n"
				<< "vol: " << e.area << "\n"
				<< "N1: " << pcl.N1 << "\n"
				<< "N2: " << pcl.N2 << "\n"
				<< "N3: " << pcl.N3 << "\n"
				<< "dN1_dx: " << e.dN1_dx << "\n"
				<< "dN1_dy: " << e.dN1_dy << "\n"
				<< "dN2_dx: " << e.dN2_dx << "\n"
				<< "dN2_dy: " << e.dN2_dy << "\n"
				<< "dN3_dx: " << e.dN3_dx << "\n"
				<< "dN3_dy: " << e.dN3_dy << "\n";
		}
	}

}

void test_Model_T2D_ME_s_display(int argc, char** argv)
{
	Model_T2D_ME_s model;

	// load mesh
	double node_coords[6] = {
		0.0, 0.0, // 0
		1.0, 0.0, // 1
		0.0, 1.0  // 2
	};
	size_t node_num = sizeof(node_coords) / (sizeof(node_coords[0]) * 2);
	size_t elem_indices[] = {
		0, 1, 2 // 0
	};
	size_t elem_num = sizeof(elem_indices) / (sizeof(elem_indices[0]) * 3);
	model.init_mesh(node_coords, node_num, elem_indices, elem_num);
	//model.load_mesh_from_hdf5("..\\..\\Asset\\square_mesh.h5");

	model.init_search_grid(0.05, 0.05);
	
	elem_num = model.get_elem_num();
	Model_T2D_ME_s::Element* elems = model.get_elems();
	std::cout << "elem num: " << elem_num << "\n";
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Model_T2D_ME_s::Element &e = elems[e_id];
		std::cout << e.id << " " 
				  << e.n1 << " "
				  << e.n2 << " "
				  << e.n3 << "\n";
	}

	size_t edge_num = model.get_edge_num();
	Model_T2D_ME_s::Edge* edges = model.get_edges();
	std::cout << "edge num: " << edge_num << "\n";
	for (size_t e_id = 0; e_id < edge_num; ++e_id)
	{
		Model_T2D_ME_s::Edge &e = edges[e_id];
		std::cout << e.n1 << " " << e.n2 << "\n";
	}

	Point2D cen = model.get_centre();
	std::cout << "cen: " << cen.x << ", " << cen.y << "\n";

	Rect bbox = model.get_bounding_box();
	std::cout << "bbox: "
		<< bbox.xl << ", " << bbox.xu << ", "
		<< bbox.yl << ", " << bbox.yu << "\n";

	// generate particles
	ParticleGenerator2D<Model_T2D_ME_s> pcl_gen;
	//pcl_gen.generate_pcls_at_1st_gauss(model);
	pcl_gen.generate_pcls_at_2nd_gauss(model);
	//pcl_gen.generate_pcls_in_grid_layout(Rect(0.0f, 1.0f, 0.0f, 1.0f), 0.5f, 0.5f);
	model.init_pcls(pcl_gen, 10.0);
	std::cout << "pcl num: " << model.get_pcl_num() << "\n";

	std::fstream psf_file; // pcl shape function
	psf_file.open("pcl_shape_func.txt", std::ios::out | std::ios::binary);
	Particle* pcls = model.get_pcls();
	size_t pcl_num = model.get_pcl_num();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle& pcl = pcls[pcl_id];
		print_pcl_shape_func(model, pcl, psf_file);
	}
	psf_file.close();

	Point2D pts[2] = {
		{ 0.2, 0.2 },
		{ 0.8, 0.75 }
	};

	QtApp_Prep_2DMPM view_app(argc, argv);
	view_app.set_win_size(1000, 1000);
	view_app.set_model(model);
	//view_app.set_pts(pts, 2, 0.15);
	view_app.start();
}
