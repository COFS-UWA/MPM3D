#include "Simulations_pcp.h"

#include "Model_T3D_ME_s.h"
#include "Model_T3D_ME_s_hdf5_utilities.h"

#include "ModelData_T3D_ME_s.h"

int model_data_output_func_t3d_me_s_to_xml_res_file(ModelData &_self)
{
	ModelData_T3D_ME_s &md = static_cast<ModelData_T3D_ME_s &>(_self);
	Model_T3D_ME_s &model = static_cast<Model_T3D_ME_s &>(*md.model);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(*md.res_file);
	std::fstream &file = rf.get_file();

	char str_buffer[512];

	// model data
	const char *model_data_info = ""
		"<ModelData>\n";
	file.write(model_data_info, strlen(model_data_info));

	// mesh
	size_t node_num = model.get_node_num();
	Model_T3D_ME_s::Node *nodes = model.get_nodes();
	size_t elem_num = model.get_elem_num();
	Model_T3D_ME_s::Element *elems = model.get_elems();

	const char *mesh_info1 = ""
		"    <BackGroundMesh type = \"T3D\">\n"
		"        <Nodes num=%zu>\n"
		"            <coordinates>\n"
		"            <!-- id, x, y, z -->\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info1, node_num);
	file.write(str_buffer, strlen(str_buffer));
	const char *mesh_coords = "            %zu, %le, %le, %le\n";
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Model_T3D_ME_s::Node &n = nodes[n_id];
		snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]),
				 mesh_coords, n.id, n.x, n.y, n.z);
		file.write(str_buffer, strlen(str_buffer));
	}

	const char *mesh_info2 = ""
		"            </coordinates>\n"
		"        </Nodes>\n"
		"        <Elements num=%zu>\n"
		"            <node_indices>\n";
		"            <!-- id, n1, n2, n3, n4 -->\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]), mesh_info2, elem_num);
	file.write(str_buffer, strlen(str_buffer));
	const char *mesh_indices = "            %zu, %zu, %zu, %zu, %zu\n";
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Model_T3D_ME_s::Element &e = elems[e_id];
		snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]),
				 mesh_indices, e.id, e.n1, e.n2, e.n3, e.n4);
		file.write(str_buffer, strlen(str_buffer));
	}
	const char *mesh_info3 = ""
		"            </node_indices>\n"
		"        </Elements>\n"
		"    </BackGroundMesh>\n";
	file.write(mesh_info3, strlen(mesh_info3));

	// material point object
	const char *mps_info = ""
		"    <MaterialPoints type = \"ME_3D\" num=%zu>\n"
		"    </MaterialPoints>\n";
	snprintf(str_buffer, sizeof(str_buffer) / sizeof(str_buffer[0]),
			 mps_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));

	// ending
	file << "</ModelData>\n";
	
	return 0;
}

int model_data_output_func_t3d_me_s_to_hdf5_res_file(ModelData &_self)
{
	ModelData_T3D_ME_s &md = static_cast<ModelData_T3D_ME_s &>(_self);
	Model_T3D_ME_s &model = static_cast<Model_T3D_ME_s &>(*md.model);
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*md.res_file);
	Model_T3D_ME_s_hdf5_utilities::output_model_to_hdf5_file(model, rf);
	return 0;
}
