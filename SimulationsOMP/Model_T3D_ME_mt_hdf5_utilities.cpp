#include "SimulationsOMP_pcp.h"

#include <exception>

#include "MatModelIdToPointerMap.h"
#include "Model_hdf5_utilities.h"
#include "Model_T3D_ME_mt_hdf5_utilities.h"
#include "RigidObject/RigidObject_hdf5_utilities.h"

namespace Model_T3D_ME_mt_hdf5_utilities
{
using namespace Model_hdf5_utilities;

int output_background_mesh_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bg_mesh_grp_id = rf.create_group(grp_id, "BackgroundMesh");

	size_t node_num = md.get_node_num();
	size_t elem_num = md.get_elem_num();

	// bg mesh properties
	const char* bg_mesh_type = "T2D_mt";
	rf.write_attribute(bg_mesh_grp_id, "type", strlen(bg_mesh_type), bg_mesh_type);
	rf.write_attribute(bg_mesh_grp_id, "node_num", node_num);
	rf.write_attribute(bg_mesh_grp_id, "element_num", elem_num);

	// bg grid properties
	rf.write_attribute(bg_mesh_grp_id, "bg_grid_hx", md.get_bg_grid_hx());
	rf.write_attribute(bg_mesh_grp_id, "bg_grid_hy", md.get_bg_grid_hy());
	rf.write_attribute(bg_mesh_grp_id, "bg_grid_hz", md.get_bg_grid_hz());

	// node coordinates
	NodeData *nodes_data = new NodeData[node_num];
	Model_T3D_ME_mt::Position *node_pos = md.node_pos;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		NodeData &node_data = nodes_data[n_id];
		node_data.id = n_id;
		Model_T3D_ME_mt::Position &np = node_pos[n_id];
		node_data.x = np.x;
		node_data.y = np.y;
		node_data.z = np.z;
	}
	hid_t nd_dt_id = get_node_dt_id();
	rf.write_dataset(
		bg_mesh_grp_id,
		"NodeData",
		node_num,
		nodes_data,
		nd_dt_id
		);
	H5Tclose(nd_dt_id);
	delete[] nodes_data;

	// element indices
	ElementData *elems_data = new ElementData[elem_num];
	Model_T3D_ME_mt::ElemNodeIndex *e_node_index = md.elem_node_id;
	double* e_vol = md.elem_vol;
	Model_T3D_ME_mt::DShapeFuncABC* e_dN_abc = md.elem_dN_abc;
	Model_T3D_ME_mt::DShapeFuncD* e_dN_d = md.elem_dN_d;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		ElementData &elem_data = elems_data[e_id];
		elem_data.id = e_id;
		Model_T3D_ME_mt::ElemNodeIndex &eni = e_node_index[e_id];
		elem_data.n1 = eni.n1;
		elem_data.n2 = eni.n2;
		elem_data.n3 = eni.n3;
		elem_data.n4 = eni.n4;
		elem_data.vol = e_vol[e_id];
		Model_T3D_ME_mt::DShapeFuncABC& edNabc = e_dN_abc[e_id];
		elem_data.a1 = edNabc.a1;
		elem_data.b1 = edNabc.b1;
		elem_data.c1 = edNabc.c1;
		elem_data.a2 = edNabc.a2;
		elem_data.b2 = edNabc.b2;
		elem_data.c2 = edNabc.c2;
		elem_data.a3 = edNabc.a3;
		elem_data.b3 = edNabc.b3;
		elem_data.c3 = edNabc.c3;
		elem_data.a4 = edNabc.a4;
		elem_data.b4 = edNabc.b4;
		elem_data.c4 = edNabc.c4;
		Model_T3D_ME_mt::DShapeFuncD &edNd = e_dN_d[e_id];
		elem_data.d1 = edNd.d1;
		elem_data.d2 = edNd.d2;
		elem_data.d3 = edNd.d3;
		elem_data.d4 = edNd.d4;
	}
	hid_t ed_dt_id = get_element_dt_id();
	rf.write_dataset(
		bg_mesh_grp_id,
		"ElementData",
		elem_num,
		elems_data,
		ed_dt_id
		);
	H5Tclose(ed_dt_id);
	delete[] elems_data;

	rf.close_group(bg_mesh_grp_id);
	return 0;
}

int load_background_mesh_from_hdf5_file(
	Model_T3D_ME_mt &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bg_mesh_grp_id = rf.open_group(grp_id, "BackgroundMesh");

	size_t elem_num, node_num;
	rf.read_attribute(bg_mesh_grp_id, "node_num", node_num);
	rf.read_attribute(bg_mesh_grp_id, "element_num", elem_num);
	md.alloc_mesh(node_num, elem_num);

	// nodes
	NodeData *nodes_data = new NodeData[node_num];
	hid_t nd_dt_id = get_node_dt_id();
	rf.read_dataset(
		bg_mesh_grp_id,
		"NodeData",
		node_num,
		nodes_data,
		nd_dt_id
		);
	H5Tclose(nd_dt_id);
	Model_T3D_ME_mt::Position *node_pos = md.node_pos;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Model_T3D_ME_mt::Position& np = node_pos[n_id];
		NodeData &node_data = nodes_data[n_id];
		np.x = node_data.x;
		np.y = node_data.y;
		np.z = node_data.z;
	}
	delete[] nodes_data;

	// elements
	ElementData *elems_data = new ElementData[elem_num];
	hid_t ed_dt_id = get_element_dt_id();
	rf.read_dataset(
		bg_mesh_grp_id,
		"ElementData",
		elem_num,
		elems_data,
		ed_dt_id
		);
	H5Tclose(ed_dt_id);
	double *elem_vol = md.elem_vol;
	Model_T3D_ME_mt::ElemNodeIndex* e_node_id = md.elem_node_id;
	Model_T3D_ME_mt::DShapeFuncABC* e_dN_abc = md.elem_dN_abc;
	Model_T3D_ME_mt::DShapeFuncD* e_dN_d = md.elem_dN_d;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		ElementData &elem_data = elems_data[e_id];
		Model_T3D_ME_mt::ElemNodeIndex& eni = e_node_id[e_id];
		eni.n1 = elem_data.n1;
		eni.n2 = elem_data.n2;
		eni.n3 = elem_data.n3;
		eni.n4 = elem_data.n4;
		elem_vol[e_id] = elem_data.vol;
		Model_T3D_ME_mt::DShapeFuncABC &edNabc = e_dN_abc[e_id];
		edNabc.a1 = elem_data.a1;
		edNabc.a2 = elem_data.a2;
		edNabc.a3 = elem_data.a3;
		edNabc.a4 = elem_data.a4;
		edNabc.b1 = elem_data.b1;
		edNabc.b2 = elem_data.b2;
		edNabc.b3 = elem_data.b3;
		edNabc.b4 = elem_data.b4;
		edNabc.c1 = elem_data.c1;
		edNabc.c2 = elem_data.c2;
		edNabc.c3 = elem_data.c3;
		edNabc.c4 = elem_data.c4;
		Model_T3D_ME_mt::DShapeFuncD& edNd = e_dN_d[e_id];
		edNd.d1 = elem_data.d1;
		edNd.d2 = elem_data.d2;
		edNd.d3 = elem_data.d3;
		edNd.d4 = elem_data.d4;
	}
	delete[] elems_data;
	
	rf.close_group(bg_mesh_grp_id);
	return 0;
}

int output_boundary_condition_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0 || md.pcl_num == 0)
		return -1;

	hid_t bc_grp_id = rf.create_group(grp_id, "BoundaryCondition");
	
	rf.write_attribute(bc_grp_id, "node_num", md.node_num);

	hid_t n_vbc_dt_id = get_node_vbc_dt_id();
	NodeVBCData* nvbc_data = new NodeVBCData[md.node_num];
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		NodeVBCData& nvbc = nvbc_data[n_id];
		Model_T3D_ME_mt::NodeHasVBC& md_nvbc = md.node_has_vbc[n_id];
		nvbc.has_vx_bc = md_nvbc.has_vx_bc;
		nvbc.has_vy_bc = md_nvbc.has_vy_bc;
		nvbc.has_vz_bc = md_nvbc.has_vz_bc;
	}
	rf.write_dataset(bc_grp_id, "NodeVelocityBC", md.node_num, nvbc_data, n_vbc_dt_id);
	delete[] nvbc_data;
	H5Tclose(n_vbc_dt_id);

	hid_t n_vbc_vec_dt_id = get_node_vbc_vec_dt_id();
	NodeVBCVecData* nvbc_vec = new NodeVBCVecData[md.node_num];
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		NodeVBCVecData& nvbc = nvbc_vec[n_id];
		Model_T3D_ME_mt::NodeVBCVec& md_nvbc = md.node_vbc_vec[n_id];
		nvbc.x = md_nvbc.x;
		nvbc.y = md_nvbc.y;
		nvbc.z = md_nvbc.z;
	}
	rf.write_dataset(bc_grp_id, "NodeVBCVector", md.node_num, nvbc_vec, n_vbc_vec_dt_id);
	delete[] nvbc_vec;
	H5Tclose(n_vbc_vec_dt_id);
	
	rf.close_group(bc_grp_id);
	return 0;
}

int load_boundary_condition_from_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;
	
	hid_t md_id = rf.get_model_data_grp_id();

	hid_t bc_id = rf.open_group(md_id, "BoundaryCondition");
	
	rf.read_attribute(bc_id, "node_num", md.node_num);
	
	hid_t n_vbc_dt_id = get_node_vbc_dt_id();
	NodeVBCData *nvbc_data = new NodeVBCData[md.node_num];
	rf.read_dataset(bc_id, "NodeVelocityBC", md.node_num, nvbc_data, n_vbc_dt_id);
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		NodeVBCData& nvbc = nvbc_data[n_id];
		Model_T3D_ME_mt::NodeHasVBC &md_nvbc = md.node_has_vbc[n_id];
		md_nvbc.has_vx_bc = nvbc.has_vx_bc;
		md_nvbc.has_vy_bc = nvbc.has_vy_bc;
		md_nvbc.has_vz_bc = nvbc.has_vz_bc;
	}
	delete[] nvbc_data;
	H5Tclose(n_vbc_dt_id);

	hid_t n_vbc_vec_dt_id = get_node_vbc_vec_dt_id();
	NodeVBCVecData* nvbc_vec = new NodeVBCVecData[md.node_num];
	rf.read_dataset(bc_id, "NodeVBCVector", md.node_num, nvbc_vec, n_vbc_vec_dt_id);
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		NodeVBCVecData& nvbc = nvbc_vec[n_id];
		Model_T3D_ME_mt::NodeVBCVec& md_nvbc = md.node_vbc_vec[n_id];
		md_nvbc.x = nvbc.x;
		md_nvbc.y = nvbc.y;
		md_nvbc.z = nvbc.z;
	}
	delete[] nvbc_vec;
	H5Tclose(n_vbc_vec_dt_id);
	
	rf.close_group(bc_id);
	return 0;
}

int output_search_mesh_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t search_mesh_grp_id = rf.create_group(grp_id, "SearchMesh");

	rf.write_attribute(search_mesh_grp_id, "grid_xl", md.grid_xl);
	rf.write_attribute(search_mesh_grp_id, "grid_yl", md.grid_yl);
	rf.write_attribute(search_mesh_grp_id, "grid_zl", md.grid_zl);
	rf.write_attribute(search_mesh_grp_id, "grid_xu", md.grid_xu);
	rf.write_attribute(search_mesh_grp_id, "grid_yu", md.grid_yu);
	rf.write_attribute(search_mesh_grp_id, "grid_zu", md.grid_zu);
	rf.write_attribute(search_mesh_grp_id, "grid_hx", md.grid_hx);
	rf.write_attribute(search_mesh_grp_id, "grid_hy", md.grid_hy);
	rf.write_attribute(search_mesh_grp_id, "grid_hz", md.grid_hz);
	rf.write_attribute(search_mesh_grp_id, "grid_x_num", md.grid_x_num);
	rf.write_attribute(search_mesh_grp_id, "grid_y_num", md.grid_y_num);
	rf.write_attribute(search_mesh_grp_id, "grid_z_num", md.grid_z_num);

	int res;
	size_t grid_elem_list_len = md.grid_x_num * md.grid_y_num * md.grid_z_num + 1;
	rf.write_attribute(
		search_mesh_grp_id,
		"grid_elem_list_len",
		grid_elem_list_len
		);
	res = rf.write_dataset(
		search_mesh_grp_id,
		"grid_elem_list",
		grid_elem_list_len,
		(unsigned long long *)md.grid_elem_list
		);

	size_t grid_elem_list_id_len = md.grid_elem_list[grid_elem_list_len-1];
	rf.write_attribute(
		search_mesh_grp_id,
		"grid_elem_list_id_len",
		grid_elem_list_id_len
		);
	res = rf.write_dataset(
		search_mesh_grp_id,
		"grid_elem_list_id",
		grid_elem_list_id_len,
		(unsigned long long *)md.grid_elem_list_id_array
		);

	rf.close_group(search_mesh_grp_id);
	return res;
}	

int load_search_mesh_from_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	md.clear_search_grid();

	hid_t search_mesh_grp_id = rf.open_group(grp_id, "SearchMesh");

	rf.read_attribute(search_mesh_grp_id, "grid_xl", md.grid_xl);
	rf.read_attribute(search_mesh_grp_id, "grid_yl", md.grid_yl);
	rf.read_attribute(search_mesh_grp_id, "grid_zl", md.grid_zl);
	rf.read_attribute(search_mesh_grp_id, "grid_xu", md.grid_xu);
	rf.read_attribute(search_mesh_grp_id, "grid_yu", md.grid_yu);
	rf.read_attribute(search_mesh_grp_id, "grid_zu", md.grid_zu);
	rf.read_attribute(search_mesh_grp_id, "grid_hx", md.grid_hx);
	rf.read_attribute(search_mesh_grp_id, "grid_hy", md.grid_hy);
	rf.read_attribute(search_mesh_grp_id, "grid_hz", md.grid_hz);
	rf.read_attribute(search_mesh_grp_id, "grid_x_num", md.grid_x_num);
	rf.read_attribute(search_mesh_grp_id, "grid_y_num", md.grid_y_num);
	rf.read_attribute(search_mesh_grp_id, "grid_z_num", md.grid_z_num);
	md.grid_xy_num = md.grid_x_num * md.grid_y_num;

	size_t grid_elem_list_len;
	rf.read_attribute(
		search_mesh_grp_id,
		"grid_elem_list_len",
		grid_elem_list_len
		);
	md.grid_elem_list = new size_t[grid_elem_list_len];
	rf.read_dataset(
		search_mesh_grp_id,
		"grid_elem_list",
		grid_elem_list_len,
		md.grid_elem_list);

	size_t grid_elem_list_id_len;
	rf.read_attribute(
		search_mesh_grp_id,
		"grid_elem_list_id_len",
		grid_elem_list_id_len);
	md.grid_elem_list_id_array = new size_t[grid_elem_list_id_len];
	rf.read_dataset(
		search_mesh_grp_id,
		"grid_elem_list_id",
		grid_elem_list_id_len,
		md.grid_elem_list_id_array);

	rf.close_group(search_mesh_grp_id);
	return 0;
}

int output_ori_pcl_data_to_hdf5_file(
	Model_T3D_ME_mt &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.create_group(grp_id, "ParticleData");

	size_t ori_pcl_num = md.get_ori_pcl_num();
	rf.write_attribute(pcl_data_grp_id, "ori_pcl_num", ori_pcl_num);
	size_t pcl_num = md.get_pcl_num();
	rf.write_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

	ParticleData* pcl_data = new ParticleData[pcl_num];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		ParticleData& pd = pcl_data[p_id];
		pd.from_pcl(md, 0, p_id);
	}
	hid_t pcl_dt_id = get_pcl_dt_id();
	int res = rf.write_dataset(
		pcl_data_grp_id,
		"field", pcl_num,
		pcl_data, pcl_dt_id);
	H5Tclose(pcl_dt_id);
	delete[] pcl_data;

	rf.close_group(pcl_data_grp_id);
	return res;
}

int output_pcl_data_to_hdf5_file(
	Model_T3D_ME_mt &md,
	Step_T3D_ME_mt &stp,
	ResultFile_hdf5 &rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.create_group(grp_id, "ParticleData");

	size_t ori_pcl_num = md.get_ori_pcl_num();
	rf.write_attribute(pcl_data_grp_id, "ori_pcl_num", ori_pcl_num);
	size_t pcl_num = stp.get_pcl_num();
	rf.write_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

	int res = 0;
	if (pcl_num)
	{
		size_t sorted_pcl_var_id = stp.get_sorted_pcl_var_id();
		ParticleData* pcl_data = new ParticleData[pcl_num];
		for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		{
			ParticleData& pd = pcl_data[p_id];
			pd.from_pcl(stp, sorted_pcl_var_id, p_id);
		}
		hid_t pcl_dt_id = get_pcl_dt_id();
		res = rf.write_dataset(
			pcl_data_grp_id,
			"field",
			pcl_num,
			pcl_data,
			pcl_dt_id);
		H5Tclose(pcl_dt_id);
		delete[] pcl_data;
	}

	rf.close_group(pcl_data_grp_id);
	return res;
}

int output_pcl_data_to_hdf5_file(
	Model_T3D_ME_mt& md,
	Step_T3D_ME_TBB& stp,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.create_group(grp_id, "ParticleData");

	const size_t ori_pcl_num = md.get_ori_pcl_num();
	rf.write_attribute(pcl_data_grp_id, "ori_pcl_num", ori_pcl_num);

	const size_t pcl_num = stp.prev_valid_pcl_num;
	rf.write_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

	int res = 0;
	if (pcl_num)
	{
		ParticleData* pcl_data = new ParticleData[pcl_num];
		for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		{
			ParticleData& pd = pcl_data[p_id];
			pd.from_pcl(stp, p_id);
		}
		hid_t pcl_dt_id = get_pcl_dt_id();
		res = rf.write_dataset(
			pcl_data_grp_id,
			"field",
			pcl_num,
			pcl_data,
			pcl_dt_id);
		H5Tclose(pcl_dt_id);
		delete[] pcl_data;
	}

	rf.close_group(pcl_data_grp_id);
	return res;
}

int output_pcl_data_to_hdf5_file(
	Model_T3D_ME_mt& md, 
	Step_T3D_ME_mt_Geo& stp,
	ResultFile_hdf5& rf, 
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.create_group(grp_id, "ParticleData");

	size_t ori_pcl_num = md.get_ori_pcl_num();
	rf.write_attribute(pcl_data_grp_id, "ori_pcl_num", ori_pcl_num);
	size_t pcl_num = stp.get_pcl_num();
	rf.write_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

	int res = 0;
	if (pcl_num)
	{
		size_t sorted_pcl_var_id = stp.get_sorted_pcl_var_id();
		ParticleData* pcl_data = new ParticleData[pcl_num];
		for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		{
			ParticleData& pd = pcl_data[p_id];
			pd.from_pcl(stp, sorted_pcl_var_id, p_id);
		}
		hid_t pcl_dt_id = get_pcl_dt_id();
		res = rf.write_dataset(
			pcl_data_grp_id,
			"field",
			pcl_num,
			pcl_data,
			pcl_dt_id);
		H5Tclose(pcl_dt_id);
		delete[] pcl_data;
	}

	rf.close_group(pcl_data_grp_id);
	return res;
}

int load_pcl_data_from_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.open_group(grp_id, "ParticleData");

	size_t ori_pcl_num, pcl_num;
	rf.read_attribute(pcl_data_grp_id, "ori_pcl_num", ori_pcl_num);
	rf.read_attribute(pcl_data_grp_id, "pcl_num", pcl_num);
	md.ori_pcl_num = ori_pcl_num;
	md.pcl_num = pcl_num;

	if (pcl_num && rf.has_dataset(pcl_data_grp_id, "field"))
	{
		ParticleData* pcls_data = new ParticleData[pcl_num];
		hid_t pcl_dt_id = get_pcl_dt_id();
		int res = rf.read_dataset(
			pcl_data_grp_id,
			"field",
			pcl_num,
			pcls_data,
			pcl_dt_id
			);
		H5Tclose(pcl_dt_id);

		if (res)
		{
			delete[] pcls_data;
			rf.close_group(pcl_data_grp_id);
			return res;
		}

		MatModel::MatModelIdToPointerMap mm_id_map(md);
		
		md.alloc_pcls(pcl_num, ori_pcl_num);
		for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		{
			ParticleData& pcl_data = pcls_data[p_id];
			MatModel::MaterialModel* pmat = mm_id_map.get_mm_by_id(pcl_data.mat_id);
			if (!pmat)
			{
				throw std::exception("func load_pcl_data_from_hdf5_file error: "
					"particle has no material model.");
			}
			pcl_data.to_pcl(md, 0, p_id, *pmat);
		}
		delete[] pcls_data;
	}

	rf.close_group(pcl_data_grp_id);
	return 0;
}


int output_material_model_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;
	hid_t mc_grp_id = rf.create_group(grp_id, "MaterialModel");
	output_material_model_container_to_hdf5_file(md, rf, mc_grp_id);
	rf.close_group(mc_grp_id);
	return 0;
}

int load_material_model_from_hdf5_file(
	Model_T3D_ME_mt &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;
	hid_t mc_grp_id = rf.open_group(grp_id, "MaterialModel");
	load_material_model_container_from_hdf5_file(md, rf, mc_grp_id);
	rf.close_group(mc_grp_id);
	return 0;
}

int output_rigid_cylinder_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!md.has_rigid_cylinder())
		return 0;

	hid_t rc_grp_id = rf.create_group(grp_id, "RigidCylinder");

	rf.write_attribute(rc_grp_id, "Kn_cont", md.Kn_cont);
	rf.write_attribute(rc_grp_id, "Kt_cont", md.Kt_cont);
	rf.write_attribute(rc_grp_id, "fric_ratio", md.fric_ratio);
	rf.write_attribute(rc_grp_id, "shear_strength", md.shear_strength);

	RigidObject_hdf5_utilities::output_rigid_cylinder_to_hdf5_file(md.rigid_cylinder, rf, rc_grp_id);

	rf.close_group(rc_grp_id);
	return 0;
}

int load_rigid_cylinder_from_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!rf.has_group(grp_id, "RigidCylinder"))
		return 0;

	hid_t rc_grp_id = rf.open_group(grp_id, "RigidCylinder");

	double Kn_cont, Kt_cont, fric_ratio, shear_strength;
	rf.read_attribute(rc_grp_id, "Kn_cont", Kn_cont);
	rf.read_attribute(rc_grp_id, "Kt_cont", Kt_cont);
	rf.read_attribute(rc_grp_id, "fric_ratio", fric_ratio);
	rf.read_attribute(rc_grp_id, "shear_strength", shear_strength);
	md.set_contact_param(Kn_cont, Kt_cont, fric_ratio, shear_strength);

	RigidObject_hdf5_utilities::load_rigid_cylinder_from_hdf5_file(md.rigid_cylinder, rf, rc_grp_id);

	md.rigid_cylinder_is_valid = true;
	rf.close_group(rc_grp_id);
	return 0;
}

int output_rigid_cone_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!md.has_rigid_cone())
		return 0;

	hid_t rc_grp_id = rf.create_group(grp_id, "RigidCone");

	rf.write_attribute(rc_grp_id, "Kn_cont", md.Kn_cont);
	rf.write_attribute(rc_grp_id, "Kt_cont", md.Kt_cont);
	rf.write_attribute(rc_grp_id, "fric_ratio", md.fric_ratio);
	rf.write_attribute(rc_grp_id, "shear_strength", md.shear_strength);

	RigidObject_hdf5_utilities::output_rigid_cone_to_hdf5_file(md.rigid_cone, rf, rc_grp_id);

	rf.close_group(rc_grp_id);
	return 0;
}

int load_rigid_cone_from_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!rf.has_group(grp_id, "RigidCone"))
		return 0;

	hid_t rc_grp_id = rf.open_group(grp_id, "RigidCone");

	double Kn_cont, Kt_cont, fric_ratio, shear_strength;
	rf.read_attribute(rc_grp_id, "Kn_cont", Kn_cont);
	rf.read_attribute(rc_grp_id, "Kt_cont", Kt_cont);
	rf.read_attribute(rc_grp_id, "fric_ratio", fric_ratio);
	rf.read_attribute(rc_grp_id, "shear_strength", shear_strength);
	md.set_contact_param(Kn_cont, Kt_cont, fric_ratio, shear_strength);

	RigidObject_hdf5_utilities::load_rigid_cone_from_hdf5_file(md.rigid_cone, rf, rc_grp_id);

	md.rigid_cone_is_valid = true;
	rf.close_group(rc_grp_id);
	return 0;
}

int output_rigid_cube_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!md.has_rigid_cube())
		return 0;

	hid_t rc_grp_id = rf.create_group(grp_id, "RigidCube");

	rf.write_attribute(rc_grp_id, "Kn_cont", md.Kn_cont);
	rf.write_attribute(rc_grp_id, "Kt_cont", md.Kt_cont);
	rf.write_attribute(rc_grp_id, "fric_ratio", md.fric_ratio);
	rf.write_attribute(rc_grp_id, "shear_strength", md.shear_strength);

	RigidObject_hdf5_utilities::output_rigid_cube_to_hdf5_file(md.rigid_cube, rf, rc_grp_id);

	rf.close_group(rc_grp_id);
	return 0;
}

int load_rigid_cube_from_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!rf.has_group(grp_id, "RigidCube"))
		return 0;

	hid_t rc_grp_id = rf.open_group(grp_id, "RigidCube");

	double Kn_cont, Kt_cont, fric_ratio, shear_strength;
	rf.read_attribute(rc_grp_id, "Kn_cont", Kn_cont);
	rf.read_attribute(rc_grp_id, "Kt_cont", Kt_cont);
	rf.read_attribute(rc_grp_id, "fric_ratio", fric_ratio);
	rf.read_attribute(rc_grp_id, "shear_strength", shear_strength);
	md.set_contact_param(Kn_cont, Kt_cont, fric_ratio, shear_strength);

	RigidObject_hdf5_utilities::load_rigid_cube_from_hdf5_file(md.rigid_cube, rf, rc_grp_id);

	md.rigid_cube_is_valid = true;
	rf.close_group(rc_grp_id);
	return 0;
}

int output_t3d_rigid_mesh_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!md.has_t3d_rigid_mesh())
		return 0;

	hid_t rm_grp_id = rf.create_group(grp_id, "RigidObjectByT3DMesh");

	rf.write_attribute(rm_grp_id, "Kn_cont", md.Kn_cont);
	rf.write_attribute(rm_grp_id, "Kt_cont", md.Kt_cont);
	rf.write_attribute(rm_grp_id, "fric_ratio", md.fric_ratio);
	rf.write_attribute(rm_grp_id, "shear_strength", md.shear_strength);

	RigidObject_hdf5_utilities::output_rigid_object_by_3dmesh_to_hdf5_file(
		md.get_t3d_rigid_mesh(), rf, rm_grp_id);

	rf.close_group(rm_grp_id);
	return 0;
}

int output_t3d_rigid_mesh_state_to_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!md.has_t3d_rigid_mesh())
		return 0;

	hid_t rm_grp_id = rf.create_group(grp_id, "RigidObjectByT3DMesh");
	RigidObject_hdf5_utilities::output_rigid_object_by_3dmesh_state_to_hdf5_file(
		md.get_t3d_rigid_mesh(), rf, rm_grp_id);
	rf.close_group(rm_grp_id);
	return 0;
}

int load_t3d_rigid_mesh_from_hdf5_file(
	Model_T3D_ME_mt& md,
	ResultFile_hdf5& rf,
	hid_t grp_id)
{
	if (grp_id < 0)
		return -1;

	if (!rf.has_group(grp_id, "RigidObjectByT3DMesh"))
		return 0;

	hid_t rm_grp_id = rf.open_group(grp_id, "RigidObjectByT3DMesh");

	double Kn_cont, Kt_cont, fric_ratio, shear_strength;
	rf.read_attribute(rm_grp_id, "Kn_cont", Kn_cont);
	rf.read_attribute(rm_grp_id, "Kt_cont", Kt_cont);
	rf.read_attribute(rm_grp_id, "fric_ratio", fric_ratio);
	rf.read_attribute(rm_grp_id, "shear_strength", shear_strength);
	md.set_contact_param(Kn_cont, Kt_cont, fric_ratio, shear_strength);

	RigidObject_hdf5_utilities::load_rigid_object_by_3dmesh_from_hdf5_file(
		md.get_t3d_rigid_mesh(), rf, rm_grp_id);

	md.rigid_t3d_mesh_is_valid = true;
	rf.close_group(rm_grp_id);
	return 0;
}

// output the whole model to ModelData
int output_model_to_hdf5_file(
	Model_T3D_ME_mt &md,
	ResultFile_hdf5 &rf)
{
	hid_t md_grp_id = rf.get_model_data_grp_id();
	// background mesh
	output_background_mesh_to_hdf5_file(md, rf, md_grp_id);
	// boundary condition
	output_boundary_condition_to_hdf5_file(md, rf, md_grp_id);
	// search mesh
	output_search_mesh_to_hdf5_file(md, rf, md_grp_id);

	// particle data
	output_ori_pcl_data_to_hdf5_file(md, rf, md_grp_id);
	// material model
	output_material_model_to_hdf5_file(md, rf, md_grp_id);

	// rigid object
	output_rigid_cylinder_to_hdf5_file(md, rf, md_grp_id);
	output_rigid_cone_to_hdf5_file(md, rf, md_grp_id);
	output_rigid_cube_to_hdf5_file(md, rf, md_grp_id);
	output_t3d_rigid_mesh_to_hdf5_file(md, rf, md_grp_id);
	return 0;
}

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(
	Model_T3D_ME_mt &md,
	Step_T3D_ME_mt &stp,
	ResultFile_hdf5 &rf,
	hid_t frame_grp_id)
{
	// particle data
	output_pcl_data_to_hdf5_file(md, stp, rf, frame_grp_id);
	// material model
	output_material_model_to_hdf5_file(md, rf, frame_grp_id);
	// rigid object
	output_rigid_cylinder_to_hdf5_file(md, rf, frame_grp_id);
	output_rigid_cone_to_hdf5_file(md, rf, frame_grp_id);
	output_rigid_cube_to_hdf5_file(md, rf, frame_grp_id);
	output_t3d_rigid_mesh_state_to_hdf5_file(md, rf, frame_grp_id);
	return 0;
}

int time_history_complete_output_to_hdf5_file(
	Model_T3D_ME_mt& md,
	Step_T3D_ME_TBB& stp,
	ResultFile_hdf5& rf,
	hid_t frame_grp_id)
{
	// particle data
	output_pcl_data_to_hdf5_file(md, stp, rf, frame_grp_id);
	// material model
	output_material_model_to_hdf5_file(md, rf, frame_grp_id);
	// rigid object
	output_rigid_cylinder_to_hdf5_file(md, rf, frame_grp_id);
	output_rigid_cone_to_hdf5_file(md, rf, frame_grp_id);
	output_rigid_cube_to_hdf5_file(md, rf, frame_grp_id);
	output_t3d_rigid_mesh_state_to_hdf5_file(md, rf, frame_grp_id);
	return 0;
}

int time_history_complete_output_to_hdf5_file(
	Model_T3D_ME_mt& md,
	Step_T3D_ME_mt_Geo& stp,
	ResultFile_hdf5& rf,
	hid_t frame_grp_id)
{
	// particle data
	output_pcl_data_to_hdf5_file(md, stp, rf, frame_grp_id);
	// material model
	output_material_model_to_hdf5_file(md, rf, frame_grp_id);
	// rigid object
	output_rigid_cylinder_to_hdf5_file(md, rf, frame_grp_id);
	output_rigid_cone_to_hdf5_file(md, rf, frame_grp_id);
	output_rigid_cube_to_hdf5_file(md, rf, frame_grp_id);
	output_t3d_rigid_mesh_state_to_hdf5_file(md, rf, frame_grp_id);
	return 0;
}

// load model data from hdf5
int load_me_mt_model_from_hdf5_file(
	Model_T3D_ME_mt &md,
	const char* hdf5_name)
{
	ResultFile_hdf5 rf;
	rf.open(hdf5_name);
	hid_t file_id = rf.get_file_id();
	if (file_id < 0)
		return -1;

	// model data
	hid_t md_grp_id = rf.get_model_data_grp_id();
	// background mesh
	load_background_mesh_from_hdf5_file(md, rf, md_grp_id);
	// search mesh
	load_search_mesh_from_hdf5_file(md, rf, md_grp_id);
	// boundary condition
	load_boundary_condition_from_hdf5_file(md, rf, md_grp_id);

	// material model
	load_material_model_from_hdf5_file(md, rf, md_grp_id);
	// particle data
	load_pcl_data_from_hdf5_file(md, rf, md_grp_id);

	// rigid object
	load_rigid_cylinder_from_hdf5_file(md, rf, md_grp_id);
	load_rigid_cone_from_hdf5_file(md, rf, md_grp_id);
	load_rigid_cube_from_hdf5_file(md, rf, md_grp_id);
	load_t3d_rigid_mesh_from_hdf5_file(md, rf, md_grp_id);

	return 0;
}

// load model data from hdf5
int load_me_mt_model_from_hdf5_file(
	Model_T3D_ME_mt& md,
	Step_T3D_ME_mt& step,
	const char* hdf5_name,
	const char* th_name,
	size_t frame_id)
{
	ResultFile_hdf5 rf;
	rf.open(hdf5_name);
	hid_t file_id = rf.get_file_id();
	if (file_id < 0)
		return -1;

	// model data
	hid_t md_grp_id = rf.get_model_data_grp_id();
	// background mesh
	load_background_mesh_from_hdf5_file(md, rf, md_grp_id);
	// search mesh
	load_search_mesh_from_hdf5_file(md, rf, md_grp_id);
	// boundary condition
	load_boundary_condition_from_hdf5_file(md, rf, md_grp_id);

	// time history
	hid_t th_grp_id = rf.get_time_history_grp_id();
	hid_t th_id = rf.open_group(th_grp_id, th_name);
	char th_frame_name[30];
	snprintf(th_frame_name, 30, "frame_%zu", frame_id);
	hid_t th_frame_id = rf.open_group(th_id, th_frame_name);

	// material model
	load_material_model_from_hdf5_file(md, rf, th_frame_id);
	// particle data
	load_pcl_data_from_hdf5_file(md, rf, th_frame_id);
	// rigid object
	load_rigid_cylinder_from_hdf5_file(md, rf, th_frame_id);
	load_rigid_cone_from_hdf5_file(md, rf, th_frame_id);
	load_rigid_cube_from_hdf5_file(md, rf, th_frame_id);;

	step.set_model(md);
	step.is_first_step = false;
	rf.read_attribute(th_frame_id, "total_time", step.start_time);
	rf.read_attribute(th_frame_id, "total_substep_num", step.prev_substep_num);

	rf.close_group(th_frame_id);
	rf.close_group(th_id);
	return 0;
}

// load model data from hdf5
int load_me_mt_model_from_hdf5_file(
	Model_T3D_ME_mt &md,
	Step_T3D_ME_TBB &step,
	const char* hdf5_name,
	const char* th_name,
	size_t frame_id)
{
	ResultFile_hdf5 rf;
	rf.open(hdf5_name);
	hid_t file_id = rf.get_file_id();
	if (file_id < 0)
		return -1;

	// model data
	hid_t md_grp_id = rf.get_model_data_grp_id();
	// background mesh
	load_background_mesh_from_hdf5_file(md, rf, md_grp_id);
	// search mesh
	load_search_mesh_from_hdf5_file(md, rf, md_grp_id);
	// boundary condition
	load_boundary_condition_from_hdf5_file(md, rf, md_grp_id);

	// time history
	hid_t th_grp_id = rf.get_time_history_grp_id();
	hid_t th_id = rf.open_group(th_grp_id, th_name);
	char th_frame_name[30];
	snprintf(th_frame_name, 30, "frame_%zu", frame_id);
	hid_t th_frame_id = rf.open_group(th_id, th_frame_name);

	// material model
	load_material_model_from_hdf5_file(md, rf, th_frame_id);
	// particle data
	load_pcl_data_from_hdf5_file(md, rf, th_frame_id);
	// rigid object
	load_rigid_cylinder_from_hdf5_file(md, rf, th_frame_id);
	load_rigid_cone_from_hdf5_file(md, rf, th_frame_id);
	load_rigid_cube_from_hdf5_file(md, rf, th_frame_id);;

	step.set_model(md);
	step.is_first_step = false;
	rf.read_attribute(th_frame_id, "total_time", step.start_time);
	rf.read_attribute(th_frame_id, "total_substep_num", step.prev_substep_num);

	rf.close_group(th_frame_id);
	rf.close_group(th_id);
	return 0;
}

};
