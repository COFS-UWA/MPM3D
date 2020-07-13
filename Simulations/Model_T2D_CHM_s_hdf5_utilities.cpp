#include "Simulations_pcp.h"

#include "Model_hdf5_utilities.h"
#include "Model_T2D_CHM_s_hdf5_utilities.h"

namespace Model_T2D_CHM_s_hdf5_utilities
{
using namespace Model_hdf5_utilities;

int output_background_mesh_to_hdf5_file(
	Model_T2D_CHM_s &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bg_mesh_grp_id = rf.create_group(grp_id, "BackgroundMesh");

	size_t node_num = md.get_node_num();
	Model_T2D_CHM_s::Node* nodes = md.get_nodes();
	size_t elem_num = md.get_elem_num();
	Model_T2D_CHM_s::Element* elems = md.get_elems();

	// bg mesh properties
	const char* bg_mesh_type = "T2D";
	rf.write_attribute(bg_mesh_grp_id, "type", strlen(bg_mesh_type), bg_mesh_type);
	rf.write_attribute(bg_mesh_grp_id, "node_num", node_num);
	rf.write_attribute(bg_mesh_grp_id, "element_num", elem_num);

	// fluid properties
	rf.write_attribute(bg_mesh_grp_id, "Kf", md.Kf);
	rf.write_attribute(bg_mesh_grp_id, "k", md.k);
	rf.write_attribute(bg_mesh_grp_id, "miu", md.miu);

	// bg grid properties
	rf.write_attribute(bg_mesh_grp_id, "bg_grid_hx", md.get_bg_grid_hx());
	rf.write_attribute(bg_mesh_grp_id, "bg_grid_hy", md.get_bg_grid_hy());
	
	// node coordinates
	Node2DData *nodes_data = new Node2DData[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node2DData &node_data = nodes_data[n_id];
		Model_T2D_CHM_s::Node n = nodes[n_id];
		node_data.id = n.id;
		node_data.x = n.x;
		node_data.y = n.y;
	}
	hid_t nd_dt_id = get_nd_2d_dt_id();
	rf.write_dataset(
		bg_mesh_grp_id,
		"NodeCoordinate",
		node_num,
		nodes_data,
		nd_dt_id
	);
	H5Tclose(nd_dt_id);
	delete[] nodes_data;
	// element indices
	Elem2DData *elems_data = new Elem2DData[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Elem2DData &elem_data = elems_data[e_id];
		Model_T2D_CHM_s::Element &e = elems[e_id];
		elem_data.id = e.id;
		elem_data.n1 = e.n1;
		elem_data.n2 = e.n2;
		elem_data.n3 = e.n3;
	}
	hid_t ed_dt_id = get_ed_2d_dt_id();
	rf.write_dataset(
		bg_mesh_grp_id,
		"ElementTopology",
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
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bg_mesh_grp_id = rf.open_group(grp_id, "BackgroundMesh");
	
	// liquid properties
	rf.read_attribute(bg_mesh_grp_id, "Kf", md.Kf);
	rf.read_attribute(bg_mesh_grp_id, "k", md.k);
	rf.read_attribute(bg_mesh_grp_id, "miu", md.miu);

	// mesh properties
	size_t elem_num, node_num;
	rf.read_attribute(bg_mesh_grp_id, "node_num", node_num);
	rf.read_attribute(bg_mesh_grp_id, "element_num", elem_num);

	// nodes
	Node2DData *nodes_data = new Node2DData[node_num];
	hid_t nd_dt_id = get_nd_2d_dt_id();
	rf.read_dataset(
		bg_mesh_grp_id,
		"NodeCoordinate",
		node_num,
		nodes_data,
		nd_dt_id
	);
	H5Tclose(nd_dt_id);
	md.alloc_nodes(node_num);
	Model_T2D_CHM_s::Node* nodes = md.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node2DData &node_data = nodes_data[n_id];
		Model_T2D_CHM_s::Node &n = nodes[n_id];
		n.id = node_data.id;
		n.x = node_data.x;
		n.y = node_data.y;
	}
	delete[] nodes_data;

	// elements
	Elem2DData *elems_data = new Elem2DData[elem_num];
	hid_t ed_dt_id = get_ed_2d_dt_id();
	rf.read_dataset(
		bg_mesh_grp_id,
		"ElementTopology",
		elem_num,
		elems_data,
		ed_dt_id
	);
	H5Tclose(ed_dt_id);
	md.alloc_elements(elem_num);
	Model_T2D_CHM_s::Element* elems = md.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Elem2DData &elem_data = elems_data[e_id];
		Model_T2D_CHM_s::Element &elem = elems[e_id];
		elem.id = elem_data.id;
		elem.n1 = elem_data.n1;
		elem.n2 = elem_data.n2;
		elem.n3 = elem_data.n3;
	}
	delete[] elems_data;

	md.init_mesh_properties_after_loading();

	// init bg_grid
	double bg_grid_hx, bg_grid_hy;
	rf.read_attribute(bg_mesh_grp_id, "bg_grid_hx", bg_grid_hx);
	rf.read_attribute(bg_mesh_grp_id, "bg_grid_hy", bg_grid_hy);

	md.init_search_grid(bg_grid_hx, bg_grid_hy);

	rf.close_group(bg_mesh_grp_id);
	return 0;
}

int output_boundary_condition_to_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	hid_t md_id = rf.get_model_data_grp_id();

	hid_t bc_id = rf.create_group(md_id, "BoundaryCondition");
	rf.write_attribute(bc_id, "asx_num", md.asx_num);
	rf.write_attribute(bc_id, "asy_num", md.asy_num);
	rf.write_attribute(bc_id, "afx_num", md.afx_num);
	rf.write_attribute(bc_id, "afy_num", md.afy_num);
	rf.write_attribute(bc_id, "vsx_num", md.vsx_num);
	rf.write_attribute(bc_id, "vsy_num", md.vsy_num);
	rf.write_attribute(bc_id, "vfx_num", md.vfx_num);
	rf.write_attribute(bc_id, "vfy_num", md.vfy_num);
	rf.write_attribute(bc_id, "tx_num", md.tx_num);
	rf.write_attribute(bc_id, "ty_num", md.ty_num);
	rf.write_attribute(bc_id, "bfx_num", md.bfx_num);
	rf.write_attribute(bc_id, "bfy_num", md.bfy_num);
	
	// Acceleration bcs
	hid_t abc_dt_id = get_abc_dt_id();
	AccelerationBCData* abcds;
	if (md.asx_num)
	{
		abcds = new AccelerationBCData[md.asx_num];
		for (size_t a_id = 0; a_id < md.asx_num; ++a_id)
		{
			AccelerationBC& abc = md.asxs[a_id];
			AccelerationBCData& abcd = abcds[a_id];
			abcd.from_abc(abc);
		}
		rf.write_dataset(bc_id, "AccelerationBC_sx", md.asx_num, abcds, abc_dt_id);
		delete[] abcds;
	}
	if (md.asy_num)
	{
		abcds = new AccelerationBCData[md.asy_num];
		for (size_t a_id = 0; a_id < md.asy_num; ++a_id)
		{
			AccelerationBC &abc = md.asys[a_id];
			AccelerationBCData& abcd = abcds[a_id];
			abcd.from_abc(abc);
		}
		rf.write_dataset(bc_id, "AccelerationBC_sy", md.asx_num, abcds, abc_dt_id);
		delete[] abcds;
	}
	if (md.afx_num)
	{
		abcds = new AccelerationBCData[md.afx_num];
		for (size_t a_id = 0; a_id < md.afx_num; ++a_id)
		{
			AccelerationBC& abc = md.afxs[a_id];
			AccelerationBCData& abcd = abcds[a_id];
			abcd.from_abc(abc);
		}
		rf.write_dataset(bc_id, "AccelerationBC_fx", md.afx_num, abcds, abc_dt_id);
		delete[] abcds;
	}
	if (md.afy_num)
	{
		abcds = new AccelerationBCData[md.afy_num];
		for (size_t a_id = 0; a_id < md.afy_num; ++a_id)
		{
			AccelerationBC& abc = md.afys[a_id];
			AccelerationBCData& abcd = abcds[a_id];
			abcd.from_abc(abc);
		}
		rf.write_dataset(bc_id, "AccelerationBC_fy", md.afy_num, abcds, abc_dt_id);
		delete[] abcds;
	}
	H5Tclose(abc_dt_id);

	// Velocity bcs
	hid_t vbc_dt_id = get_vbc_dt_id();
	VelocityBCData *vbcds;
	if (md.vsx_num)
	{
		vbcds = new VelocityBCData[md.vsx_num];
		for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
		{
			VelocityBC& vbc = md.vsxs[v_id];
			VelocityBCData& vbcd = vbcds[v_id];
			vbcd.from_vbc(vbc);
		}
		rf.write_dataset(bc_id, "VelocityBC_sx", md.vsx_num, vbcds, vbc_dt_id);
		delete[] vbcds;
	}
	if (md.vsy_num)
	{
		vbcds = new VelocityBCData[md.vsy_num];
		for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
		{
			VelocityBC& vbc = md.vsys[v_id];
			VelocityBCData& vbcd = vbcds[v_id];
			vbcd.from_vbc(vbc);
		}
		rf.write_dataset(bc_id, "VelocityBC_sy", md.vsy_num, vbcds, vbc_dt_id);
		delete[] vbcds;
	}
	if (md.vfx_num)
	{
		vbcds = new VelocityBCData[md.vfx_num];
		for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
		{
			VelocityBC& vbc = md.vfxs[v_id];
			VelocityBCData& vbcd = vbcds[v_id];
			vbcd.from_vbc(vbc);
		}
		rf.write_dataset(bc_id, "VelocityBC_fx", md.vfx_num, vbcds, vbc_dt_id);
		delete[] vbcds;
	}
	if (md.vfy_num)
	{
		vbcds = new VelocityBCData[md.vfy_num];
		for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
		{
			VelocityBC& vbc = md.vfys[v_id];
			VelocityBCData& vbcd = vbcds[v_id];
			vbcd.from_vbc(vbc);
		}
		rf.write_dataset(bc_id, "VelocityBC_fy", md.vfy_num, vbcds, vbc_dt_id);
		delete[] vbcds;
	}
	H5Tclose(vbc_dt_id);

	// traction bc
	hid_t tbc_dt_id = get_tbc_pcl_dt_id();
	TractionBCAtPclData* tbcds;
	if (md.tx_num)
	{
		tbcds = new TractionBCAtPclData[md.tx_num];
		for (size_t t_id = 0; t_id < md.tx_num; ++t_id)
		{
			TractionBCAtPcl& tbc = md.txs[t_id];
			TractionBCAtPclData& tbcd = tbcds[t_id];
			tbcd.from_tbc(tbc);
		}
		rf.write_dataset(bc_id, "TractionBCAtPcl_x", md.tx_num, tbcds, tbc_dt_id);
		delete[] tbcds;
	}
	if (md.ty_num)
	{
		tbcds = new TractionBCAtPclData[md.ty_num];
		for (size_t t_id = 0; t_id < md.ty_num; ++t_id)
		{
			TractionBCAtPcl& tbc = md.tys[t_id];
			TractionBCAtPclData& tbcd = tbcds[t_id];
			tbcd.from_tbc(tbc);
		}
		rf.write_dataset(bc_id, "TractionBCAtPcl_y", md.ty_num, tbcds, tbc_dt_id);
		delete[] tbcds;
	}
	H5Tclose(tbc_dt_id);

	// body force bc
	hid_t bf_dt_id = get_bf_pcl_dt_id();
	BodyForceAtPclData *bfds;
	if (md.bfx_num)
	{
		bfds = new BodyForceAtPclData[md.bfx_num];
		for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
		{
			BodyForceAtPcl& bf = md.bfxs[bf_id];
			BodyForceAtPclData& bfd = bfds[bf_id];
			bfd.from_bf(bf);
		}
		rf.write_dataset(bc_id, "BodyForceAtPcl_x", md.bfx_num, bfds, bf_dt_id);
		delete[] bfds;
	}
	if (md.bfy_num)
	{
		bfds = new BodyForceAtPclData[md.bfy_num];
		for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
		{
			BodyForceAtPcl& bf = md.bfys[bf_id];
			BodyForceAtPclData& bfd = bfds[bf_id];
			bfd.from_bf(bf);
		}
		rf.write_dataset(bc_id, "BodyForceAtPcl_y", md.bfy_num, bfds, bf_dt_id);
		delete[] bfds;
	}
	H5Tclose(bf_dt_id);

	rf.close_group(bc_id);
	return 0;
}

int load_boundary_condition_from_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bc_id = rf.open_group(grp_id, "BoundaryCondition");
	rf.read_attribute(bc_id, "asx_num", md.asx_num);
	rf.read_attribute(bc_id, "asy_num", md.asy_num);
	rf.read_attribute(bc_id, "afx_num", md.afx_num);
	rf.read_attribute(bc_id, "afy_num", md.afy_num);
	rf.read_attribute(bc_id, "vsx_num", md.vsx_num);
	rf.read_attribute(bc_id, "vsy_num", md.vsy_num);
	rf.read_attribute(bc_id, "vfx_num", md.vfx_num);
	rf.read_attribute(bc_id, "vfy_num", md.vfy_num);
	rf.read_attribute(bc_id, "tx_num", md.tx_num);
	rf.read_attribute(bc_id, "ty_num", md.ty_num);
	rf.read_attribute(bc_id, "bfx_num", md.bfx_num);
	rf.read_attribute(bc_id, "bfy_num", md.bfy_num);

	// acceleration bcs
	hid_t abc_dt_id = get_abc_dt_id();
	AccelerationBCData *abcds;
	if (md.asx_num)
	{
		md.init_asxs(md.asx_num);
		abcds = new AccelerationBCData[md.asx_num];
		rf.read_dataset(bc_id, "AccelerationBC_sx", md.asx_num, abcds, abc_dt_id);
		for (size_t a_id = 0; a_id < md.asx_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& abc = md.asxs[a_id];
			abcd.to_abc(abc);
		}
		delete[] abcds;
	}
	if (md.asy_num)
	{
		md.init_asys(md.asy_num);
		abcds = new AccelerationBCData[md.asy_num];
		rf.read_dataset(bc_id, "AccelerationBC_sy", md.asy_num, abcds, abc_dt_id);
		for (size_t a_id = 0; a_id < md.asy_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& abc = md.asys[a_id];
			abcd.to_abc(abc);
		}
		delete[] abcds;
	}
	if (md.afx_num)
	{
		md.init_afxs(md.afx_num);
		abcds = new AccelerationBCData[md.afx_num];
		rf.read_dataset(bc_id, "AccelerationBC_fx", md.afx_num, abcds, abc_dt_id);
		for (size_t a_id = 0; a_id < md.afx_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& abc = md.afxs[a_id];
			abcd.to_abc(abc);
		}
		delete[] abcds;
	}
	if (md.afy_num)
	{
		md.init_afys(md.afy_num);
		abcds = new AccelerationBCData[md.afy_num];
		rf.read_dataset(bc_id, "AccelerationBC_fy", md.afy_num, abcds, abc_dt_id);
		for (size_t a_id = 0; a_id < md.afy_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& abc = md.afys[a_id];
			abcd.to_abc(abc);
		}
		delete[] abcds;
	}
	H5Tclose(abc_dt_id);

	// velocity bcs
	hid_t vbc_dt_id = get_vbc_dt_id();
	VelocityBCData* vbcds;
	if (md.vsx_num)
	{
		md.init_vsxs(md.vsx_num);
		vbcds = new VelocityBCData[md.vsx_num];
		rf.read_dataset(bc_id, "VelocityBC_sx", md.vsx_num, vbcds, vbc_dt_id);
		for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
		{
			VelocityBCData& vbcd = vbcds[v_id];
			VelocityBC &vbc = md.vsxs[v_id];
			vbcd.to_vbc(vbc);
		}
		delete[] vbcds;
	}
	if (md.vsy_num)
	{
		md.init_vsys(md.vsy_num);
		vbcds = new VelocityBCData[md.vsy_num];
		rf.read_dataset(bc_id, "VelocityBC_sy", md.vsy_num, vbcds, vbc_dt_id);
		for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
		{
			VelocityBCData &vbcd = vbcds[v_id];
			VelocityBC &vbc = md.vsys[v_id];
			vbcd.to_vbc(vbc);
		}
		delete[] vbcds;
	}
	if (md.vfx_num)
	{
		md.init_vfxs(md.vfx_num);
		vbcds = new VelocityBCData[md.vfx_num];
		rf.read_dataset(bc_id, "VelocityBC_fx", md.vfx_num, vbcds, vbc_dt_id);
		for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
		{
			VelocityBCData &vbcd = vbcds[v_id];
			VelocityBC &vbc = md.vfxs[v_id];
			vbcd.to_vbc(vbc);
		}
		delete[] vbcds;
	}
	if (md.vfy_num)
	{
		md.init_vfys(md.vfy_num);
		vbcds = new VelocityBCData[md.vfy_num];
		rf.read_dataset(bc_id, "VelocityBC_fy", md.vfy_num, vbcds, vbc_dt_id);
		for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
		{
			VelocityBCData& vbcd = vbcds[v_id];
			VelocityBC& vbc = md.vfys[v_id];
			vbcd.to_vbc(vbc);
		}
		delete[] vbcds;
	}
	H5Tclose(vbc_dt_id);

	// traction bc
	hid_t tbc_dt_id = get_tbc_pcl_dt_id();
	TractionBCAtPclData *tbcds;
	if (md.tx_num)
	{
		md.init_txs(md.tx_num);
		tbcds = new TractionBCAtPclData[md.tx_num];
		rf.read_dataset(bc_id, "TractionBCAtPcl_x", md.tx_num, tbcds, tbc_dt_id);
		for (size_t t_id = 0; t_id < md.tx_num; ++t_id)
		{
			TractionBCAtPclData& tbcd = tbcds[t_id];
			TractionBCAtPcl& tbc = md.txs[t_id];
			tbcd.to_tbc(tbc);
		}
		delete[] tbcds;
	}
	if (md.ty_num)
	{
		md.init_tys(md.ty_num);
		tbcds = new TractionBCAtPclData[md.ty_num];
		rf.read_dataset(bc_id, "TractionBCAtPcl_y", md.ty_num, tbcds, tbc_dt_id);
		for (size_t t_id = 0; t_id < md.ty_num; ++t_id)
		{
			TractionBCAtPclData &tbcd = tbcds[t_id];
			TractionBCAtPcl &tbc = md.tys[t_id];
			tbcd.to_tbc(tbc);
		}
		delete[] tbcds;
	}
	H5Tclose(tbc_dt_id);

	// body force
	hid_t bf_dt_id = get_bf_pcl_dt_id();
	BodyForceAtPclData* bfds;
	if (md.bfx_num)
	{
		md.init_bfxs(md.bfx_num);
		bfds = new BodyForceAtPclData[md.bfx_num];
		rf.read_dataset(bc_id, "BodyForceAtPcl_x", md.bfx_num, bfds, bf_dt_id);
		for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
		{
			BodyForceAtPclData& bfd = bfds[bf_id];
			BodyForceAtPcl& bf = md.bfxs[bf_id];
			bfd.to_bf(bf);
		}
		delete[] bfds;
	}
	if (md.bfy_num)
	{
		md.init_bfys(md.bfy_num);
		bfds = new BodyForceAtPclData[md.bfy_num];
		rf.read_dataset(bc_id, "BodyForceAtPcl_y", md.bfy_num, bfds, bf_dt_id);
		for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
		{
			BodyForceAtPclData& bfd = bfds[bf_id];
			BodyForceAtPcl& bf = md.bfys[bf_id];
			bfd.to_bf(bf);
		}
		delete[] bfds;
	}
	H5Tclose(bf_dt_id);

	rf.close_group(bc_id);
	return 0;
}


int output_pcl_data_to_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.create_group(grp_id, "ParticleData");

	size_t pcl_num = md.get_pcl_num();
	Model_T2D_CHM_s::Particle* pcls = md.get_pcls();
	rf.write_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

	ParticleData* pcl_data = new ParticleData[md.pcl_num];
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		ParticleData& pd = pcl_data[p_id];
		Model_T2D_CHM_s::Particle& pcl = md.pcls[p_id];
		pd.from_pcl(pcl);
	}
	hid_t pcl_dt_id = get_pcl_dt_id();
	int res = rf.write_dataset(
		pcl_data_grp_id,
		"field",
		pcl_num,
		pcl_data,
		pcl_dt_id
	);
	H5Tclose(pcl_dt_id);
	delete[] pcl_data;

	rf.close_group(pcl_data_grp_id);
	return res;
}

int load_pcl_data_from_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.open_group(grp_id, "ParticleData");

	size_t pcl_num;
	rf.read_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

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
		return res;
	}

	md.alloc_pcls(pcl_num);
	Model_T2D_CHM_s::Particle* pcls = md.get_pcls();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		ParticleData& pcl_data = pcls_data[p_id];
		Model_T2D_CHM_s::Particle& pcl = md.pcls[p_id];
		pcl_data.to_pcl(pcl);
	}
	delete[] pcls_data;

	rf.close_group(pcl_data_grp_id);
	return 0;
}

namespace
{
	template <typename Particle>
	inline size_t get_pcl_id(void *ext_data)
	{ return static_cast<Particle *>(ext_data)->id; }
};

int output_material_model_to_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	size_t mm_id, mm_num;
	hid_t mm_grp_id = rf.create_group(grp_id, "MaterialModel");

	// linear elasticity
	mm_num = md.get_num_LinearElasticity();
	if (mm_num)
	{
		rf.write_attribute(mm_grp_id, "LinearElasticity_num", mm_num);

		LinearElasticityStateData* mm_data = new LinearElasticityStateData[mm_num];
		mm_id = 0;
		for (MatModel::LinearElasticity* iter = md.first_LinearElasticity();
			md.is_not_end_LinearElasticity(iter);
			iter = md.next_LinearElasticity(iter))
		{
			mm_data[mm_id].pcl_id = get_pcl_id<Model_T2D_CHM_s::Particle>(iter->ext_data);
			mm_data[mm_id].from_mm(*iter);
			++mm_id;
		}
		hid_t le_dt_id = get_le_hdf5_dt_id();
		rf.write_dataset(mm_grp_id, "LinearElasticity", mm_num, mm_data, le_dt_id);
		H5Tclose(le_dt_id);
		delete[] mm_data;
	}

	mm_num = md.get_num_ModifiedCamClay();
	if (mm_num)
	{
		rf.write_attribute(mm_grp_id, "ModifiedCamClay_num", mm_num);

		ModifiedCamClayStateData* mm_data = new ModifiedCamClayStateData[mm_num];
		mm_id = 0;
		for (MatModel::ModifiedCamClay* iter = md.first_ModifiedCamClay();
			md.is_not_end_ModifiedCamClay(iter);
			iter = md.next_ModifiedCamClay(iter))
		{
			mm_data[mm_id].pcl_id = get_pcl_id<Model_T2D_CHM_s::Particle>(iter->ext_data);
			mm_data[mm_id].from_mm(*iter);
			++mm_id;
		}
		hid_t mcc_dt_id = get_mcc_hdf5_dt_id();
		rf.write_dataset(mm_grp_id, "ModifiedCamClay", mm_num, mm_data, mcc_dt_id);
		H5Tclose(mcc_dt_id);
		delete[] mm_data;
	}

	rf.close_group(mm_grp_id);
	return 0;
}

int load_material_model_from_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t mm_dset_id;
	size_t mm_num;
	hid_t mm_grp_id = rf.open_group(grp_id, "MaterialModel");

	// linear elasticity
	if (rf.has_dataset(mm_grp_id, "LinearElasticity"))
	{
		rf.read_attribute(mm_grp_id, "LinearElasticity_num", mm_num);

		// get data
		LinearElasticityStateData* mm_data = new LinearElasticityStateData[mm_num];
		hid_t le_dt_id = get_le_hdf5_dt_id();
		rf.read_dataset(
			mm_grp_id,
			"LinearElasticity",
			mm_num,
			mm_data,
			le_dt_id
		);
		H5Tclose(le_dt_id);
		MatModel::LinearElasticity* mms = md.add_LinearElasticity(mm_num);
		for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
		{
			LinearElasticityStateData& mmd = mm_data[mm_id];
			MatModel::LinearElasticity& mm = mms[mm_id];
			mmd.to_mm(mm);
			md.pcls[mmd.pcl_id].set_mat_model(mm);
		}
		delete[] mm_data;
	}

	// modified cam clay
	if (rf.has_dataset(mm_grp_id, "ModifiedCamClay"))
	{
		rf.read_attribute(mm_grp_id, "ModifiedCamClay_num", mm_num);

		// get data
		ModifiedCamClayStateData* mm_data = new ModifiedCamClayStateData[mm_num];
		hid_t mcc_dt_id = get_mcc_hdf5_dt_id();
		rf.read_dataset(
			mm_grp_id,
			"ModifiedCamClay",
			mm_num,
			mm_data,
			mcc_dt_id
		);
		H5Tclose(mcc_dt_id);
		MatModel::ModifiedCamClay* mms = md.add_ModifiedCamClay(mm_num);
		for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
		{
			ModifiedCamClayStateData& mmd = mm_data[mm_id];
			MatModel::ModifiedCamClay& mm = mms[mm_id];
			mmd.to_mm(mm);
			md.pcls[mmd.pcl_id].set_mat_model(mm);
		}
		delete[] mm_data;
	}

	rf.close_group(mm_grp_id);
	return 0;
}

int output_rigid_circle_to_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	if (!md.rigid_circle_is_valid())
		return 0;
	
	hid_t rb_grp_id = rf.create_group(grp_id, "RigidCircle");

	rf.write_attribute(rb_grp_id, "Ks_cont", md.Ks_cont);
	rf.write_attribute(rb_grp_id, "Kf_cont", md.Kf_cont);

	RigidCircle& rc = md.get_rigid_circle();
	rf.write_attribute(rb_grp_id, "radius", rc.get_radius());
	rf.write_attribute(rb_grp_id, "density", rc.get_density());
	rf.write_attribute(rb_grp_id, "rfx", rc.get_rfx());
	rf.write_attribute(rb_grp_id, "rfy", rc.get_rfy());
	rf.write_attribute(rb_grp_id, "rm", rc.get_rm());
	rf.write_attribute(rb_grp_id, "ax", rc.get_ax());
	rf.write_attribute(rb_grp_id, "ay", rc.get_ay());
	rf.write_attribute(rb_grp_id, "a_angle", rc.get_a_ang());
	rf.write_attribute(rb_grp_id, "vx", rc.get_vx());
	rf.write_attribute(rb_grp_id, "vy", rc.get_vy());
	rf.write_attribute(rb_grp_id, "v_angle", rc.get_v_ang());
	rf.write_attribute(rb_grp_id, "x", rc.get_x());
	rf.write_attribute(rb_grp_id, "y", rc.get_y());
	rf.write_attribute(rb_grp_id, "angle", rc.get_ang());
	
	// boundary conditions
	rf.write_attribute(rb_grp_id, "rfx_bc", rc.get_rfx_bc());
	rf.write_attribute(rb_grp_id, "rfy_bc", rc.get_rfy_bc());
	rf.write_attribute(rb_grp_id, "rm_bc", rc.get_rm_bc());
	if (rc.has_ax_bc())
		rf.write_attribute(rb_grp_id, "ax_bc", rc.get_ax_bc());
	if (rc.has_ay_bc())
		rf.write_attribute(rb_grp_id, "ay_bc", rc.get_ay_bc());
	if (rc.has_a_ang_bc())
		rf.write_attribute(rb_grp_id, "a_ang_bc", rc.get_a_ang_bc());
	if (rc.has_vx_bc())
		rf.write_attribute(rb_grp_id, "vx_bc", rc.get_vx_bc());
	if (rc.has_vy_bc())
		rf.write_attribute(rb_grp_id, "vy_bc", rc.get_vy_bc());
	if (rc.has_v_ang_bc())
		rf.write_attribute(rb_grp_id, "v_ang_bc", rc.get_v_ang());

	rf.close_group(rb_grp_id);
	return 0;
}

int load_rigid_circle_from_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	if (!rf.has_group(grp_id, "RigidCircle"))
		return 0;
	
	hid_t rb_grp_id = rf.open_group(grp_id, "RigidCircle");

	double Ks_cont, Kf_cont;
	rf.read_attribute(rb_grp_id, "Ks_cont", Ks_cont);
	rf.read_attribute(rb_grp_id, "Kf_cont", Kf_cont);

	RigidCircle& rc = md.get_rigid_circle();
	double rc_radius, rc_density;
	double rc_ax, rc_ay, rc_a_ang;
	double rc_vx, rc_vy, rc_v_ang;
	double rc_x, rc_y, rc_ang;
	rf.read_attribute(rb_grp_id, "radius", rc_radius);
	rf.read_attribute(rb_grp_id, "density", rc_density);
	rf.read_attribute(rb_grp_id, "ax", rc_ax);
	rf.read_attribute(rb_grp_id, "ay", rc_ay);
	rf.read_attribute(rb_grp_id, "a_angle", rc_a_ang);
	rf.read_attribute(rb_grp_id, "vx", rc_vx);
	rf.read_attribute(rb_grp_id, "vy", rc_vy);
	rf.read_attribute(rb_grp_id, "v_angle", rc_v_ang);
	rf.read_attribute(rb_grp_id, "x", rc_x);
	rf.read_attribute(rb_grp_id, "y", rc_y);
	rf.read_attribute(rb_grp_id, "angle", rc_ang);
	md.set_rigid_circle_state(
		Ks_cont, Kf_cont,
		rc_x, rc_y, rc_ang,
		rc_radius, rc_density,
		rc_ax, rc_ay, rc_a_ang,
		rc_vx, rc_vy, rc_v_ang
		);

	rf.close_group(rb_grp_id);
	return 0;
}

// output the whole model to ModelData
int output_model_to_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf
	)
{
	hid_t md_grp_id = rf.get_model_data_grp_id();
	// background mesh
	output_background_mesh_to_hdf5_file(md, rf, md_grp_id);
	// boundary condition
	output_boundary_condition_to_hdf5_file(md, rf, md_grp_id);
	// particle data
	output_pcl_data_to_hdf5_file(md, rf, md_grp_id);
	// material model
	output_material_model_to_hdf5_file(md, rf, md_grp_id);
	// rigid object
	output_rigid_circle_to_hdf5_file(md, rf, md_grp_id);
	return 0;
}

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(
	Model_T2D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t frame_grp_id
	)
{
	// particle data
	output_pcl_data_to_hdf5_file(md, rf, frame_grp_id);
	// material model
	output_material_model_to_hdf5_file(md, rf, frame_grp_id);
	// rigid object
	output_rigid_circle_to_hdf5_file(md, rf, frame_grp_id);
	return 0;
}

// load model data from hdf5 to model data
int load_CHM_s_model_from_hdf5_file(
	Model_T2D_CHM_s& md,
	const char* hdf5_name,
	const char* th_name,
	size_t frame_id
	)
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
	// boundary condition
	load_boundary_condition_from_hdf5_file(md, rf, md_grp_id);

	// time history
	hid_t th_grp_id = rf.get_time_history_grp_id();
	hid_t th_id = rf.open_group(th_grp_id, th_name);
	char th_frame_name[30];
	snprintf(th_frame_name, 30, "frame_%zu", frame_id);
	hid_t th_frame_id = rf.open_group(th_id, th_frame_name);
	// particle data
	load_pcl_data_from_hdf5_file(md, rf, th_frame_id);
	// material model
	load_material_model_from_hdf5_file(md, rf, th_frame_id);
	// rigid object
	load_rigid_circle_from_hdf5_file(md, rf, th_frame_id);
	rf.close_group(th_frame_id);
	rf.close_group(th_id);

	return 0;
}

};