#include "Simulations_pcp.h"

#include <exception>

#include "MatModelIdToPointerMap.h"
#include "Model_hdf5_utilities.h"
#include "RigidBody/RigidBody_hdf5_utilities.h"
#include "Model_T2D_ME_s_hdf5_utilities.h"

namespace Model_T2D_ME_s_hdf5_utilities
{
using namespace Model_hdf5_utilities;

int output_background_mesh_to_hdf5_file(
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bg_mesh_grp_id = rf.create_group(grp_id, "BackgroundMesh");

	size_t node_num = md.get_node_num();
	Model_T2D_ME_s::Node* nodes = md.get_nodes();
	size_t elem_num = md.get_elem_num();
	Model_T2D_ME_s::Element* elems = md.get_elems();

	// bg mesh properties
	const char* bg_mesh_type = "T2D";
	rf.write_attribute(bg_mesh_grp_id, "type", strlen(bg_mesh_type), bg_mesh_type);
	rf.write_attribute(bg_mesh_grp_id, "node_num", node_num);
	rf.write_attribute(bg_mesh_grp_id, "element_num", elem_num);

	// bg grid properties
	rf.write_attribute(bg_mesh_grp_id, "bg_grid_hx", md.get_bg_grid_hx());
	rf.write_attribute(bg_mesh_grp_id, "bg_grid_hy", md.get_bg_grid_hy());

	// node coordinates
	Node2DData *nodes_data = new Node2DData[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node2DData &node_data = nodes_data[n_id];
		Model_T2D_ME_s::Node n = nodes[n_id];
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
		Model_T2D_ME_s::Element &e = elems[e_id];
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
	Model_T2D_ME_s &md,
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
	Model_T2D_ME_s::Node* nodes = md.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node2DData &node_data = nodes_data[n_id];
		Model_T2D_ME_s::Node &n = nodes[n_id];
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
	Model_T2D_ME_s::Element* elems = md.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Elem2DData &elem_data = elems_data[e_id];
		Model_T2D_ME_s::Element &elem = elems[e_id];
		elem.id = elem_data.id;
		elem.n1 = elem_data.n1;
		elem.n2 = elem_data.n2;
		elem.n3 = elem_data.n3;
	}
	delete[] elems_data;

	md.init_mesh_properties_after_loading();
	md.init_mesh_shape_funcs();

	// init bg_grid
	double bg_grid_hx, bg_grid_hy;
	rf.read_attribute(bg_mesh_grp_id, "bg_grid_hx", bg_grid_hx);
	rf.read_attribute(bg_mesh_grp_id, "bg_grid_hy", bg_grid_hy);
	
	md.init_search_grid(bg_grid_hx, bg_grid_hy);

	rf.close_group(bg_mesh_grp_id);
	return 0;
}

int output_boundary_condition_to_hdf5_file(
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bc_grp_id = rf.create_group(grp_id, "BoundaryCondition");

	rf.write_attribute(bc_grp_id, "bfx_num", md.bfx_num);
	rf.write_attribute(bc_grp_id, "bfy_num", md.bfy_num);
	rf.write_attribute(bc_grp_id, "tx_num", md.tx_num);
	rf.write_attribute(bc_grp_id, "ty_num", md.ty_num);
	rf.write_attribute(bc_grp_id, "ax_num", md.ax_num);
	rf.write_attribute(bc_grp_id, "ay_num", md.ay_num);
	rf.write_attribute(bc_grp_id, "vx_num", md.vx_num);
	rf.write_attribute(bc_grp_id, "vy_num", md.vy_num);
	
	// body force at particles
	hid_t bf_dt_id = get_bf_pcl_dt_id();
	BodyForceAtPclData* bfds;
	// bfx
	if (md.bfx_num)
	{
		bfds = new BodyForceAtPclData[md.bfx_num];
		for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
		{
			BodyForceAtPclData& bfd = bfds[bf_id];
			BodyForceAtPcl& bf = md.bfxs[bf_id];
			bfd.pcl_id = bf.pcl_id;
			bfd.bf = bf.bf;
		}
		rf.write_dataset(
			bc_grp_id,
			"BodyForceAtPcl_x",
			md.bfx_num,
			bfds,
			bf_dt_id
		);
		delete[] bfds;
	}
	// bfy
	if (md.bfy_num)
	{
		bfds = new BodyForceAtPclData[md.bfy_num];
		for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
		{
			BodyForceAtPclData& bfd = bfds[bf_id];
			BodyForceAtPcl& bf = md.bfys[bf_id];
			bfd.pcl_id = bf.pcl_id;
			bfd.bf = bf.bf;
		}
		rf.write_dataset(
			bc_grp_id,
			"BodyForceAtPcl_y",
			md.bfy_num,
			bfds,
			bf_dt_id
		);
		delete[] bfds;
	}
	H5Tclose(bf_dt_id);

	// traction at particles
	hid_t tbc_dt_id = get_tbc_pcl_dt_id();
	TractionBCAtPclData* tbcds;
	// tx
	if (md.tx_num)
	{
		tbcds = new TractionBCAtPclData[md.tx_num];
		for (size_t t_id = 0; t_id < md.tx_num; ++t_id)
		{
			TractionBCAtPclData& tbcd = tbcds[t_id];
			TractionBCAtPcl& t = md.txs[t_id];
			tbcd.pcl_id = t.pcl_id;
			tbcd.t = t.t;
		}
		rf.write_dataset(
			bc_grp_id,
			"TractionBCAtPcl_x",
			md.tx_num,
			tbcds,
			tbc_dt_id
		);
		delete[] tbcds;
	}
	// ty
	if (md.ty_num)
	{
		tbcds = new TractionBCAtPclData[md.ty_num];
		for (size_t t_id = 0; t_id < md.ty_num; ++t_id)
		{
			TractionBCAtPclData& tbcd = tbcds[t_id];
			TractionBCAtPcl& t = md.tys[t_id];
			tbcd.pcl_id = t.pcl_id;
			tbcd.t = t.t;
		}
		rf.write_dataset(
			bc_grp_id,
			"TractionBCAtPcl_y",
			md.ty_num,
			tbcds,
			tbc_dt_id
		);
		delete[] tbcds;
	}
	H5Tclose(tbc_dt_id);

	// acceleration bc
	hid_t abc_dt_id = get_abc_dt_id();
	AccelerationBCData* abcds;
	// ax
	if (md.ax_num)
	{
		abcds = new AccelerationBCData[md.ax_num];
		for (size_t a_id = 0; a_id < md.ax_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& a = md.axs[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_x",
			md.ax_num,
			abcds,
			abc_dt_id
		);
		delete[] abcds;
	}
	// ay
	if (md.ay_num)
	{
		abcds = new AccelerationBCData[md.ay_num];
		for (size_t a_id = 0; a_id < md.ay_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& a = md.ays[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_y",
			md.ay_num,
			abcds,
			abc_dt_id
		);
		delete[] abcds;
	}
	H5Tclose(abc_dt_id);

	// velocity bc
	hid_t vbc_dt_id = get_vbc_dt_id();
	VelocityBCData* vbcds;
	// vx
	if (md.vx_num)
	{
		vbcds = new VelocityBCData[md.vx_num];
		for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
		{
			VelocityBCData& vbcd = vbcds[v_id];
			VelocityBC& v = md.vxs[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_x",
			md.vx_num,
			vbcds,
			vbc_dt_id
		);
		delete[] vbcds;
	}
	// vy
	if (md.vy_num)
	{
		vbcds = new VelocityBCData[md.vy_num];
		for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
		{
			VelocityBCData& vbcd = vbcds[v_id];
			VelocityBC& v = md.vys[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_y",
			md.vy_num,
			vbcds,
			vbc_dt_id
		);
		delete[] vbcds;
	}
	H5Tclose(vbc_dt_id);

	rf.close_group(bc_grp_id);
	return 0;
}

int load_boundary_condition_from_hdf5_file(
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;
	
	hid_t md_id = rf.get_model_data_grp_id();

	hid_t bc_id = rf.open_group(md_id, "BoundaryCondition");
	rf.read_attribute(bc_id, "ax_num", md.ax_num);
	rf.read_attribute(bc_id, "ay_num", md.ay_num);
	rf.read_attribute(bc_id, "vx_num", md.vx_num);
	rf.read_attribute(bc_id, "vy_num", md.vy_num);
	rf.read_attribute(bc_id, "tx_num", md.tx_num);
	rf.read_attribute(bc_id, "ty_num", md.ty_num);
	rf.read_attribute(bc_id, "bfx_num", md.bfx_num);
	rf.read_attribute(bc_id, "bfy_num", md.bfy_num);

	// acceleration bcs
	hid_t abc_dt_id = get_abc_dt_id();
	AccelerationBCData* abcds;

	if (md.ax_num)
	{
		md.init_axs(md.ax_num);
		abcds = new AccelerationBCData[md.ax_num];
		rf.read_dataset(bc_id, "AccelerationBC_x", md.ax_num, abcds, abc_dt_id);
		for (size_t a_id = 0; a_id < md.ax_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& abc = md.axs[a_id];
			abcd.to_abc(abc);
		}
		delete[] abcds;
	}
	if (md.ay_num)
	{
		md.init_ays(md.ay_num);
		abcds = new AccelerationBCData[md.ay_num];
		rf.read_dataset(bc_id, "AccelerationBC_y", md.ay_num, abcds, abc_dt_id);
		for (size_t a_id = 0; a_id < md.ay_num; ++a_id)
		{
			AccelerationBCData& abcd = abcds[a_id];
			AccelerationBC& abc = md.ays[a_id];
			abcd.to_abc(abc);
		}
		delete[] abcds;
	}
	H5Tclose(abc_dt_id);

	// velocity bcs
	hid_t vbc_dt_id = get_vbc_dt_id();
	VelocityBCData* vbcds;

	if (md.vx_num)
	{
		md.init_vxs(md.vx_num);
		vbcds = new VelocityBCData[md.vx_num];
		rf.read_dataset(bc_id, "VelocityBC_x", md.vx_num, vbcds, vbc_dt_id);
		for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
		{
			VelocityBCData& vbcd = vbcds[v_id];
			VelocityBC& vbc = md.vxs[v_id];
			vbcd.to_vbc(vbc);
		}
		delete[] vbcds;
	}
	if (md.vy_num)
	{
		md.init_vys(md.vy_num);
		vbcds = new VelocityBCData[md.vy_num];
		rf.read_dataset(bc_id, "VelocityBC_y", md.vy_num, vbcds, vbc_dt_id);
		for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
		{
			VelocityBCData& vbcd = vbcds[v_id];
			VelocityBC& vbc = md.vys[v_id];
			vbcd.to_vbc(vbc);
		}
		delete[] vbcds;
	}
	H5Tclose(vbc_dt_id);

	// tbc
	hid_t tbc_dt_id = get_tbc_pcl_dt_id();
	TractionBCAtPclData* tbcds;
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
			TractionBCAtPclData& tbcd = tbcds[t_id];
			TractionBCAtPcl& tbc = md.tys[t_id];
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
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.create_group(grp_id, "ParticleData");
	
	size_t pcl_num = md.get_pcl_num();
	Model_T2D_ME_s::Particle* pcls = md.get_pcls();
	rf.write_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

	ParticleData* pcl_data = new ParticleData[md.pcl_num];
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		ParticleData& pd = pcl_data[p_id];
		Model_T2D_ME_s::Particle& pcl = md.pcls[p_id];
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
	Model_T2D_ME_s& md,
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
	
	if (res) // if detect error
	{
		delete[] pcls_data;
		return res;
	}

	MatModel::MatModelIdToPointerMap mm_id_map(md);
	md.alloc_pcls(pcl_num);
	Model_T2D_ME_s::Particle* pcls = md.get_pcls();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		ParticleData& pcl_data = pcls_data[p_id];
		Model_T2D_ME_s::Particle& pcl = pcls[p_id];
		MatModel::MaterialModel* pmat = mm_id_map.get_mm_by_id(pcl_data.mat_id);
		if (!pmat)
		{
			throw std::exception("func load_pcl_data_from_hdf5_file error: "
				"particle has no material model.");
		}
		pcl_data.to_pcl(pcl, *pmat);
	}
	delete[] pcls_data;

	rf.close_group(pcl_data_grp_id);
	return 0;
}


int output_material_model_to_hdf5_file(
	Model_T2D_ME_s& md,
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
	Model_T2D_ME_s &md,
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

int output_rigid_circle_to_hdf5_file(
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	if (!md.rigid_circle_is_valid())
		return 0;

	hid_t rc_grp_id = rf.create_group(grp_id, "RigidCircle");

	rf.write_attribute(rc_grp_id, "K_cont", md.K_cont);

	using RigidBody_hdf5_utilities::output_rigid_circle_to_hdf5_file;
	output_rigid_circle_to_hdf5_file(md.rigid_circle, rf, rc_grp_id);

	rf.close_group(rc_grp_id);
	return 0;
}

int load_rigid_circle_from_hdf5_file(
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;
	
	if (!rf.has_group(grp_id, "RigidCircle"))
		return 0;

	hid_t rc_grp_id = rf.open_group(grp_id, "RigidCircle");

	double K_cont;
	rf.read_attribute(rc_grp_id, "K_cont", K_cont);
	md.K_cont = K_cont;

	using RigidBody_hdf5_utilities::load_rigid_circle_from_hdf5_file;
	load_rigid_circle_from_hdf5_file(md.rigid_circle, rf, rc_grp_id);

	md.rigid_circle_is_init = true;
	
	rf.close_group(rc_grp_id);
	return 0;
}

int output_rigid_rect_to_hdf5_file(
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	if (!md.rigid_rect_is_valid())
		return 0;

	hid_t rr_grp_id = rf.create_group(grp_id, "RigidRect");

	rf.write_attribute(rr_grp_id, "K_cont", md.K_cont);

	using RigidBody_hdf5_utilities::output_rigid_rect_to_hdf5_file;
	output_rigid_rect_to_hdf5_file(md.rigid_rect, rf, rr_grp_id);

	rf.close_group(rr_grp_id);
	return 0;
}

int load_rigid_rect_from_hdf5_file(
	Model_T2D_ME_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	if (!rf.has_group(grp_id, "RigidRect"))
		return 0;

	hid_t rr_grp_id = rf.open_group(grp_id, "RigidRect");

	double K_cont;
	rf.read_attribute(rr_grp_id, "K_cont", K_cont);
	md.K_cont = K_cont;

	using RigidBody_hdf5_utilities::load_rigid_rect_from_hdf5_file;
	load_rigid_rect_from_hdf5_file(md.rigid_rect, rf, rr_grp_id);

	md.rigid_rect_is_init = true;

	rf.close_group(rr_grp_id);
	return 0;
}

// output the whole model to ModelData
int output_model_to_hdf5_file(
	Model_T2D_ME_s& md,
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
	output_rigid_rect_to_hdf5_file(md, rf, md_grp_id);
	return 0;
}

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(
	Model_T2D_ME_s& md,
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
	output_rigid_rect_to_hdf5_file(md, rf, frame_grp_id);
	return 0;
}

// load model data from hdf5 to model data
int load_me_s_model_from_hdf5_file(
	Model_T2D_ME_s& md,
	Step_T2D_ME_s& step,
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

	// material model
	load_material_model_from_hdf5_file(md, rf, th_frame_id);
	// particle data
	load_pcl_data_from_hdf5_file(md, rf, th_frame_id);

	// rigid object
	load_rigid_circle_from_hdf5_file(md, rf, th_frame_id);
	load_rigid_rect_from_hdf5_file(md, rf, th_frame_id);

	step.is_first_step = false;
	rf.read_attribute(th_frame_id, "total_time", step.start_time);
	rf.read_attribute(th_frame_id, "total_substep_num", step.prev_substep_num);

	rf.close_group(th_frame_id);
	rf.close_group(th_id);

	return 0;
}

};
