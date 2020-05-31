#include "Simulations_pcp.h"

#include "Model_hdf5_utilities.h"
#include "Model_T3D_CHM_s_hdf5_utilities.h"

namespace Model_T3D_CHM_s_hdf5_utilities
{
using namespace Model_hdf5_utilities;

int output_background_mesh_to_hdf5_file(
	Model_T3D_CHM_s &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	hid_t md_id = rf.get_model_data_grp_id();
	
	// fluid properties
	rf.write_attribute(md_id, "Kf", md.Kf);
	rf.write_attribute(md_id, "k", md.k);
	rf.write_attribute(md_id, "miu", md.miu);

	hid_t bg_mesh_id = rf.create_group(md_id, "BackgroundMesh");

	// bg mesh attributes
	size_t node_num = md.get_node_num();
	Model_T3D_CHM_s::Node *nodes = md.get_nodes();
	size_t elem_num = md.get_elem_num();
	Model_T3D_CHM_s::Element *elems = md.get_elems();
	const char *bg_mesh_type = "T3D";
	rf.write_attribute(bg_mesh_id, "type", strlen(bg_mesh_type), bg_mesh_type);
	rf.write_attribute(bg_mesh_id, "node_num", node_num);
	rf.write_attribute(bg_mesh_id, "element_num", elem_num);
	// node coordinates
	Node3DData *nodes_data = new Node3DData[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node3DData &node_data = nodes_data[n_id];
		Model_T3D_CHM_s::Node n = nodes[n_id];
		node_data.id = n.id;
		node_data.x = n.x;
		node_data.y = n.y;
		node_data.z = n.z;
	}
	hid_t nd_dt_id = get_nd_3d_dt_id();
	rf.write_dataset(
		bg_mesh_id,
		"NodeCoordinate",
		node_num,
		nodes_data,
		nd_dt_id
	);
	H5Tclose(nd_dt_id);
	delete[] nodes_data;
	// element indices
	Elem3DData *elems_data = new Elem3DData[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Elem3DData &elem_data = elems_data[e_id];
		Model_T3D_CHM_s::Element &e = elems[e_id];
		elem_data.id = e.id;
		elem_data.n1 = e.n1;
		elem_data.n2 = e.n2;
		elem_data.n3 = e.n3;
		elem_data.n4 = e.n4;
	}
	hid_t ed_dt_id = get_ed_3d_dt_id();
	rf.write_dataset(
		bg_mesh_id,
		"ElementTopology",
		elem_num,
		elems_data,
		ed_dt_id
	);
	H5Tclose(ed_dt_id);
	delete[] elems_data;
	rf.close_group(bg_mesh_id);

	return 0;
}

int load_background_mesh_from_hdf5_file(
	Model_T3D_CHM_s &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	// model data output
	hid_t md_out_id = rf.get_model_data_grp_id();

	// liquid properties
	rf.read_attribute(md_out_id, "Kf", md.Kf);
	rf.read_attribute(md_out_id, "k", md.k);
	rf.read_attribute(md_out_id, "miu", md.miu);

	hid_t bg_mesh_id = rf.open_group(md_out_id, "BackgroundMesh");
	size_t node_num, elem_num;
	rf.read_attribute(md_out_id, "node_num", node_num);
	md.alloc_nodes(node_num);
	Model_T3D_CHM_s::Node *nodes = md.get_nodes();
	rf.read_attribute(md_out_id, "elem_num", elem_num);
	md.alloc_elements(elem_num);
	Model_T3D_CHM_s::Element *elems = md.get_elems();

	// nodes
	Node3DData *nodes_data = new Node3DData[node_num];
	hid_t nd_dt_id = get_nd_3d_dt_id();
	rf.read_dataset(
		bg_mesh_id,
		"NodeCoordinate",
		node_num,
		nodes_data,
		nd_dt_id
	);
	H5Tclose(nd_dt_id);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node3DData &node_data = nodes_data[n_id];
		Model_T3D_CHM_s::Node &n = nodes[n_id];
		n.id = node_data.id;
		n.x = node_data.x;
		n.y = node_data.y;
		n.z = node_data.z;
	}
	delete[] nodes_data;

	// elements
	Elem3DData *elems_data = new Elem3DData[elem_num];
	hid_t ed_dt_id = get_ed_3d_dt_id();
	rf.read_dataset(
		bg_mesh_id,
		"ElementTopology",
		elem_num,
		elems_data,
		ed_dt_id
	);
	H5Tclose(ed_dt_id);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Elem3DData &elem_data = elems_data[e_id];
		Model_T3D_CHM_s::Element &elem = elems[e_id];
		elem.id = elem_data.id;
		elem.n1 = elem_data.n1;
		elem.n2 = elem_data.n2;
		elem.n3 = elem_data.n3;
		elem.n4 = elem_data.n4;
	}
	delete[] elems_data;

	rf.close_group(bg_mesh_id);
	return 0;
}

int output_boundary_condition_to_hdf5_file(
	Model_T3D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t bc_grp_id = rf.create_group(grp_id, "BoundaryCondition");

	rf.write_attribute(bc_grp_id, "bfx_num", md.bfx_num);
	rf.write_attribute(bc_grp_id, "bfy_num", md.bfy_num);
	rf.write_attribute(bc_grp_id, "bfz_num", md.bfz_num);
	rf.write_attribute(bc_grp_id, "tx_num", md.tx_num);
	rf.write_attribute(bc_grp_id, "ty_num", md.ty_num);
	rf.write_attribute(bc_grp_id, "tz_num", md.tz_num);
	rf.write_attribute(bc_grp_id, "asx_num", md.asx_num);
	rf.write_attribute(bc_grp_id, "asy_num", md.asy_num);
	rf.write_attribute(bc_grp_id, "asz_num", md.asz_num);
	rf.write_attribute(bc_grp_id, "vsx_num", md.vsx_num);
	rf.write_attribute(bc_grp_id, "vsy_num", md.vsy_num);
	rf.write_attribute(bc_grp_id, "vsz_num", md.vsz_num);
	rf.write_attribute(bc_grp_id, "afx_num", md.afx_num);
	rf.write_attribute(bc_grp_id, "afy_num", md.afy_num);
	rf.write_attribute(bc_grp_id, "afz_num", md.afz_num);
	rf.write_attribute(bc_grp_id, "vfx_num", md.vfx_num);
	rf.write_attribute(bc_grp_id, "vfy_num", md.vfy_num);
	rf.write_attribute(bc_grp_id, "vfz_num", md.vfz_num);

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
	// bfz
	if (md.bfz_num)
	{
		bfds = new BodyForceAtPclData[md.bfz_num];
		for (size_t bf_id = 0; bf_id < md.bfz_num; ++bf_id)
		{
			BodyForceAtPclData& bfd = bfds[bf_id];
			BodyForceAtPcl& bf = md.bfzs[bf_id];
			bfd.pcl_id = bf.pcl_id;
			bfd.bf = bf.bf;
		}
		rf.write_dataset(
			bc_grp_id,
			"BodyForceAtPcl_z",
			md.bfz_num,
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
	// tz
	if (md.tz_num)
	{
		tbcds = new TractionBCAtPclData[md.tz_num];
		for (size_t t_id = 0; t_id < md.tz_num; ++t_id)
		{
			TractionBCAtPclData& tbcd = tbcds[t_id];
			TractionBCAtPcl& t = md.tzs[t_id];
			tbcd.pcl_id = t.pcl_id;
			tbcd.t = t.t;
		}
		rf.write_dataset(
			bc_grp_id,
			"TractionBCAtPcl_z",
			md.tz_num,
			tbcds,
			tbc_dt_id
		);
		delete[] tbcds;
	}
	H5Tclose(tbc_dt_id);

	hid_t abc_dt_id = get_abc_dt_id();
	hid_t vbc_dt_id = get_vbc_dt_id();

	// solid acceleration bc
	AccelerationBCData* asbcds;
	// ax
	if (md.asx_num)
	{
		asbcds = new AccelerationBCData[md.asx_num];
		for (size_t a_id = 0; a_id < md.asx_num; ++a_id)
		{
			AccelerationBCData& abcd = asbcds[a_id];
			AccelerationBC& a = md.asxs[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_sx",
			md.asx_num,
			asbcds,
			abc_dt_id
		);
		delete[] asbcds;
	}
	// ay
	if (md.asy_num)
	{
		asbcds = new AccelerationBCData[md.asy_num];
		for (size_t a_id = 0; a_id < md.asy_num; ++a_id)
		{
			AccelerationBCData& abcd = asbcds[a_id];
			AccelerationBC& a = md.asys[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_sy",
			md.asy_num,
			asbcds,
			abc_dt_id
		);
		delete[] asbcds;
	}
	// az
	if (md.asz_num)
	{
		asbcds = new AccelerationBCData[md.asz_num];
		for (size_t a_id = 0; a_id < md.asz_num; ++a_id)
		{
			AccelerationBCData& abcd = asbcds[a_id];
			AccelerationBC& a = md.aszs[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_sz",
			md.asz_num,
			asbcds,
			abc_dt_id
		);
		delete[] asbcds;
	}

	// solid velocity bc
	VelocityBCData *vsbcds;
	// vx
	if (md.vsx_num)
	{
		vsbcds = new VelocityBCData[md.vsx_num];
		for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
		{
			VelocityBCData& vbcd = vsbcds[v_id];
			VelocityBC& v = md.vsxs[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_sx",
			md.vsx_num,
			vsbcds,
			vbc_dt_id
		);
		delete[] vsbcds;
	}
	// vy
	if (md.vsy_num)
	{
		vsbcds = new VelocityBCData[md.vsy_num];
		for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
		{
			VelocityBCData& vbcd = vsbcds[v_id];
			VelocityBC& v = md.vsys[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_sy",
			md.vsy_num,
			vsbcds,
			vbc_dt_id
		);
		delete[] vsbcds;
	}
	// vz
	if (md.vsz_num)
	{
		vsbcds = new VelocityBCData[md.vsz_num];
		for (size_t v_id = 0; v_id < md.vsz_num; ++v_id)
		{
			VelocityBCData& vbcd = vsbcds[v_id];
			VelocityBC& v = md.vszs[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_sz",
			md.vsz_num,
			vsbcds,
			vbc_dt_id
		);
		delete[] vsbcds;
	}

	// fluid acceleration bc
	AccelerationBCData* afbcds;
	// ax
	if (md.afx_num)
	{
		afbcds = new AccelerationBCData[md.afx_num];
		for (size_t a_id = 0; a_id < md.afx_num; ++a_id)
		{
			AccelerationBCData& abcd = afbcds[a_id];
			AccelerationBC& a = md.afxs[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_fx",
			md.afx_num,
			afbcds,
			abc_dt_id
		);
		delete[] afbcds;
	}
	// ay
	if (md.afy_num)
	{
		afbcds = new AccelerationBCData[md.afy_num];
		for (size_t a_id = 0; a_id < md.afy_num; ++a_id)
		{
			AccelerationBCData& abcd = afbcds[a_id];
			AccelerationBC& a = md.afys[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_fy",
			md.afy_num,
			afbcds,
			abc_dt_id
		);
		delete[] afbcds;
	}
	// az
	if (md.afz_num)
	{
		afbcds = new AccelerationBCData[md.afz_num];
		for (size_t a_id = 0; a_id < md.afz_num; ++a_id)
		{
			AccelerationBCData& abcd = afbcds[a_id];
			AccelerationBC& a = md.afzs[a_id];
			abcd.node_id = a.node_id;
			abcd.a = a.a;
		}
		rf.write_dataset(
			bc_grp_id,
			"AccelerationBC_fz",
			md.afz_num,
			afbcds,
			abc_dt_id
		);
		delete[] afbcds;
	}

	// fluid velocity bc
	VelocityBCData* vfbcds;
	// vx
	if (md.vfx_num)
	{
		vfbcds = new VelocityBCData[md.vfx_num];
		for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
		{
			VelocityBCData& vbcd = vfbcds[v_id];
			VelocityBC& v = md.vfxs[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_fx",
			md.vfx_num,
			vfbcds,
			vbc_dt_id
		);
		delete[] vfbcds;
	}
	// vy
	if (md.vfy_num)
	{
		vfbcds = new VelocityBCData[md.vfy_num];
		for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
		{
			VelocityBCData& vbcd = vfbcds[v_id];
			VelocityBC& v = md.vfys[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_fy",
			md.vfy_num,
			vfbcds,
			vbc_dt_id
		);
		delete[] vfbcds;
	}
	// vz
	if (md.vfz_num)
	{
		vfbcds = new VelocityBCData[md.vfz_num];
		for (size_t v_id = 0; v_id < md.vfz_num; ++v_id)
		{
			VelocityBCData& vbcd = vfbcds[v_id];
			VelocityBC& v = md.vfzs[v_id];
			vbcd.node_id = v.node_id;
			vbcd.v = v.v;
		}
		rf.write_dataset(
			bc_grp_id,
			"VelocityBC_fz",
			md.vfz_num,
			vfbcds,
			vbc_dt_id
		);
		delete[] vfbcds;
	}

	H5Tclose(abc_dt_id);
	H5Tclose(vbc_dt_id);

	rf.close_group(bc_grp_id);
	return 0;
}

int load_boundary_condition_from_hdf5_file(
	Model_T3D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	// to be finished...

	return 0;
}

int output_pcl_data_to_hdf5_file(
	Model_T3D_CHM_s &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	if (grp_id < 0)
		return -1;

	hid_t pcl_data_grp_id = rf.create_group(grp_id, "ParticleData");

	size_t pcl_num = md.get_pcl_num();
	Model_T3D_CHM_s::Particle* pcls = md.get_pcls();
	rf.write_attribute(pcl_data_grp_id, "pcl_num", pcl_num);

	ParticleData* pcl_data = new ParticleData[pcl_num];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		ParticleData &pd = pcl_data[p_id];
		Model_T3D_CHM_s::Particle& pcl = pcls[p_id];
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
	Model_T3D_CHM_s &md,
	ResultFile_hdf5 &rf,
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
		grp_id,
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

	md.alloc_pcls(pcl_num);
	Model_T3D_CHM_s::Particle* pcls = md.get_pcls();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		ParticleData& pcl_data = pcls_data[p_id];
		Model_T3D_CHM_s::Particle& pcl = pcls[p_id];
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
	{ return static_cast<Particle *>(ext_data)->id;	}
};

int output_material_model_to_hdf5_file(
	Model_T3D_CHM_s &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	using namespace MatModel;

	size_t mm_id, mm_num;
	hid_t mm_grp_id = rf.create_group(grp_id, "MaterialModel");

	// linear elasticity
	mm_num = md.get_num_LinearElasticity();
	if (mm_num)
	{
		rf.write_attribute(mm_grp_id, "LinearElasticity_num", mm_num);
		
		LinearElasticityStateData *mm_data = new LinearElasticityStateData[mm_num];
		mm_id = 0;
		for (LinearElasticity *iter = md.first_LinearElasticity(); 
			 md.is_not_end_LinearElasticity(iter);
			 iter = md.next_LinearElasticity(iter))
		{
			mm_data[mm_id].pcl_id = get_pcl_id<Model_T3D_CHM_s::Particle>(iter->ext_data);
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
		
		ModifiedCamClayStateData *mm_data = new ModifiedCamClayStateData[mm_num];
		mm_id = 0;
		for (ModifiedCamClay *iter = md.first_ModifiedCamClay();
			 md.is_not_end_ModifiedCamClay(iter);
			 iter = md.next_ModifiedCamClay(iter))
		{
			mm_data[mm_id].pcl_id = get_pcl_id<Model_T3D_CHM_s::Particle>(iter->ext_data);
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
	Model_T3D_CHM_s &md,
	ResultFile_hdf5 &rf,
	hid_t grp_id
	)
{
	using namespace MatModel;
	hid_t mm_dset_id;
	size_t mm_num;

	hid_t mm_grp_id = rf.open_group(grp_id, "MaterialModel");

	// linear elasticity
	if (rf.has_dataset(mm_grp_id, "LinearElasticity"))
	{
		rf.read_attribute(mm_grp_id, "LinearElasticity_num", mm_num);

		// get data
		LinearElasticityStateData *mm_data = new LinearElasticityStateData[mm_num];
		hid_t le_dt_id = get_le_hdf5_dt_id();
		rf.read_dataset(
			mm_grp_id,
			"LinearElasticity",
			mm_num,
			mm_data,
			le_dt_id
			);
		H5Tclose(le_dt_id);
		LinearElasticity *mms = md.add_LinearElasticity(mm_num);
		for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
		{
			LinearElasticityStateData &mmd = mm_data[mm_id];
			LinearElasticity &mm = mms[mm_id];
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
		ModifiedCamClayStateData *mm_data = new ModifiedCamClayStateData[mm_num];
		hid_t mcc_dt_id = get_mcc_hdf5_dt_id();
		rf.read_dataset(
			mm_grp_id,
			"ModifiedCamClay",
			mm_num,
			mm_data,
			mcc_dt_id
			);
		H5Tclose(mcc_dt_id);
		ModifiedCamClay *mms = md.add_ModifiedCamClay(mm_num);
		for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
		{
			ModifiedCamClayStateData &mmd = mm_data[mm_id];
			ModifiedCamClay &mm = mms[mm_id];
			mmd.to_mm(mm);
			md.pcls[mmd.pcl_id].set_mat_model(mm);
		}
		delete[] mm_data;
	}

	rf.close_group(mm_grp_id);
	return 0;
}

int output_model_to_hdf5_file(
	Model_T3D_CHM_s& md,
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
	return 0;
}

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(
	Model_T3D_CHM_s& md,
	ResultFile_hdf5& rf,
	hid_t frame_grp_id
	)
{
	// particle data
	output_pcl_data_to_hdf5_file(md, rf, frame_grp_id);
	// material model
	output_material_model_to_hdf5_file(md, rf, frame_grp_id);
	return 0;
}

int load_model_from_hdf5_file(
	Model_T3D_CHM_s &md,
	const char *hdf5_name,
	const char *th_name,
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
	// constitutive model data
	load_material_model_from_hdf5_file(md, rf, th_frame_id);
	rf.close_group(th_frame_id);
	rf.close_group(th_id);

	return 0;
}

};