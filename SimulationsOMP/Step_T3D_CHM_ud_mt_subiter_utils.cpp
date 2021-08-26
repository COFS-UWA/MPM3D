#include "SimulationsOMP_pcp.h"

#include <omp.h>

#include "Step_T3D_CHM_ud_mt_subiter.h"

#include <fstream>
extern std::fstream t3d_chm_ud_mt_subit_db_file;

#define one_fourth (0.25)
#define one_third (1.0/3.0)
#define N_min (1.0e-10)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

int Step_T3D_CHM_ud_mt_subiter::subiteration(
	size_t my_th_id,
	size_t p_id0, size_t p_id1,
	size_t ve_id0, size_t ve_id1,
	const size_t* my_valid_elem_ids, size_t my_valid_elem_num,
	SortedPclVarArrays& spva0,
	const size_t *pcl_in_elem0,
	const size_t *node_has_elem0, const size_t *node_elem_pair0)
{
	Model_T3D_CHM_mt& md = *(Model_T3D_CHM_mt*)(model);

	ThreadData& thd = thread_datas[my_th_id];
	size_t* const pcl_index0 = spva0.pcl_index;
	double* const pcl_n0 = spva0.pcl_n;
	double* const pcl_density_f0 = spva0.pcl_density_f;
	Displacement* const pcl_u_s0 = spva0.pcl_u_s;
	Displacement* const pcl_u_f0 = spva0.pcl_u_f;
	Velocity* const pcl_v_s0 = spva0.pcl_v_s;
	Velocity* const pcl_v_f0 = spva0.pcl_v_f;
	Stress* const pcl_stress0 = spva0.pcl_stress;
	double* const pcl_p0 = spva0.pcl_p;
	Strain* const pcl_strain0 = spva0.pcl_strain;
	Strain* const pcl_estrain0 = spva0.pcl_estrain;
	Strain* const pcl_pstrain0 = spva0.pcl_pstrain;
	ShapeFunc* const pcl_N0 = spva0.pcl_N;

#pragma omp master
	{
		// init rigid body
		cf_tmp.reset();

		// init subiteration
		cur_e_kin = 0.0;
	}

#pragma omp barrier

	size_t e_id, p_id, ori_p_id;
	double p_vol;
	double e_p_vol = 0.0;
	double e_s11 = 0.0;
	double e_s22 = 0.0;
	double e_s33 = 0.0;
	double e_s12 = 0.0;
	double e_s23 = 0.0;
	double e_s31 = 0.0;
	double e_p;
	double one_fourth_bfx_s;
	double one_fourth_bfy_s;
	double one_fourth_bfz_s;
	double en1_fx_ext = 0.0;
	double en1_fy_ext = 0.0;
	double en1_fz_ext = 0.0;
	double en2_fx_ext = 0.0;
	double en2_fy_ext = 0.0;
	double en2_fz_ext = 0.0;
	double en3_fx_ext = 0.0;
	double en3_fy_ext = 0.0;
	double en3_fz_ext = 0.0;
	double en4_fx_ext = 0.0;
	double en4_fy_ext = 0.0;
	double en4_fz_ext = 0.0;
	e_id = pcl_in_elem0[p_id0];
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		// solid external load
		ori_p_id = pcl_index0[p_id];
		ShapeFunc& p_N0 = pcl_N0[p_id];
		const Force& p_bf_s = pcl_bf_s[ori_p_id];
		const Force& p_bf_f = pcl_bf_f[ori_p_id];
		one_fourth_bfx_s = one_fourth * (p_bf_s.fx + p_bf_f.fx);
		one_fourth_bfy_s = one_fourth * (p_bf_s.fy + p_bf_f.fy);
		one_fourth_bfz_s = one_fourth * (p_bf_s.fz + p_bf_f.fz);
		const Force& p_t = pcl_t[ori_p_id];
		en1_fx_ext += one_fourth_bfx_s + p_N0.N1 * p_t.fx;
		en1_fy_ext += one_fourth_bfy_s + p_N0.N1 * p_t.fy;
		en1_fz_ext += one_fourth_bfz_s + p_N0.N1 * p_t.fz;
		en2_fx_ext += one_fourth_bfx_s + p_N0.N2 * p_t.fx;
		en2_fy_ext += one_fourth_bfy_s + p_N0.N2 * p_t.fy;
		en2_fz_ext += one_fourth_bfz_s + p_N0.N2 * p_t.fz;
		en3_fx_ext += one_fourth_bfx_s + p_N0.N3 * p_t.fx;
		en3_fy_ext += one_fourth_bfy_s + p_N0.N3 * p_t.fy;
		en3_fz_ext += one_fourth_bfz_s + p_N0.N3 * p_t.fz;
		en4_fx_ext += one_fourth_bfx_s + p_N0.N4 * p_t.fx;
		en4_fy_ext += one_fourth_bfy_s + p_N0.N4 * p_t.fy;
		en4_fz_ext += one_fourth_bfz_s + p_N0.N4 * p_t.fz;

		// map pcl mass and volume
		p_vol = pcl_vol[p_id];
		e_p_vol += p_vol;

		// map stress
		Stress& p_s0 = pcl_stress0[p_id];
		e_s11 += p_s0.s11 * p_vol;
		e_s22 += p_s0.s22 * p_vol;
		e_s33 += p_s0.s33 * p_vol;
		e_s12 += p_s0.s12 * p_vol;
		e_s23 += p_s0.s23 * p_vol;
		e_s31 += p_s0.s31 * p_vol;

		if (e_id != pcl_in_elem0[p_id + 1])
		{
			// external nodal force
			// node 1
			Force& en1_f_ext = elem_node_f_ext[e_id * 4];
			en1_f_ext.fx = en1_fx_ext;
			en1_f_ext.fy = en1_fy_ext;
			en1_f_ext.fz = en1_fz_ext;
			// node 2
			Force& en2_f_ext = elem_node_f_ext[e_id * 4 + 1];
			en2_f_ext.fx = en2_fx_ext;
			en2_f_ext.fy = en2_fy_ext;
			en2_f_ext.fz = en2_fz_ext;
			// node 3
			Force& en3_f_ext = elem_node_f_ext[e_id * 4 + 2];
			en3_f_ext.fx = en3_fx_ext;
			en3_f_ext.fy = en3_fy_ext;
			en3_f_ext.fz = en3_fz_ext;
			// node 4
			Force& en4_f_ext = elem_node_f_ext[e_id * 4 + 3];
			en4_f_ext.fx = en4_fx_ext;
			en4_f_ext.fy = en4_fy_ext;
			en4_f_ext.fz = en4_fz_ext;
			
			e_s11 /= e_p_vol;
			e_s22 /= e_p_vol;
			e_s33 /= e_p_vol;
			e_s12 /= e_p_vol;
			e_s23 /= e_p_vol;
			e_s31 /= e_p_vol;
			e_p = elem_p[e_id];

			const double e_p_int_vol = elem_pcl_int_vol[e_id] < elem_vol[e_id] ? elem_pcl_int_vol[e_id] : elem_vol[e_id];
			const DShapeFuncABC& e_dN = elem_N_abc[e_id];
			// node 1
			Force& en1_f_int = elem_node_f_int[e_id * 4];
			en1_f_int.fx = (e_dN.dN1_dx * (e_s11 - e_p) + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_int_vol;
			en1_f_int.fy = (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * (e_s22 - e_p) + e_dN.dN1_dz * e_s23) * e_p_int_vol;
			en1_f_int.fz = (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * (e_s33 - e_p)) * e_p_int_vol;
			// node 2
			Force& en2_f_int = elem_node_f_int[e_id * 4 + 1];
			en2_f_int.fx = (e_dN.dN2_dx * (e_s11 - e_p) + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_int_vol;
			en2_f_int.fy = (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * (e_s22 - e_p) + e_dN.dN2_dz * e_s23) * e_p_int_vol;
			en2_f_int.fz = (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * (e_s33 - e_p)) * e_p_int_vol;
			// node 3
			Force& en3_f_int = elem_node_f_int[e_id * 4 + 2];
			en3_f_int.fx = (e_dN.dN3_dx * (e_s11 - e_p) + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_int_vol;
			en3_f_int.fy = (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * (e_s22 - e_p) + e_dN.dN3_dz * e_s23) * e_p_int_vol;
			en3_f_int.fz = (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * (e_s33 - e_p)) * e_p_int_vol;
			// node 4
			Force& en4_f_int = elem_node_f_int[e_id * 4 + 3];
			en4_f_int.fx = (e_dN.dN4_dx * (e_s11 - e_p) + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_int_vol;
			en4_f_int.fy = (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * (e_s22 - e_p) + e_dN.dN4_dz * e_s23) * e_p_int_vol;
			en4_f_int.fz = (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * (e_s33 - e_p)) * e_p_int_vol;

			e_id = pcl_in_elem0[p_id + 1];
			assert(e_id < elem_num || e_id == SIZE_MAX);

			e_p_vol = 0.0;
			e_s11 = 0.0;
			e_s22 = 0.0;
			e_s33 = 0.0;
			e_s12 = 0.0;
			e_s23 = 0.0;
			e_s31 = 0.0;
			en1_fx_ext = 0.0;
			en1_fy_ext = 0.0;
			en1_fz_ext = 0.0;
			en2_fx_ext = 0.0;
			en2_fy_ext = 0.0;
			en2_fz_ext = 0.0;
			en3_fx_ext = 0.0;
			en3_fy_ext = 0.0;
			en3_fz_ext = 0.0;
			en4_fx_ext = 0.0;
			en4_fy_ext = 0.0;
			en4_fz_ext = 0.0;
		}
	}

	Force3D rc_force;

	if (md.has_rigid_cylinder())
	{
		rc_force.reset();
		apply_rigid_cylinder(
			p_id0, p_id1,
			pcl_in_elem0,
			spva0, rc_force,
			substep_index, thd);

#pragma omp critical
		cf_tmp.combine(rc_force);
	}

	if (md.has_t3d_rigid_mesh())
	{
		rc_force.reset();
		apply_t3d_rigid_mesh(
			p_id0, p_id1,
			pcl_in_elem0,
			spva0, rc_force,
			substep_index, thd);

#pragma omp critical
		cf_tmp.combine(rc_force);
	}

#pragma omp barrier
	// update node variables
	size_t n_id, ve_id, ne_id, bc_mask;
	double n_am, n_pam;
	double n_fx = 0.0;
	double n_fy = 0.0;
	double n_fz = 0.0;
	Acceleration n_pa;
	double vbc_len;
	double e_kin = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		ne_id = node_elem_pair0[ve_id];
		assert(ne_id < elem_num * 4);

		e_id = ne_id / 4;
		const Force& nf_ext = elem_node_f_ext[ne_id];
		const Force& nf_int = elem_node_f_int[ne_id];
		n_fx += nf_ext.fx - nf_int.fx;
		n_fy += nf_ext.fy - nf_int.fy;
		n_fz += nf_ext.fz - nf_int.fz;

		if (n_id != node_has_elem0[ve_id + 1])
		{
			// solid
			n_am = node_am[n_id];
			n_pam = n_am * mass_factor;
			Acceleration& n_a = node_a[n_id];
			n_pa.ax = (n_fx - n_am * n_a.ax) / n_pam;
			n_pa.ay = (n_fy - n_am * n_a.ay) / n_pam;
			n_pa.az = (n_fz - n_am * n_a.az) / n_pam;
			Velocity& n_pv = node_pv[n_id];
			n_pv.vx += n_pa.ax * pdt;
			n_pv.vy += n_pa.ay * pdt;
			n_pv.vz += n_pa.az * pdt;
			// apply bcs
			NodeVBCVec& n_vbc_v = node_vbc_vec[n_id];
			vbc_len = n_pa.ax * n_vbc_v.x + n_pa.ay * n_vbc_v.y + n_pa.az * n_vbc_v.z;
			n_pa.ax -= vbc_len * n_vbc_v.x;
			n_pa.ay -= vbc_len * n_vbc_v.y;
			n_pa.az -= vbc_len * n_vbc_v.z;
			vbc_len = n_pv.vx * n_vbc_v.x + n_pv.vy * n_vbc_v.y + n_pv.vz * n_vbc_v.z;
			n_pv.vx -= vbc_len * n_vbc_v.x;
			n_pv.vy -= vbc_len * n_vbc_v.y;
			n_pv.vz -= vbc_len * n_vbc_v.z;
			NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
			bc_mask = SIZE_MAX + size_t(n_has_vbc.has_vx_bc);
			n_pa.iax &= bc_mask;
			n_pv.ivx &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc.has_vy_bc);
			n_pa.iay &= bc_mask;
			n_pv.ivy &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc.has_vz_bc);
			n_pa.iaz &= bc_mask;
			n_pv.ivz &= bc_mask;

			e_kin += n_pam * (n_pv.vx * n_pv.vx + n_pv.vy * n_pv.vy + n_pv.vz * n_pv.vz);

			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < node_num || n_id == SIZE_MAX);

			n_fx = 0.0;
			n_fy = 0.0;
			n_fz = 0.0;
		}
	}

#pragma omp critical
	cur_e_kin += e_kin;

#pragma omp barrier
	// judge convergence
#pragma omp master
	{
		if (cur_e_kin < prev_e_kin)
		{
			cal_status = 1;
			//if (prev_e_kin < converge_e_kin_ratio * max_e_kin)
			if (prev_e_kin < 1.0e-11)
				cal_status = 2; // converge
			if (max_e_kin < prev_e_kin)
				max_e_kin = prev_e_kin;
			//cur_e_kin = 0.0;
			prev_e_kin = 0.0;
		}
		else
		{
			cal_status = 0;
			prev_e_kin = cur_e_kin;
		}
		//t3d_chm_ud_mt_subit_db_file << cur_e_kin << ", " << prev_e_kin << ", "
		//	<< max_e_kin << "\n";
	}

#pragma omp barrier
	if (cal_status == 0) // normal cal, update displacement
	{
		n_id = node_has_elem0[ve_id0];
		assert(n_id < node_num);
		for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			if (n_id != node_has_elem0[ve_id + 1])
			{
				Velocity& n_pv = node_pv[n_id];
				Displacement& n_pdu = node_pdu[n_id];
				n_pdu.ux = n_pv.vx * pdt;
				n_pdu.uy = n_pv.vy * pdt;
				n_pdu.uz = n_pv.vz * pdt;

				// modify velocity and acceleration
				Displacement& n_du = node_du[n_id];
				n_du.ux += n_pdu.ux;
				n_du.uy += n_pdu.uy;
				n_du.uz += n_pdu.uz;
				Velocity& n_v = node_v[n_id];
				n_v.vx = n_du.ux / dtime;
				n_v.vy = n_du.uy / dtime;
				n_v.vz = n_du.uz / dtime;
				Acceleration& n_a = node_a[n_id];
				Velocity& n_vn = node_vn[n_id];
				n_a.ax = (n_v.vx - n_vn.vx) / dtime;
				n_a.ay = (n_v.vy - n_vn.vy) / dtime;
				n_a.az = (n_v.vz - n_vn.vz) / dtime;

				n_id = node_has_elem0[ve_id + 1];
				assert(n_id < node_num || n_id == SIZE_MAX);
			}
		}
	}
	else if (cal_status == 1) // reset pseudo velocity and continue
	{
		n_id = node_has_elem0[ve_id0];
		assert(n_id < node_num);
		for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			if (n_id != node_has_elem0[ve_id + 1])
			{
				Velocity& n_pv = node_pv[n_id];
				n_pv.vx = 0.0;
				n_pv.vy = 0.0;
				n_pv.vz = 0.0;
				n_id = node_has_elem0[ve_id + 1];
				assert(n_id < node_num || n_id == SIZE_MAX);
			}
		}
		return 0;
	}
	else if (cal_status == 2) // converge and exit
	{
		return 1;
	}

#pragma omp barrier
	// cal element strain and "enhancement"
	double e_de_vol, e_pde_vol;
	for (ve_id = 0; ve_id < my_valid_elem_num; ++ve_id)
	{
		e_id = my_valid_elem_ids[ve_id];
		assert(e_id < elem_num);

		const ElemNodeIndex& eni = elem_node_id[e_id];
		const Displacement &n1_du = node_du[eni.n1];
		const Displacement &n2_du = node_du[eni.n2];
		const Displacement &n3_du = node_du[eni.n3];
		const Displacement &n4_du = node_du[eni.n4];
		const DShapeFuncABC& e_dN = elem_N_abc[e_id];
		StrainInc& e_de = elem_de[e_id];
		e_de.de11 = e_dN.dN1_dx * n1_du.ux + e_dN.dN2_dx * n2_du.ux + e_dN.dN3_dx * n3_du.ux + e_dN.dN4_dx * n4_du.ux;
		e_de.de22 = e_dN.dN1_dy * n1_du.uy + e_dN.dN2_dy * n2_du.uy + e_dN.dN3_dy * n3_du.uy + e_dN.dN4_dy * n4_du.uy;
		e_de.de33 = e_dN.dN1_dz * n1_du.uz + e_dN.dN2_dz * n2_du.uz + e_dN.dN3_dz * n3_du.uz + e_dN.dN4_dz * n4_du.uz;
		e_de.de12 = (e_dN.dN1_dx * n1_du.uy + e_dN.dN2_dx * n2_du.uy + e_dN.dN3_dx * n3_du.uy + e_dN.dN4_dx * n4_du.uy
				   + e_dN.dN1_dy * n1_du.ux + e_dN.dN2_dy * n2_du.ux + e_dN.dN3_dy * n3_du.ux + e_dN.dN4_dy * n4_du.ux) * 0.5;
		e_de.de23 = (e_dN.dN1_dy * n1_du.uz + e_dN.dN2_dy * n2_du.uz + e_dN.dN3_dy * n3_du.uz + e_dN.dN4_dy * n4_du.uz
				   + e_dN.dN1_dz * n1_du.uy + e_dN.dN2_dz * n2_du.uy + e_dN.dN3_dz * n3_du.uy + e_dN.dN4_dz * n4_du.uy) * 0.5;
		e_de.de31 = (e_dN.dN1_dz * n1_du.ux + e_dN.dN2_dz * n2_du.ux + e_dN.dN3_dz * n3_du.ux + e_dN.dN4_dz * n4_du.ux
				   + e_dN.dN1_dx * n1_du.uz + e_dN.dN2_dx * n2_du.uz + e_dN.dN3_dx * n3_du.uz + e_dN.dN4_dx * n4_du.uz) * 0.5;
		e_de_vol = e_de.de11 + e_de.de22 + e_de.de33;
		elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;
		e_de_vol *= one_third;
		e_de.de11 -= e_de_vol;
		e_de.de22 -= e_de_vol;
		e_de.de33 -= e_de_vol;

		const Displacement& n1_pdu = node_pdu[eni.n1];
		const Displacement& n2_pdu = node_pdu[eni.n2];
		const Displacement& n3_pdu = node_pdu[eni.n3];
		const Displacement& n4_pdu = node_pdu[eni.n4];
		e_pde_vol = (e_dN.dN1_dx * n1_pdu.ux + e_dN.dN2_dx * n2_pdu.ux + e_dN.dN3_dx * n3_pdu.ux + e_dN.dN4_dx * n4_pdu.ux)
				  + (e_dN.dN1_dy * n1_pdu.uy + e_dN.dN2_dy * n2_pdu.uy + e_dN.dN3_dy * n3_pdu.uy + e_dN.dN4_dy * n4_pdu.uy)
				  + (e_dN.dN1_dz * n1_pdu.uz + e_dN.dN2_dz * n2_pdu.uz + e_dN.dN3_dz * n3_pdu.uz + e_dN.dN4_dz * n4_pdu.uz);
		elem_m_pde_vol[e_id] = elem_pcl_m[e_id] * e_pde_vol;
	}

#pragma omp barrier

	double n_am_de_vol = 0.0;
	double n_am_pde_vol = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		e_id = node_elem_pair0[ve_id] / 4;
		assert(e_id < elem_num);
		n_am_de_vol += elem_m_de_vol[e_id];
		n_am_pde_vol += elem_m_pde_vol[e_id];
		if (n_id != node_has_elem0[ve_id + 1])
		{
			node_de_vol[n_id] = n_am_de_vol * one_fourth / node_am[n_id];
			node_pde_vol[n_id] = n_am_pde_vol * one_fourth / node_am[n_id];
			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < node_num || n_id == SIZE_MAX);
			n_am_de_vol = 0.0;
			n_am_pde_vol = 0.0;
		}
	}

#pragma omp barrier
	const Displacement* pn1_pdu, * pn2_pdu, * pn3_pdu, * pn4_pdu;
	double e_pde_vol_f;
	StrainInc* pe_de;
	const double* dstress, *estrain, * pstrain;
	e_id = SIZE_MAX;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elem0[p_id])
		{
			e_id = pcl_in_elem0[p_id];
			assert(e_id < elem_num);

			const ElemNodeIndex& eni = elem_node_id[e_id];
			pn1_pdu = node_pdu + eni.n1;
			pn2_pdu = node_pdu + eni.n2;
			pn3_pdu = node_pdu + eni.n3;
			pn4_pdu = node_pdu + eni.n4;

			e_de_vol = (node_de_vol[eni.n1] + node_de_vol[eni.n2]
					  + node_de_vol[eni.n3] + node_de_vol[eni.n4]) * one_fourth * one_third;
			pe_de = elem_de + e_id;
			pe_de->de11 += e_de_vol;
			pe_de->de22 += e_de_vol;
			pe_de->de33 += e_de_vol;

			e_pde_vol = (node_pde_vol[eni.n1] + node_pde_vol[eni.n2]
					   + node_pde_vol[eni.n3] + node_pde_vol[eni.n4]) * one_fourth;
			e_pde_vol_f = -e_pde_vol / elem_pcl_n[e_id];
			//elem_pcl_n[e_id] = (e_pde_vol + elem_pcl_n[e_id]) / (1.0 + e_pde_vol);
			//elem_pcl_int_vol[e_id] *= (1.0 + e_pde_vol);
			//elem_density_f[e_id] /= 1.0 - e_pde_vol_f;
			elem_p[e_id] += Kf * e_pde_vol_f;
		}

		// update displacement (for contact)
		ShapeFunc& p_N = pcl_N0[p_id];
		Displacement& p_u_s0 = pcl_u_s0[p_id];
		p_u_s0.ux += p_N.N1 * pn1_pdu->ux + p_N.N2 * pn2_pdu->ux + p_N.N3 * pn3_pdu->ux + p_N.N4 * pn4_pdu->ux;
		p_u_s0.uy += p_N.N1 * pn1_pdu->uy + p_N.N2 * pn2_pdu->uy + p_N.N3 * pn3_pdu->uy + p_N.N4 * pn4_pdu->uy;
		p_u_s0.uz += p_N.N1 * pn1_pdu->uz + p_N.N2 * pn2_pdu->uz + p_N.N3 * pn3_pdu->uz + p_N.N4 * pn4_pdu->uz;

		// update stress
		const size_t ori_p_id = pcl_index0[p_id];
		MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
		char* pcl_mm_copy = mat_model_copy + pcl_mat_model_copy_offset[ori_p_id];
		//t3d_chm_ud_mt_subit_db_file << ori_p_id << ", " << pcl_mat_model_copy_offset[ori_p_id] << "\n";
		pcl_mm.retrieve_from(pcl_mm_copy);
		pcl_mm.integrate(pe_de->de);
		dstress = pcl_mm.get_dstress();
		const Stress& p_ps = pcl_prev_stress[p_id];
		Stress& p_s0 = pcl_stress0[p_id];
		p_s0.s11 = p_ps.s11 + dstress[0];
		p_s0.s22 = p_ps.s22 + dstress[1];
		p_s0.s33 = p_ps.s33 + dstress[2];
		p_s0.s12 = p_ps.s12 + dstress[3];
		p_s0.s23 = p_ps.s23 + dstress[4];
		p_s0.s31 = p_ps.s31 + dstress[5];

		estrain = pcl_mm.get_dstrain_e();
		StrainInc& p_dee = pcl_destrain[p_id];
		p_dee.de11 = estrain[0];
		p_dee.de22 = estrain[1];
		p_dee.de33 = estrain[2];
		p_dee.de12 = estrain[3];
		p_dee.de23 = estrain[4];
		p_dee.de31 = estrain[5];

		pstrain = pcl_mm.get_dstrain_p();
		StrainInc& p_dpe = pcl_dpstrain[p_id];
		p_dpe.de11 = pstrain[0];
		p_dpe.de22 = pstrain[1];
		p_dpe.de33 = pstrain[2];
		p_dpe.de12 = pstrain[3];
		p_dpe.de23 = pstrain[4];
		p_dpe.de31 = pstrain[5];

		//if (ori_p_id == 200 && substep_index % 400 == 0 /* && cur_time > 0.30 && cur_time < 0.35 && */)
		//{
		//	Stress& ps = pcl_stress0[p_id];
		//	t3d_chm_ud_mt_subit_db_file << substep_index << ", "
		//		<< pe_de->de11 << ", " << pe_de->de22 << ", " << pe_de->de33 << ", "
		//		//<< pe_de->de12 << ", " << pe_de->de23 << ", " << pe_de->de31 << ", "
		//		<< pe_de->de11 + pe_de->de22 + pe_de->de33 << ", "
		//		//<< ps.s11 << ", " << ps.s22 << ", " << ps.s33 << ", "
		//		//<< ps.s12 << ", " << ps.s23 << ", " << ps.s31
		//		<< ",\n";
		//}
	}
	
	return 0;
}

int Step_T3D_CHM_ud_mt_subiter::apply_rigid_cylinder(
	size_t p_id0, size_t p_id1,
	const size_t* pcl_in_elem,
	const SortedPclVarArrays& cur_spva,
	Force3D& rc_cf,
	size_t substp_id,
	ThreadData& thd
	) noexcept
{
	double dist;
	Vector3D lnorm, gnorm;
	Point3D cur_cont_pos;
	Force lcont_fs, gcont_fs;
	Force lcont_ff, gcont_ff;
	const size_t* pcl_index = cur_spva.pcl_index;
	const Displacement* pcl_u_s = cur_spva.pcl_u_s;
	const ShapeFunc* pcl_N = cur_spva.pcl_N;
	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
	{
		const size_t ori_p_id = pcl_index[p_id];
		const Position& p_p = pcl_pos[ori_p_id];
		const Displacement& p_u = pcl_u_s[p_id];
		const double p_x = p_p.x + p_u.ux;
		const double p_y = p_p.y + p_u.uy;
		const double p_z = p_p.z + p_u.uz;
		const double p_r = 0.5 * pow(pcl_vol[p_id], one_third);
		const ShapeFunc& p_N = pcl_N[p_id];
		const size_t e_id = pcl_in_elem[p_id];
		if (prcy->detect_collision_with_point(
			p_x, p_y, p_z, p_r, dist, lnorm, cur_cont_pos))
		{
			prcy->get_global_vector(lnorm, gnorm);
			// solid pcl
			pcm_s->cal_contact_force(
				substp_id,
				ori_p_id,
				dist,
				lnorm,
				cur_cont_pos,
				p_r + p_r,
				pv_place_holder,
				lcont_fs.vec);
			prcy->get_global_vector(lcont_fs.vec, gcont_fs.vec);
			// fluid pcl
			pcm_f->cal_contact_force(
				substp_id,
				ori_p_id,
				dist,
				lnorm,
				cur_cont_pos,
				p_r + p_r,
				pv_place_holder,
				lcont_ff.vec);
			prcy->get_global_vector(lcont_ff.vec, gcont_ff.vec);
			// apply contact forc to bg mesh
			gcont_fs.fx += gcont_ff.fx;
			gcont_fs.fy += gcont_ff.fy;
			gcont_fs.fz += gcont_ff.fz;
			Force& en_fe1 = elem_node_f_ext[e_id * 4];
			en_fe1.fx += p_N.N1 * gcont_fs.fx;
			en_fe1.fy += p_N.N1 * gcont_fs.fy;
			en_fe1.fz += p_N.N1 * gcont_fs.fz;
			Force& en_fe2 = elem_node_f_ext[e_id * 4 + 1];
			en_fe2.fx += p_N.N2 * gcont_fs.fx;
			en_fe2.fy += p_N.N2 * gcont_fs.fy;
			en_fe2.fz += p_N.N2 * gcont_fs.fz;
			Force& en_fe3 = elem_node_f_ext[e_id * 4 + 2];
			en_fe3.fx += p_N.N3 * gcont_fs.fx;
			en_fe3.fy += p_N.N3 * gcont_fs.fy;
			en_fe3.fz += p_N.N3 * gcont_fs.fz;
			Force& en_fe4 = elem_node_f_ext[e_id * 4 + 3];
			en_fe4.fx += p_N.N4 * gcont_fs.fx;
			en_fe4.fy += p_N.N4 * gcont_fs.fy;
			en_fe4.fz += p_N.N4 * gcont_fs.fz;
			// apply contact force to rigid body
			const Point3D& rc_cen = prcy->get_centre();
			rc_cf.add_force(p_x, p_y, p_z,
				-gcont_fs.fx, -gcont_fs.fy, -gcont_fs.fz,
				 rc_cen.x, rc_cen.y, rc_cen.z);
		}
	}

	return 0;
}

int Step_T3D_CHM_ud_mt_subiter::apply_t3d_rigid_mesh(
	size_t p_id0, size_t p_id1,
	const size_t* pcl_in_elem,
	const SortedPclVarArrays& cur_spva,
	Force3D& rc_cf,
	size_t substp_id,
	ThreadData& thd
	) noexcept
{
	double dist;
	Vector3D lnorm, gnorm;
	Point3D cur_cont_pos;
	Force lcont_fs, gcont_fs;
	Force lcont_ff, gcont_ff;
	const size_t* pcl_index = cur_spva.pcl_index;
	const Displacement* pcl_u_s = cur_spva.pcl_u_s;
	const ShapeFunc* pcl_N = cur_spva.pcl_N;
	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
	{
		const size_t ori_p_id = pcl_index[p_id];
		const Position& p_p = pcl_pos[ori_p_id];
		const Displacement& p_u = pcl_u_s[p_id];
		const double p_x = p_p.x + p_u.ux;
		const double p_y = p_p.y + p_u.uy;
		const double p_z = p_p.z + p_u.uz;
		const double p_r = 0.5 * pow(pcl_vol[p_id], one_third);
		const ShapeFunc& p_N = pcl_N[p_id];
		const size_t e_id = pcl_in_elem[p_id];
		if (prm->detect_collision_with_point(
			p_x, p_y, p_z, p_r,
			dist, lnorm, cur_cont_pos))
		{
			prm->get_global_vector(lnorm, gnorm);
			// solid pcl
			pcm_s->cal_contact_force(
				substp_id,
				ori_p_id,
				dist,
				lnorm,
				cur_cont_pos,
				p_r + p_r,
				pv_place_holder,
				lcont_fs.vec);
			prm->get_global_vector(lcont_fs.vec, gcont_fs.vec);
			// fluid pcl
			pcm_f->cal_contact_force(
				substp_id,
				ori_p_id,
				dist,
				lnorm,
				cur_cont_pos,
				p_r + p_r,
				pv_place_holder,
				lcont_ff.vec);
			prm->get_global_vector(lcont_ff.vec, gcont_ff.vec);
			// map to bg mesh
			gcont_fs.fx += gcont_ff.fx;
			gcont_fs.fy += gcont_ff.fy;
			gcont_fs.fz += gcont_ff.fz;
			Force& en_fe1 = elem_node_f_ext[e_id * 4];
			en_fe1.fx += p_N.N1 * gcont_fs.fx;
			en_fe1.fy += p_N.N1 * gcont_fs.fy;
			en_fe1.fz += p_N.N1 * gcont_fs.fz;
			Force& en_fe2 = elem_node_f_ext[e_id * 4 + 1];
			en_fe2.fx += p_N.N2 * gcont_fs.fx;
			en_fe2.fy += p_N.N2 * gcont_fs.fy;
			en_fe2.fz += p_N.N2 * gcont_fs.fz;
			Force& en_fe3 = elem_node_f_ext[e_id * 4 + 2];
			en_fe3.fx += p_N.N3 * gcont_fs.fx;
			en_fe3.fy += p_N.N3 * gcont_fs.fy;
			en_fe3.fz += p_N.N3 * gcont_fs.fz;
			Force& en_fe4 = elem_node_f_ext[e_id * 4 + 3];
			en_fe4.fx += p_N.N4 * gcont_fs.fx;
			en_fe4.fy += p_N.N4 * gcont_fs.fy;
			en_fe4.fz += p_N.N4 * gcont_fs.fz;
			// apply contact force to rigid body
			const Point3D& rm_cen = prm->get_pos();
			rc_cf.add_force(p_x, p_y, p_z,
				-gcont_fs.fx, -gcont_fs.fy, -gcont_fs.fz,
				rm_cen.x, rm_cen.y, rm_cen.z);
		}
	}
	return 0;
}
