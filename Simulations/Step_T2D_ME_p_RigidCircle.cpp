#include "Simulations_pcp.h"

#include <cmath>
#include <iostream>

#include "Step_T2D_ME_p.h"

void Step_T2D_ME_p::cal_thread_func_RigidCircle(unsigned int th_id)
{
	step_barrier.wait();

	while (not_yet_completed)
	{
		init_cal_vars(th_id);
		cal_barrier.wait();
		
		find_pcls_in_which_elems(th_id);
		cal_barrier.wait();

		map_pcl_vars_to_nodes_at_elems(th_id);
		cal_barrier.wait();

		// rigid circle
		cal_contact_force(th_id);
		cal_barrier.wait();

		update_node_a_and_v(th_id);
		cal_barrier.wait();

		apply_a_and_v_bcs(th_id);
		cal_barrier.wait();

		cal_de_at_elem(th_id);
		cal_barrier.wait();

		map_de_vol_from_elem_to_node(th_id);
		cal_barrier.wait();

		update_pcl_vars(th_id);
		step_barrier.wait();
	}
}

int solve_substep_T2D_ME_p_RigidCircle(void *_self)
{
	typedef Model_T2D_ME_p::Node Node;
	typedef Model_T2D_ME_p::Element Element;
	typedef Model_T2D_ME_p::Particle Particle;
	typedef Step_T2D_ME_p::ThreadData ThreadData;
	typedef Step_T2D_ME_p::ThreadElemData ThreadElemData;
	
	Step_T2D_ME_p& self = *static_cast<Step_T2D_ME_p*>(_self);
	Model_T2D_ME_p& md = static_cast<Model_T2D_ME_p&>(self.get_model());
	RigidCircle& rc = md.rigid_circle;

	self.cur_node_id.store(0);
	self.cur_elem_id.store(0);
	self.step_barrier.lift_barrier();
	self.init_cal_vars(0);
	self.cal_barrier.wait_for_others();

	self.cur_pcl_id.store(0);
	self.cal_barrier.lift_barrier();
	self.find_pcls_in_which_elems(0);
	self.cal_barrier.wait_for_others();

	self.cur_elem_id.store(0);
	self.cal_barrier.lift_barrier();
	self.map_pcl_vars_to_nodes_at_elems(0);
	self.cal_barrier.wait_for_others();

	// rigid circle
	self.cal_contact_force(0);
	// accumulate contact force on rigid circle
	rc.reset_rf();
	for (size_t th_id = 0; th_id < self.thread_num; ++th_id)
	{
		ThreadData &th_data = *(self.pth_datas[th_id]);
		rc.add_rc_f(th_data.rcf);
	}
	// update motion
	rc.update_motion(self.dtime);
	self.cal_barrier.wait();

	self.cur_node_id.store(0);
	self.cal_barrier.lift_barrier();
	self.update_node_a_and_v(0);
	self.cal_barrier.wait_for_others();

	self.cur_ax_bc_id.store(0);
	self.cur_ay_bc_id.store(0);
	self.cur_vx_bc_id.store(0);
	self.cur_vy_bc_id.store(0);
	self.cal_barrier.lift_barrier();
	self.apply_a_and_v_bcs(0);
	self.cal_barrier.wait_for_others();

	self.cur_elem_id.store(0);
	self.cal_barrier.lift_barrier();
	self.cal_de_at_elem(0);
	self.cal_barrier.wait_for_others();

	self.cur_node_id.store(0);
	self.cal_barrier.lift_barrier();
	self.map_de_vol_from_elem_to_node(0);
	self.cal_barrier.wait_for_others();

	self.cur_pcl_id.store(0);
	self.cal_barrier.lift_barrier();
	self.update_pcl_vars(0);

	self.step_barrier.wait_for_others();

	return 0;
}

void Step_T2D_ME_p::cal_contact_force(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);
	RigidCircle& rc = md.rigid_circle;
	RigidCircleForce &rcf = pth_datas[th_id]->rcf;

	double rc_x = rc.get_x(), rc_y = rc.get_y();
	rcf.reset_rf();
	double dist, norm_x, norm_y, f_cont, fx_cont, fy_cont;
	for (size_t e_id = cur_elem_id++; e_id < md.elem_num; e_id = cur_elem_id++)
	{
		Element& e = md.elems[e_id];

		if (!e.has_mp /* || not contact circle */)
			continue;
	
		NodeVarAtElem& nv1 = e.node_vars[0];
		NodeVarAtElem& nv2 = e.node_vars[1];
		NodeVarAtElem& nv3 = e.node_vars[2];

		for (Particle *ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
			 ppcl = e.next_pcl(ppcl))
		{
			Particle& pcl = *ppcl;
			if (rc.detect_collision_with_point(pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
			{
				f_cont = md.K_cont * dist;
				fx_cont = f_cont * norm_x;
				fy_cont = f_cont * norm_y;
				// add reaction force to rigid object
				rcf.add_rf(pcl.x, pcl.y, -fx_cont, -fy_cont, rc_x, rc_y);
				// node 1
				nv1.fx += pcl.N1 * fx_cont;
				nv1.fy += pcl.N1 * fy_cont;
				// node 2
				nv2.fx += pcl.N2 * fx_cont;
				nv2.fy += pcl.N2 * fy_cont;
				// node 3
				nv3.fx += pcl.N3 * fx_cont;
				nv3.fy += pcl.N3 * fy_cont;
			}
		}
	}
}
