#include "Simulations_pcp.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/tick_count.h>

#include "Step_T2D_ME_p_tbb.h"

Step_T2D_ME_p_tbb::Step_T2D_ME_p_tbb(const char* _name) : 
	Step_T2D_ME_p(_name),
	thread_pool_config(tbb::task_scheduler_init::deferred)
{
	cal_substep_func = &solve_substep_T2D_ME_p_tbb;
}

Step_T2D_ME_p_tbb::~Step_T2D_ME_p_tbb()
{
	
}

int Step_T2D_ME_p_tbb::init_calculation()
{
	Model_T2D_ME_p& md = *model;

	if (is_first_step)
	{
		md.init_bcs();
	}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle& pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol = pcl.m / pcl.density;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}

	if (thread_num == 0)
		thread_num = std::thread::hardware_concurrency();
	thread_pool_config.initialize(int(thread_num));

	plist_mem.init_mem(md.elem_num, size_t(thread_num)*2);

	init_node_time = 0.0;
	find_pcl_in_elem_time = 0.0;
	map_pcl_to_node_time = 0.0;
	rigid_circle_time = 0.0;
	update_node_a_v_time = 0.0;
	apply_a_v_bc_time = 0.0;
	cal_strain_at_elem_time = 0.0;
	map_de_vol_from_elem_to_node_time = 0.0;
	update_pcl_time = 0.0;
	total_parallel_time = 0.0;

	return 0;
}

int Step_T2D_ME_p_tbb::finalize_calculation()
{
	thread_pool_config.terminate();

	total_parallel_time = 
		  init_node_time
		+ find_pcl_in_elem_time
		+ map_pcl_to_node_time
		+ rigid_circle_time
		+ update_node_a_v_time
		+ apply_a_v_bc_time
		+ cal_strain_at_elem_time
		+ map_de_vol_from_elem_to_node_time
		+ update_pcl_time;
	
	std::cout << "\ninit_node, " << init_node_time << "\n"
		"find_pcl_int_elem, " << find_pcl_in_elem_time << "\n"
		"map_pcl_to_node, " << map_pcl_to_node_time << "\n"
		"rigid_circle, " << rigid_circle_time << "\n"
		"update_node_a_v, " << update_node_a_v_time << "\n"
		"apply_a_v_bc, " << apply_a_v_bc_time << "\n"
		"cal_strain_at_elem, " << cal_strain_at_elem_time << "\n"
		"map_de_vol_from_elem_to_node, " << map_de_vol_from_elem_to_node_time << "\n"
		"update_pcl, " << update_pcl_time << "\n"
		"total_parallel, " << total_parallel_time;

	return 0;
}

namespace
{
	typedef Step_T2D_ME_p_tbb::NodeToElem NodeToElem;
	typedef Step_T2D_ME_p_tbb::Node Node;
	typedef Step_T2D_ME_p_tbb::NodeVarAtElem NodeVarAtElem;
	typedef Step_T2D_ME_p_tbb::Element Element;
	typedef Step_T2D_ME_p_tbb::Particle Particle;
	typedef Step_T2D_ME_p_tbb::PclListAtElem PclListAtElem;

	using namespace tbb;

	class InitNodeOp
	{
	protected:
		Node* nodes;

	public:
		InitNodeOp(Node* _nodes) : nodes(_nodes) {}

		void operator()(const blocked_range<size_t>& br) const
		{
			for (size_t n_id = br.begin(); n_id < br.end(); ++n_id)
			{
				Node& n = nodes[n_id];
				n.has_mp = false;
			}
		}
	};

	class FindPclInWhichElemOp
	{
	protected:
		Model_T2D_ME_p& md;
		size_t elem_num;
		Element* elems;
		Particle* pcls;

		FixedSizeMemAllocator<PclListAtElem> &plist_mem;
		PclListAtElem* elem_plists;

		inline void init_teds()
		{
			elem_plists = plist_mem.alloc_mem();
			for (size_t e_id = 0; e_id < elem_num; ++e_id)
				elem_plists[e_id].init();
		}

	public:
		FindPclInWhichElemOp(
			Model_T2D_ME_p& _md,
			FixedSizeMemAllocator<PclListAtElem> &_plist_mem
			) :
			md(_md),
			elems(_md.get_elems()), elem_num(_md.get_elem_num()),
			pcls(md.get_pcls()),
			plist_mem(_plist_mem)
		{
			init_teds();
		}

		FindPclInWhichElemOp(const FindPclInWhichElemOp& other, split) :
			md(other.md),
			elem_num(other.elem_num), elems(other.elems),
			pcls(other.pcls),
			plist_mem(other.plist_mem)
		{
			init_teds();
		}

		void operator() (const blocked_range<size_t>& br)
		{
			for (size_t pcl_id = br.begin(); pcl_id < br.end(); ++pcl_id)
			{
				Particle& pcl = pcls[pcl_id];

				if (pcl.pe == nullptr)
					continue;
				pcl.pe = md.find_in_which_element(pcl);
				if (pcl.pe == nullptr)
					continue;

				elem_plists[pcl.pe->id].add_pcl(pcl);
			}
		}

		void join(const FindPclInWhichElemOp& other)
		{
			for (size_t e_id = 0; e_id < elem_num; ++e_id)
				elem_plists[e_id].combine(other.elem_plists[e_id]);
		}

		inline PclListAtElem *get_plists() noexcept { return elem_plists; }
	};

	class MapPclToNodeOp
	{
	protected:
		Node* nodes;
		Element* elems;
		Particle* pcls;
		PclListAtElem *elem_plists;

	public:
		MapPclToNodeOp(Model_T2D_ME_p& md,
			PclListAtElem* _elem_plists) :
			nodes(md.get_nodes()),
			elems(md.get_elems()),
			pcls(md.get_pcls()),
			elem_plists(_elem_plists) {}

		void operator() (const blocked_range<size_t>& br)  const
		{
			for (size_t e_id = br.begin(); e_id < br.end(); ++e_id)
			{
				Element &e = elems[e_id];

				e.pcls = elem_plists[e_id].get_top();
				if (!e.pcls)
					continue;

				e.s11 = 0.0;
				e.s22 = 0.0;
				e.s12 = 0.0;
				e.mi_pcl_vol = 0.0;

				Node& n1 = nodes[e.n1];
				Node& n2 = nodes[e.n2];
				Node& n3 = nodes[e.n3];
				n1.has_mp = true;
				n2.has_mp = true;
				n3.has_mp = true;

				NodeVarAtElem& nv1 = e.node_vars[0];
				NodeVarAtElem& nv2 = e.node_vars[1];
				NodeVarAtElem& nv3 = e.node_vars[2];
				nv1.init();
				nv2.init();
				nv3.init();

				for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
					 ppcl = e.next_pcl(ppcl))
				{
					Particle& pcl = *ppcl;

					// mixed integration
					pcl.vol = pcl.m / pcl.density;
					e.s11 += pcl.s11 * pcl.vol;
					e.s22 += pcl.s22 * pcl.vol;
					e.s12 += pcl.s12 * pcl.vol;
					e.mi_pcl_vol += pcl.vol;

					// map velocity
					double mvx = pcl.m * pcl.vx;
					double mvy = pcl.m * pcl.vy;
					// node 1
					nv1.m += pcl.N1 * pcl.m;
					nv1.vx += pcl.N1 * mvx;
					nv1.vy += pcl.N1 * mvy;
					// node 2
					nv2.m += pcl.N2 * pcl.m;
					nv2.vx += pcl.N2 * mvx;
					nv2.vy += pcl.N2 * mvy;
					// node 3
					nv3.m += pcl.N3 * pcl.m;
					nv3.vx += pcl.N3 * mvx;
					nv3.vy += pcl.N3 * mvy;

					// external load
					if (pcl.has_fx_ext)
					{
						nv1.fx += pcl.N1 * pcl.fx_ext;
						nv2.fx += pcl.N2 * pcl.fx_ext;
						nv3.fx += pcl.N3 * pcl.fx_ext;
					}
					if (pcl.has_fy_ext)
					{
						nv1.fy += pcl.N1 * pcl.fy_ext;
						nv2.fy += pcl.N2 * pcl.fy_ext;
						nv3.fy += pcl.N3 * pcl.fy_ext;
					}
				}

				e.s11 /= e.mi_pcl_vol;
				e.s22 /= e.mi_pcl_vol;
				e.s12 /= e.mi_pcl_vol;
				if (e.mi_pcl_vol > e.area)
					e.mi_pcl_vol = e.area;

				// internal force
				// node 1
				nv1.fx -= (e.dN1_dx * e.s11 + e.dN1_dy * e.s12) * e.mi_pcl_vol;
				nv1.fy -= (e.dN1_dx * e.s12 + e.dN1_dy * e.s22) * e.mi_pcl_vol;
				// node 2
				nv2.fx -= (e.dN2_dx * e.s11 + e.dN2_dy * e.s12) * e.mi_pcl_vol;
				nv2.fy -= (e.dN2_dx * e.s12 + e.dN2_dy * e.s22) * e.mi_pcl_vol;
				// node 3
				nv3.fx -= (e.dN3_dx * e.s11 + e.dN3_dy * e.s12) * e.mi_pcl_vol;
				nv3.fy -= (e.dN3_dx * e.s12 + e.dN3_dy * e.s22) * e.mi_pcl_vol;
			}
		}
	};

	class UpdateNodeAAndVOp
	{
	protected:
		Node* nodes;
		Element* elems;
		double dtime;

	public:
		UpdateNodeAAndVOp(Model_T2D_ME_p& md, double dt) :
			nodes(md.get_nodes()),
			elems(md.get_elems()),
			dtime(dt) {}

		void operator() (const blocked_range<size_t>& br) const
		{
			for (size_t n_id = br.begin(); n_id < br.end(); ++n_id)
			{
				Node& n = nodes[n_id];

				if (!n.has_mp)
					continue;

				n.m = 0.0;
				n.vx_ori = 0.0;
				n.vy_ori = 0.0;
				n.fx = 0.0;
				n.fy = 0.0;

				for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
				{
					NodeToElem& n2e = n.n2es[e_id];
					Element& e = elems[n2e.e_id];
					if (e.pcls)
					{
						NodeVarAtElem& nv = e.node_vars[n2e.n_id];
						n.m += nv.m;
						n.vx_ori += nv.vx;
						n.vy_ori += nv.vy;
						n.fx += nv.fx;
						n.fy += nv.fy;
					}
				}

				n.ax = n.fx / n.m;
				n.ay = n.fy / n.m;
				n.vx_ori /= n.m;
				n.vy_ori /= n.m;
				n.vx = n.vx_ori + n.ax * dtime;
				n.vy = n.vy_ori + n.ay * dtime;
				n.dux = n.vx * dtime;
				n.duy = n.vy * dtime;
			}
		}
	};

	class ApplyAxBcOp
	{
	protected:
		Node* nodes;
		AccelerationBC* abcs;
	public:
		ApplyAxBcOp(Model_T2D_ME_p &md) :
			nodes(md.get_nodes()),
			abcs(md.get_axs()) {}

		void operator() (const blocked_range<size_t>& br) const
		{
			for (size_t a_id = br.begin(); a_id < br.end(); ++a_id)
			{
				AccelerationBC& abc = abcs[a_id];
				nodes[abc.node_id].ax = abc.a;
			}
		}
	};

	class ApplyAyBcOp
	{
	protected:
		Node* nodes;
		AccelerationBC* abcs;
	public:
		ApplyAyBcOp(Model_T2D_ME_p& md) :
			nodes(md.get_nodes()),
			abcs(md.get_ays()) {}

		void operator() (const blocked_range<size_t>& br) const
		{
			for (size_t a_id = br.begin(); a_id < br.end(); ++a_id)
			{
				AccelerationBC& abc = abcs[a_id];
				nodes[abc.node_id].ay = abc.a;
			}
		}
	};

	class ApplyVxBcOp
	{
	protected:
		Node* nodes;
		VelocityBC* vbcs;
	public:
		ApplyVxBcOp(Model_T2D_ME_p& md) :
			nodes(md.get_nodes()),
			vbcs(md.get_vxs()) {}

		void operator() (const blocked_range<size_t>& br) const
		{
			for (size_t v_id = br.begin(); v_id < br.end(); ++v_id)
			{
				VelocityBC& vbc = vbcs[v_id];
				Node& n = nodes[vbc.node_id];
				n.ax = 0.0;
				n.vx = vbc.v;
			}
		}
	};

	class ApplyVyBcOp
	{
	protected:
		Node* nodes;
		VelocityBC* vbcs;
	public:
		ApplyVyBcOp(Model_T2D_ME_p& md) :
			nodes(md.get_nodes()),
			vbcs(md.get_vys()) {}

		void operator() (const blocked_range<size_t>& br) const
		{
			for (size_t v_id = br.begin(); v_id < br.end(); ++v_id)
			{
				VelocityBC& vbc = vbcs[v_id];
				Node& n = nodes[vbc.node_id];
				n.ay = 0.0;
				n.vy = vbc.v;
			}
		}
	};

	class CalElemStrainOp
	{
	protected:
		Node* nodes;
		Element* elems;

	public:
		CalElemStrainOp(Model_T2D_ME_p& md) :
			nodes(md.get_nodes()),
			elems(md.get_elems()) {}

		void operator() (const blocked_range<size_t>& br) const
		{
			double de11, de22;
			for (size_t e_id = br.begin(); e_id < br.end(); ++e_id)
			{
				Element& e = elems[e_id];

				if (!e.pcls)
					continue;

				Node& n1 = nodes[e.n1];
				Node& n2 = nodes[e.n2];
				Node& n3 = nodes[e.n3];
				de11 = n1.dux * e.dN1_dx + n2.dux * e.dN2_dx + n3.dux * e.dN3_dx;
				de22 = n1.duy * e.dN1_dy + n2.duy * e.dN2_dy + n3.duy * e.dN3_dy;
				e.de_vol_by_3 = (de11 + de22) / 3.0;
				e.dde11 = de11 - e.de_vol_by_3;
				e.dde22 = de22 - e.de_vol_by_3;
				e.de12 = (n1.dux * e.dN1_dy + n2.dux * e.dN2_dy + n3.dux * e.dN3_dy
						+ n1.duy * e.dN1_dx + n2.duy * e.dN2_dx + n3.duy * e.dN3_dx) * 0.5;

				// strain enhancement
				NodeVarAtElem& nv1 = e.node_vars[0];
				NodeVarAtElem& nv2 = e.node_vars[1];
				NodeVarAtElem& nv3 = e.node_vars[2];
				double N_vol;
				for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
					ppcl = e.next_pcl(ppcl))
				{
					Particle& pcl = *ppcl;
					// node 1
					N_vol = pcl.N1 * pcl.vol;
					nv1.de_vol_by_3 += e.de_vol_by_3 * N_vol;
					nv1.se_pcl_vol += N_vol;
					// node 2
					N_vol = pcl.N2 * pcl.vol;
					nv2.de_vol_by_3 += e.de_vol_by_3 * N_vol;
					nv2.se_pcl_vol += N_vol;
					// node 3
					N_vol = pcl.N3 * pcl.vol;
					nv3.de_vol_by_3 += e.de_vol_by_3 * N_vol;
					nv3.se_pcl_vol += N_vol;
				}
			}
		}
	};

	class MapVolStrainFromElemToNodeOp
	{
	protected:
		Node* nodes;
		Element* elems;

	public:
		MapVolStrainFromElemToNodeOp(Model_T2D_ME_p &md) : 
			nodes(md.get_nodes()),
			elems(md.get_elems()) {}

		void operator() (const blocked_range<size_t>& br) const
		{
			for (size_t n_id = br.begin(); n_id < br.end(); ++n_id)
			{
				Node& n = nodes[n_id];

				if (!n.has_mp)
					continue;

				n.de_vol_by_3 = 0.0;
				n.se_pcl_vol = 0.0;
				for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
				{
					NodeToElem& n2e = n.n2es[e_id];
					Element& e = elems[n2e.e_id];
					if (e.pcls)
					{
						NodeVarAtElem& nv = e.node_vars[n2e.n_id];
						n.de_vol_by_3 += nv.de_vol_by_3;
						n.se_pcl_vol += nv.se_pcl_vol;
					}
				}
				n.de_vol_by_3 /= n.se_pcl_vol;
			}
		}
	};

	class UpdatePclVarOp
	{
	protected:
		Node *nodes;
		Particle *pcls;
		double dtime;

	public:
		UpdatePclVarOp(Model_T2D_ME_p &md, double dt) :
			nodes(md.get_nodes()),
			pcls(md.get_pcls()),
			dtime(dt) {}

		void operator() (const blocked_range<size_t> &br) const
		{
			double de_vol_by_3, de11, de22, de12;
			for (size_t pcl_id = br.begin(); pcl_id < br.end(); ++pcl_id)
			{
				Particle& pcl = pcls[pcl_id];
				if (!pcl.pe)
					return;

				Element& e = *pcl.pe;
				Node& n1 = nodes[e.n1];
				Node& n2 = nodes[e.n2];
				Node& n3 = nodes[e.n3];

				// velocity
				pcl.vx += (n1.ax * pcl.N1 + n2.ax * pcl.N2 + n3.ax * pcl.N3) * dtime;
				pcl.vy += (n1.ay * pcl.N1 + n2.ay * pcl.N2 + n3.ay * pcl.N3) * dtime;

				// displacement
				pcl.ux += n1.dux * pcl.N1 + n2.dux * pcl.N2 + n3.dux * pcl.N3;
				pcl.uy += n1.duy * pcl.N1 + n2.duy * pcl.N2 + n3.duy * pcl.N3;

				// position
				pcl.x = pcl.x_ori + pcl.ux;
				pcl.y = pcl.y_ori + pcl.uy;

				// strain
				de_vol_by_3 = (n1.de_vol_by_3 + n2.de_vol_by_3 + n3.de_vol_by_3) / 3.0;
				de11 = e.dde11 + de_vol_by_3;
				de22 = e.dde22 + de_vol_by_3;
				de12 = e.de12;
				pcl.e11 += de11;
				pcl.e22 += de22;
				pcl.e12 += de12;

				// stress
				double dstrain[6] = { de11, de22, 0.0, de12, 0.0, 0.0 };
				pcl.mm->integrate(dstrain);
				const double* dstress = pcl.mm->get_dstress();
				pcl.s11 += dstress[0];
				pcl.s22 += dstress[1];
				pcl.s12 += dstress[3];

				// density
				pcl.density /= 1.0 + de_vol_by_3 * 3.0;
			}
		}
	};

	// rigid circle
	class RigidCircleContactForceOp
	{
	protected:
		Model_T2D_ME_p &md;
		Element* elems;
		Node* nodes;

		RigidCircle &rc;
		double K_cont;

		RigidCircleForce cf;

	public:
		RigidCircleContactForceOp(Model_T2D_ME_p &_md) :
			md(_md), elems(md.get_elems()), nodes(md.get_nodes()),
			rc(md.get_rigid_circle()), K_cont(md.get_K_cont()) {}

		RigidCircleContactForceOp(const RigidCircleContactForceOp &other, split):
			md(other.md), elems(other.elems), nodes(other.nodes),
			rc(other.rc), K_cont(other.K_cont) {}

		void join(const RigidCircleContactForceOp &other)
		{
			cf.combine(other.cf);
		}

		void operator()(const blocked_range<size_t> &rb)
		{
			cf.reset_rf();
			double dist, norm_x, norm_y, f_cont, fx_cont, fy_cont;
			double rc_x = rc.get_x(), rc_y = rc.get_y();
			for (size_t e_id = rb.begin(); e_id < rb.end(); ++e_id)
			{
				Element& e = elems[e_id];

				if (!e.pcls /* || not contact circle */)
					continue;

				NodeVarAtElem& nv1 = e.node_vars[0];
				NodeVarAtElem& nv2 = e.node_vars[1];
				NodeVarAtElem& nv3 = e.node_vars[2];

				for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
					ppcl = e.next_pcl(ppcl))
				{
					Particle& pcl = *ppcl;
					if (rc.detect_collision_with_point(pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
					{
						f_cont = K_cont * dist;
						fx_cont = f_cont * norm_x;
						fy_cont = f_cont * norm_y;
						// add reaction force to rigid object
						cf.add_rf(pcl.x, pcl.y, -fx_cont, -fy_cont, rc_x, rc_y);
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

		inline RigidCircleForce &get_cf() noexcept { return cf; }
	};
}

int solve_substep_T2D_ME_p_tbb(void* _self)
{
	Step_T2D_ME_p_tbb& self = *static_cast<Step_T2D_ME_p_tbb*>(_self);
	Model_T2D_ME_p& md = static_cast<Model_T2D_ME_p&>(self.get_model());
	
	blocked_range<size_t> node_range(0, md.get_node_num());
	blocked_range<size_t> elem_range(0, md.get_elem_num());
	blocked_range<size_t> pcl_range(0, md.get_pcl_num());
	blocked_range<size_t> ax_range(0, md.get_ax_num());
	blocked_range<size_t> ay_range(0, md.get_ay_num());
	blocked_range<size_t> vx_range(0, md.get_vx_num());
	blocked_range<size_t> vy_range(0, md.get_vy_num());

	tick_count t0, t1;

	// init node
	t0 = tick_count::now();
	parallel_for(node_range, InitNodeOp(md.get_nodes()));
	t1 = tick_count::now();
	self.init_node_time += (t1 - t0).seconds();

	// find pcl in which element
	self.plist_mem.reset();
	FindPclInWhichElemOp find_pcl_in_which_elem_op(md, self.plist_mem);
	t0 = tick_count::now();
	parallel_reduce(pcl_range, find_pcl_in_which_elem_op);
	t1 = tick_count::now();
	self.find_pcl_in_elem_time += (t1 - t0).seconds();

	// map pcl vars to nodes
	t0 = tick_count::now();
	parallel_for(elem_range, MapPclToNodeOp(md, find_pcl_in_which_elem_op.get_plists()));
	t1 = tick_count::now();
	self.map_pcl_to_node_time += (t1 - t0).seconds();

	// rigid circle contact
	if (md.rigid_circle_is_valid())
	{
		RigidCircleContactForceOp rc_cf_op(md);
		t0 = tick_count::now();
		parallel_reduce(elem_range, rc_cf_op);
		t1 = tick_count::now();
		self.rigid_circle_time += (t1 - t0).seconds();
		md.rigid_circle.set_rcf(rc_cf_op.get_cf());
		md.rigid_circle.update_motion(self.dtime);
	}
	
	// update node a and v
	t0 = tick_count::now();
	parallel_for(node_range, UpdateNodeAAndVOp(md, self.dtime));
	t1 = tick_count::now();
	self.update_node_a_v_time += (t1 - t0).seconds();

	// apply a anv v bcs
	t0 = tick_count::now();
	parallel_for(ax_range, ApplyAxBcOp(md));
	parallel_for(ay_range, ApplyAyBcOp(md));
	parallel_for(vx_range, ApplyVxBcOp(md));
	parallel_for(vy_range, ApplyVyBcOp(md));
	t1 = tick_count::now();
	self.apply_a_v_bc_time += (t1 - t0).seconds();

	// cal strain at elems
	t0 = tick_count::now();
	parallel_for(elem_range, CalElemStrainOp(md));
	t1 = tick_count::now();
	self.cal_strain_at_elem_time += (t1 - t0).seconds();

	// map de_vol from elem to node
	t0 = tick_count::now();
	parallel_for(node_range, MapVolStrainFromElemToNodeOp(md));
	t1 = tick_count::now();
	self.map_de_vol_from_elem_to_node_time += (t1 - t0).seconds();

	// update pcl vars
	t0 = tick_count::now();
	parallel_for(pcl_range, UpdatePclVarOp(md, self.dtime));
	t1 = tick_count::now();
	self.update_pcl_time += (t1 - t0).seconds();

	return 0;
}
