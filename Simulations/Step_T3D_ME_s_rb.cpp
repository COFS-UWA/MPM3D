#include "Simulations_pcp.h"

#include <iostream>
#include <assert.h>
#include <cmath>
#include "MaterialModel.h"

#include "Step_T3D_ME_s.h"

int Step_T3D_ME_s::apply_rb_to_mesh(RigidTetrahedronMesh& rb)
{
	Model_T3D_ME_s &md = *model;
	SearchingGrid &grid = md.get_bg_grid();
	const IdCube &grid_id_box = md.get_bg_grid_id_box();

	const Cube rb_box = rb.get_cur_bbox();
	IdCube rb_id_box;
	rb_id_box.from_cube(rb_box,
		grid.get_xl(), grid.get_yl(), grid.get_zl(),
		grid.get_hx(), grid.get_hy(), grid.get_hz()
		);
	if (rb_id_box.does_not_overlap(grid_id_box))
		return 0;
	rb_id_box.trim_by(grid_id_box);

	contact_state_list.reset_all_contact_state();
	SearchingGrid::Grid* cur_grid, *grid_y_start, *grid_z_start;
	grid_z_start = &grid.get_grid(rb_id_box.xl_id, rb_id_box.yl_id, rb_id_box.zl_id);
	for (long long z_id = rb_id_box.zl_id; z_id < rb_id_box.zu_id; ++z_id)
	{
		grid_y_start = grid_z_start;
		for (long long y_id = rb_id_box.yl_id; y_id < rb_id_box.yu_id; ++y_id)
		{
			cur_grid = grid_y_start;
			for (long long x_id = rb_id_box.xl_id; x_id < rb_id_box.xu_id; ++x_id)
			{
				SearchingGrid::Grid &g = *cur_grid;
				if (g.data.pcls /*&& detect_obb_cube_collision()*/)
				{
					for (Particle* ppcl = g.data.pcls; ppcl;
						 ppcl = ppcl->next_in_grid)
					{
						if (ppcl->pe)
							apply_pcl_contact_force(*ppcl, rb);
					}
				}
				++cur_grid;
			}
			grid_y_start += grid.get_x_num();
		}
		grid_z_start += grid.get_xy_num();
	}
	contact_state_list.clear_not_in_contact();

	return 0;
}

bool Step_T3D_ME_s::apply_pcl_contact_force(Particle& pcl, RigidTetrahedronMesh& rb)
{
	Model_T3D_ME_s& md = *model;
	double dist, nx, ny, nz;
	double fn_cont, fnx_cont, fny_cont, fnz_cont;
	double dtx, dty, dtz, dtn, dt;
	double ft_cont, ft_cont_max, ft_ratio;
	double ftx_cont, fty_cont, ftz_cont;
	double fx_cont, fy_cont, fz_cont;
	if (rb.cal_dist_and_dir_to_pt(pcl, dist, nx, ny, nz))
	{
		dist += pow(pcl.vol, 0.333333333) * 0.5;
		if (dist > 0.0)
		{
			// Normal force
			fn_cont = md.Kn_cont * dist;
			fnx_cont = nx * fn_cont;
			fny_cont = ny * fn_cont;
			fnz_cont = nz * fn_cont;
			// Frictional force
			// update contact state
			ContactState& cs = pcl.contact_state;
			if (cs.prev_in_contact)
				// already in contact state list
				cs.cur_in_contact = true;
			else // not yet in contact state list
				contact_state_list.add_pcl(pcl);
			// Tangential relative movement
			dtx = (pcl.vx - rb.get_vx()) * dtime; // currently no rotation
			dty = (pcl.vy - rb.get_vy()) * dtime;
			dtz = (pcl.vz - rb.get_vz()) * dtime;
			dtn = dtx * nx + dty * ny + dtz * nz;
			dtx -= dtn * nx;
			dty -= dtn * ny;
			dtz -= dtn * nz;
			dt = sqrt(dtx * dtx + dty * dty + dtz * dtz);
			// Tangential contact force model: ft = next_ft (cur_ft, ds, ds_dt)
			// 1. ft = 0.0 (smooth contact)
			// do nothing
			// 2. ft = miu * fn (frictional contact)
			ftx_cont = 0.0;
			fty_cont = 0.0;
			ftz_cont = 0.0;
			if (dt > 1.0e-10)
			{
				// adjust ft direction
				Vector3D& ft_dir = cs.ft_dir;
				double ftn = ft_dir.dot(nx, ny, nz);
				ft_dir.x -= ftn * nx;
				ft_dir.y -= ftn * ny;
				ft_dir.z -= ftn * nz;
				ft_dir.normalize();
				// update ft
				ft_cont = md.Kt_cont * dt;
				ftx_cont = cs.ft_cont * ft_dir.x + ft_cont * dtx / dt;
				fty_cont = cs.ft_cont * ft_dir.y + ft_cont * dty / dt;
				ftz_cont = cs.ft_cont * ft_dir.z + ft_cont * dtz / dt;
				ft_cont = sqrt(ftx_cont * ftx_cont + fty_cont * fty_cont + ftz_cont * ftz_cont);
				ft_cont_max = md.miu_cont * fn_cont;
				if (ft_cont_max < ft_cont)
				{
					if (ft_cont != 0.0)
					{
						ft_ratio = ft_cont_max / ft_cont;
						ftx_cont *= ft_ratio;
						fty_cont *= ft_ratio;
						ftz_cont *= ft_ratio;
					}
					else
					{
						ftx_cont = 0.0;
						fty_cont = 0.0;
						ftz_cont = 0.0;
					}
				}
				cs.ft_cont = ft_cont;
				ft_dir.x = ftx_cont;
				ft_dir.y = fty_cont;
				ft_dir.z = ftz_cont;
				ft_dir.normalize();
			}

			fx_cont = fnx_cont + ftx_cont;
			fy_cont = fny_cont + fty_cont;
			fz_cont = fnz_cont + ftz_cont;
			// apply to bg grid
			Element &e = *pcl.pe;
			Node& n1 = md.nodes[e.n1];
			n1.fx_ext += pcl.N1 * fx_cont;
			n1.fy_ext += pcl.N1 * fy_cont;
			n1.fz_ext += pcl.N1 * fz_cont;
			Node& n2 = md.nodes[e.n2];
			n2.fx_ext += pcl.N2 * fx_cont;
			n2.fy_ext += pcl.N2 * fy_cont;
			n2.fz_ext += pcl.N2 * fz_cont;
			Node& n3 = md.nodes[e.n3];
			n3.fx_ext += pcl.N3 * fx_cont;
			n3.fy_ext += pcl.N3 * fy_cont;
			n3.fz_ext += pcl.N3 * fz_cont;
			Node& n4 = md.nodes[e.n4];
			n4.fx_ext += pcl.N4 * fx_cont;
			n4.fy_ext += pcl.N4 * fy_cont;
			n4.fz_ext += pcl.N4 * fz_cont;
			// currently not consider rotation
			rb.add_con_force(
				-fx_cont,
				-fy_cont,
				-fz_cont,
				0.0,
				0.0,
				0.0
				);
			return true;
		}
	}
	return false;
}
