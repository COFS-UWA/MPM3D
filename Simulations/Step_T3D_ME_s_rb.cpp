#include "Simulations_pcp.h"

#include <iostream>
#include <assert.h>
#include <cmath>
#include "MaterialModel.h"

#include "Step_T3D_ME_s.h"

int Step_T3D_ME_s::apply_rb_to_mesh(RigidTetrahedronMesh& rb)
{
	Model_T3D_ME_s& md = *model;
	SearchingGrid grid = md.get_bg_grid();
	const IdCube& grid_id_box = md.get_bg_grid_id_box();

	const Cube& rb_box = rb.get_grid_bbox();
	IdCube rb_id_box;
	rb_id_box.from_cube(
		rb_box,
		grid.get_xl(),
		grid.get_yl(),
		grid.get_zl(),
		grid.get_hx(),
		grid.get_hy(),
		grid.get_hz()
		);
	if (rb_id_box.does_not_overlap(grid_id_box))
		return 0;
	rb_id_box.trim_by(grid_id_box);

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
					for (Particle* ppcl = g.data.pcls; ppcl; ppcl = ppcl->next_in_grid)
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

	return 0;
}

bool Step_T3D_ME_s::apply_pcl_contact_force(Particle& pcl, RigidTetrahedronMesh& rb)
{
	Model_T3D_ME_s& md = *model;
	double dist, nx, ny, nz;
	double fn_cont, fnx_cont, fny_cont, fnz_cont;
	double ft_cont, ftx_cont, fty_cont, ftz_cont;
	double fx_cont, fy_cont, fz_cont;
	Point3D pt(pcl.x, pcl.y, pcl.z);
	if (rb.cal_distance_to_boundary(pt, dist, nx, ny, nz))
	{
		dist += pow(pcl.vol, 0.333333333);
		if (dist > 0.0)
		{
			fn_cont = md.K_cont * dist;
			fnx_cont = nx * fn_cont;
			fny_cont = ny * fn_cont;
			fnz_cont = nz * fn_cont;
			// friction force here
			ftx_cont = 0.0;
			fty_cont = 0.0;
			ftz_cont = 0.0;
			fx_cont = fnx_cont + ftx_cont;
			fy_cont = fny_cont + fty_cont;
			fz_cont = fnz_cont + ftz_cont;
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
				fx_cont,
				fy_cont,
				fz_cont,
				0.0,
				0.0,
				0.0
				);
			return true;
		}
	}
	return false;
}
