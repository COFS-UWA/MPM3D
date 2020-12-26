#ifndef __Rigid_Object_Rigid_Circle_h__
#define __Rigid_Object_Rigid_Circle_h__

#include <cmath>

#include "Geometry2D.h"
#include "Force2D.h"

namespace RigidObject
{

	class RigidCircle
	{
	protected:
		double r; // radius
		double density;

		union // position
		{
			struct { double x, y, ang; };
			Point2D centre;
		};

		union // acceleration
		{
			struct { double ax, ay, a_ang; };
			struct { size_t iax, iay, ia_ang; };
			Vector2D acceleration;
		};

		union // velocity
		{
			struct { double vx, vy, v_ang; };
			struct { size_t ivx, ivy, iv_ang; };
			Vector2D velocity;
		};

		union // contact force
		{
			struct { double fx_cont, fy_cont, m_cont; };
			Force2D cont_force;
		};

		// boundary condition
		union // external force
		{
			struct { double fx_ext, fy_ext, m_ext; };
			Force2D ext_force;
		};

		size_t ax_bc_mask, ay_bc_mask, a_ang_bc_mask;
		size_t vx_bc_mask, vy_bc_mask, v_ang_bc_mask;

		union
		{
			struct { double ax_bc, ay_bc, a_ang_bc; };
			struct { size_t iax_bc, iay_bc, ia_ang_bc; };
			Vector2D acceleration_bc;
		};
		union
		{
			struct { double vx_bc, vy_bc, v_ang_bc; };
			struct { size_t ivx_bc, ivy_bc, iv_ang_bc; };
			Vector2D velocity_bc;
		};

		// calculation variables
		Rect lbbox;
		double r2, m, moi; // moment of inertia
		double inv_m, inv_moi;

	public:
		RigidCircle() : r(0.0), density(1.0),
			x(0.0), y(0.0), ang(0.0),
			ax(0.0), ay(0.0), a_ang(0.0),
			vx(0.0), vy(0.0), v_ang(0.0),
			fx_cont(0.0), fy_cont(0.0), m_cont(0.0),
			fx_ext(0.0), fy_ext(0.0), m_ext(0.0),
			ax_bc_mask(0), ay_bc_mask(0), a_ang_bc_mask(0),
			vx_bc_mask(0), vy_bc_mask(0), v_ang_bc_mask(0),
			ax_bc(0.0), ay_bc(0.0), a_ang_bc(0.0),
			vx_bc(0.0), vy_bc(0.0), v_ang_bc(0.0),
			r2(0.0), m(0.0), moi(0.0) {}
		~RigidCircle() {}

		inline double get_radius() const noexcept { return r; }
		inline double get_density() const noexcept { return density; }
		inline double get_mass() const noexcept { return m; }
		inline double get_moi() const noexcept { return moi; }
		inline double get_x() const noexcept { return x; }
		inline double get_y() const noexcept { return y; }
		inline double get_ang() const noexcept { return ang; }
		inline const Point2D& get_centre() const noexcept { return centre; }
		inline double get_ax() const noexcept { return ax; }
		inline double get_ay() const noexcept { return ay; }
		inline double get_a_ang() const noexcept { return a_ang; }
		inline const Vector2D& get_acceleration() const noexcept { return acceleration; }
		inline double get_vx() const noexcept { return vx; }
		inline double get_vy() const noexcept { return vy; }
		inline double get_v_ang() const noexcept { return v_ang; }
		inline const Vector2D& get_velocity() const noexcept { return velocity; }
		inline double get_fx_cont() const noexcept { return fx_cont; }
		inline double get_fy_cont() const noexcept { return fy_cont; }
		inline double get_m_cont() const noexcept { return m_cont; }
		inline const Force2D& get_cont_force() const noexcept { return cont_force; }
		inline double get_fx_ext() const noexcept { return fx_ext; }
		inline double get_fy_ext() const noexcept { return fy_ext; }
		inline double get_m_ext() const noexcept { return m_ext; }
		inline const Force2D& get_ext_force() const noexcept { return ext_force; }

		inline bool has_ax_bc() const noexcept { return ax_bc_mask != 0; }
		inline bool has_ay_bc() const noexcept { return ay_bc_mask != 0; }
		inline bool has_a_ang_bc() const noexcept { return a_ang_bc_mask != 0; }
		inline bool has_vx_bc() const noexcept { return vx_bc_mask != 0; }
		inline bool has_vy_bc() const noexcept { return vy_bc_mask != 0; }
		inline bool has_v_ang_bc() const noexcept { return v_ang_bc_mask != 0; }
		inline double get_ax_bc() const noexcept { return ax_bc; }
		inline double get_ay_bc() const noexcept { return ay_bc; }
		inline double get_a_ang_bc() const noexcept { return a_ang_bc; }
		inline double get_vx_bc() const noexcept { return vx_bc; }
		inline double get_vy_bc() const noexcept { return vy_bc; }
		inline double get_v_ang_bc() const noexcept { return v_ang_bc; }

		void set_ax_bc(double _a) noexcept { ax_bc = _a; }
		void set_ay_bc(double _a) noexcept { ay_bc = _a; }
		void set_a_ang_bc(double _a) noexcept { a_ang_bc = _a; }
		void set_vx_bc(double _v) noexcept { vx_bc = _v; }
		void set_vy_bc(double _v) noexcept { vy_bc = _v; }
		void set_v_ang_bc(double _v) noexcept { v_ang_bc = _v; }
		void set_vbc(double _vx, double _vy, double _vang) noexcept;

		void init(double _x, double _y, double _r, double _density = 1.0) noexcept;
		void set_cont_force(Force2D& cont_force) noexcept;

		inline void get_global_point(const Point2D& lp, Point2D& gp) const noexcept
		{
			gp.x = lp.x + x; gp.y = lp.y + y;
		}

		inline void get_local_point(const Point2D& gp, Point2D& lp) const noexcept
		{
			lp.x = gp.x - x; lp.y = gp.y - y;
		}

		inline void get_global_vector(const Vector2D& lv, Vector2D& gv) const noexcept
		{
			gv.x = lv.x; gv.y = lv.y;
		}

		inline void get_local_vector(const Vector2D& gv, Vector2D& lv) const noexcept
		{
			lv.x = gv.x; lv.y = gv.y;
		}

		inline void get_bbox(
			Rect& bbox,
			double exp_size = 0.0
		) const noexcept
		{
			bbox.xl = x + lbbox.xl - exp_size;
			bbox.xu = x + lbbox.xu + exp_size;
			bbox.yl = y + lbbox.yl - exp_size;
			bbox.yu = y + lbbox.yu + exp_size;
		}

		inline void reset_cont_force() noexcept
		{
			fx_cont = 0.0; fy_cont = 0.0; m_cont = 0.0;
		}

		inline void add_cont_force(
			double _x,
			double _y,
			double fx,
			double fy
		) noexcept
		{
			fx_cont += fx; fy_cont += fy;
			m_cont += (_x - x) * fy - (_y - y) * fx;
		}

		inline bool detect_collision_with_point(
			double p_x,
			double p_y,
			double p_r,
			double& dist,
			Vector2D& lnorm,
			Point2D& lcontpos
		) const noexcept
		{
			double x_diff = p_x - x, y_diff = p_y - y;
			double len = sqrt(x_diff * x_diff + y_diff * y_diff);
			dist = len - p_r;
			if (dist > r) // not overlapping
				return false;
			dist = r - dist;
			lnorm.x = x_diff / len;
			lnorm.y = y_diff / len;
			lcontpos.x = x_diff;
			lcontpos.y = y_diff;
			return true;
		}

		// update rigid circle motion and position
		void update_motion(double dt) noexcept
		{
			// update a
			ax = (fx_cont + fx_ext) * inv_m;
			ay = (fy_cont + fy_ext) * inv_m;
			a_ang = (m_cont + m_ext) * inv_moi;
			iax = (iax & ~ax_bc_mask) | (iax_bc & ax_bc_mask);
			iay = (iay & ~ay_bc_mask) | (iay_bc & ay_bc_mask);
			ia_ang = (ia_ang & ~a_ang_bc_mask) | (ia_ang_bc & a_ang_bc_mask);
			// update velocity
			vx += ax * dt;
			vy += ay * dt;
			v_ang += a_ang * dt;
			ivx = (ivx & ~vx_bc_mask) | (ivx_bc & vx_bc_mask);
			ivy = (ivy & ~vy_bc_mask) | (ivy_bc & vy_bc_mask);
			iv_ang = (iv_ang & ~v_ang_bc_mask) | (iv_ang_bc & v_ang_bc_mask);
			// update position
			x += vx * dt;
			y += vy * dt;
			ang += v_ang * dt;
		}

		void set_init_state(double _r, double _den,
			double _ax, double _ay, double _a_ang,
			double _vx, double _vy, double _v_ang,
			double _x, double _y, double _ang,
			const Force2D& f_cont, const Force2D& f_ext) noexcept
		{
			r = _r;	density = _den;
			ax = _ax; ay = _ay; a_ang = _a_ang;
			vx = _vx; vy = _vy; v_ang = _v_ang;
			x = _x; y = _y; ang = _ang;
			fx_cont = f_cont.fx;
			fy_cont = f_cont.fy;
			m_cont = f_cont.m;
			fx_ext = f_ext.fx;
			fy_ext = f_ext.fy;
			m_ext = f_ext.m;
		}
	};

}

#endif