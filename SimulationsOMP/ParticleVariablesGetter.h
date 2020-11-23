#ifndef __Particle_Variables_Getter_h__
#define __Particle_Variables_Getter_h__

namespace MatModel { class MaterialModel; }

class ParticleVariablesGetter
{
public:
	virtual size_t get_index() const noexcept { return 0; }
	virtual double get_m() const noexcept { return 0.0; }
	virtual double get_bfx() const noexcept { return 0.0; }
	virtual double get_bfy() const noexcept { return 0.0; }
	virtual double get_bfz() const noexcept { return 0.0; }
	virtual double get_tx() const noexcept { return 0.0; }
	virtual double get_ty() const noexcept { return 0.0; }
	virtual double get_tz() const noexcept { return 0.0; }
	virtual double get_x() const noexcept { return 0.0; }
	virtual double get_y() const noexcept { return 0.0; }
	virtual double get_z() const noexcept { return 0.0; }
	virtual double get_vol() const noexcept { return 0.0; }
	virtual double get_square_r() const noexcept { return 0.0; }
	virtual double get_circle_r() const noexcept { return 0.0; }
	virtual MatModel::MaterialModel* get_mat_model() const noexcept
	{ return nullptr; }
	virtual double get_density() const noexcept { return 0.0; }
	virtual double get_vx() const noexcept { return 0.0; }
	virtual double get_vy() const noexcept { return 0.0; }
	virtual double get_vz() const noexcept { return 0.0; }
	virtual double get_ux() const noexcept { return 0.0; }
	virtual double get_uy() const noexcept { return 0.0; }
	virtual double get_uz() const noexcept { return 0.0; }
	virtual double get_s11() const noexcept { return 0.0; }
	virtual double get_s22() const noexcept { return 0.0; }
	virtual double get_s33() const noexcept { return 0.0; }
	virtual double get_s12() const noexcept { return 0.0; }
	virtual double get_s23() const noexcept { return 0.0; }
	virtual double get_s31() const noexcept { return 0.0; }
	virtual double get_e11() const noexcept { return 0.0; }
	virtual double get_e22() const noexcept { return 0.0; }
	virtual double get_e33() const noexcept { return 0.0; }
	virtual double get_e12() const noexcept { return 0.0; }
	virtual double get_e23() const noexcept { return 0.0; }
	virtual double get_e31() const noexcept { return 0.0; }
	virtual double get_ee11() const noexcept { return 0.0; }
	virtual double get_ee22() const noexcept { return 0.0; }
	virtual double get_ee33() const noexcept { return 0.0; }
	virtual double get_ee12() const noexcept { return 0.0; }
	virtual double get_ee23() const noexcept { return 0.0; }
	virtual double get_ee31() const noexcept { return 0.0; }
	virtual double get_pe11() const noexcept { return 0.0; }
	virtual double get_pe22() const noexcept { return 0.0; }
	virtual double get_pe33() const noexcept { return 0.0; }
	virtual double get_pe12() const noexcept { return 0.0; }
	virtual double get_pe23() const noexcept { return 0.0; }
	virtual double get_pe31() const noexcept { return 0.0; }
	virtual double get_N1() const noexcept { return 0.0; }
	virtual double get_N2() const noexcept { return 0.0; }
	virtual double get_N3() const noexcept { return 0.0; }
	virtual double get_N4() const noexcept { return 0.0; }
};

#endif