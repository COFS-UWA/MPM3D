#ifndef __Hdf5_Field_h__
#define __Hdf5_Field_h__

#include "Hdf5FieldExtraction.h"
#include "Hdf5FieldExtraction_x.h"
#include "Hdf5FieldExtraction_y.h"
#include "Hdf5FieldExtraction_z.h"
#include "Hdf5FieldExtraction_vol.h"
#include "Hdf5FieldExtraction_m.h"
#include "Hdf5FieldExtraction_density.h"
#include "Hdf5FieldExtraction_vx.h"
#include "Hdf5FieldExtraction_vy.h"
#include "Hdf5FieldExtraction_vz.h"
#include "Hdf5FieldExtraction_s11.h"
#include "Hdf5FieldExtraction_s22.h"
#include "Hdf5FieldExtraction_s33.h"
#include "Hdf5FieldExtraction_s12.h"
#include "Hdf5FieldExtraction_s23.h"
#include "Hdf5FieldExtraction_s31.h"
#include "Hdf5FieldExtraction_e11.h"
#include "Hdf5FieldExtraction_e22.h"
#include "Hdf5FieldExtraction_e33.h"
#include "Hdf5FieldExtraction_e12.h"
#include "Hdf5FieldExtraction_e23.h"
#include "Hdf5FieldExtraction_e31.h"
#include "Hdf5FieldExtraction_mises_strain_2d.h"
#include "Hdf5FieldExtraction_ee11.h"
#include "Hdf5FieldExtraction_ee22.h"
#include "Hdf5FieldExtraction_ee33.h"
#include "Hdf5FieldExtraction_ee12.h"
#include "Hdf5FieldExtraction_ee23.h"
#include "Hdf5FieldExtraction_ee31.h"
#include "Hdf5FieldExtraction_pe11.h"
#include "Hdf5FieldExtraction_pe22.h"
#include "Hdf5FieldExtraction_pe33.h"
#include "Hdf5FieldExtraction_pe12.h"
#include "Hdf5FieldExtraction_pe23.h"
#include "Hdf5FieldExtraction_pe31.h"
#include "Hdf5FieldExtraction_plastic_mises_strain_2d.h"
// chm
#include "Hdf5FieldExtraction_m_s.h"
#include "Hdf5FieldExtraction_m_f.h"
#include "Hdf5FieldExtraction_density_s.h"
#include "Hdf5FieldExtraction_density_f.h"
#include "Hdf5FieldExtraction_vol_m.h"
#include "Hdf5FieldExtraction_vx_s.h"
#include "Hdf5FieldExtraction_vy_s.h"
#include "Hdf5FieldExtraction_vz_s.h"
#include "Hdf5FieldExtraction_vx_f.h"
#include "Hdf5FieldExtraction_vy_f.h"
#include "Hdf5FieldExtraction_vz_f.h"
#include "Hdf5FieldExtraction_p.h"

namespace Hdf5Field
{
	enum FieldType : unsigned int
	{
		x = 0,
		y = 1,
		z = 2,
		vol = 3,
		m = 4,
		density = 5,
		vx = 6,
		vy = 7,
		vz = 8,
		s11 = 9,
		s22 = 10,
		s33 = 11,
		s12 = 12,
		s23 = 13,
		s31 = 14,
		e11 = 15,
		e22 = 16,
		e33 = 17,
		e12 = 18,
		e23 = 19,
		e31 = 20,
		m_s = 21,
		m_f = 22,
		density_s = 23,
		density_f = 24,
		vol_m = 25,
		vx_s = 26,
		vy_s = 27,
		vz_s = 28,
		vx_f = 29,
		vy_f = 30,
		vz_f = 31,
		p = 32,
		mises_strain_2d = 33,
		ee11 = 34,
		ee22 = 35,
		ee33 = 36,
		ee12 = 37,
		ee23 = 38,
		ee31 = 39,
		pe11 = 40,
		pe22 = 41,
		pe33 = 42,
		pe12 = 43,
		pe23 = 44,
		pe31 = 45,
		plastic_mises_strain_2d = 46
	};

	namespace Hdf5Field_internal
	{
		const char * field_name[];

		typedef Hdf5FieldExtraction* (*MakeFunc)();

		template <class Field>
		Hdf5FieldExtraction* make_Hdf5FieldExtraction_template() { return new Field; }

		const MakeFunc make_funcs[] = {
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_x>, // 0
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_y>, // 1
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_z>, // 2
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vol>, // 3
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_m>, // 4
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_density>, // 5
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vx>, // 6
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vy>, // 7
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vz>, // 8
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_s11>, // 9
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_s22>, // 10
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_s33>, // 11
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_s12>, // 12
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_s23>, // 13
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_s31>, // 14
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_e11>, // 15
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_e22>, // 16
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_e33>, // 17
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_e12>, // 18
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_e23>, // 19
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_e31>, // 20
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_m_s>, // 21
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_m_f>, // 22
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_density_s>, // 23
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_density_f>, // 24
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vol_m>, // 25
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vx_s>, // 26
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vy_s>, // 27
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vz_s>, // 28
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vx_f>, // 29
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vy_f>, // 30
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vz_f>, // 31
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_p>, // 32
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_mises_strain_2d>, // 33
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_ee11>, // 34
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_ee22>, // 35
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_ee33>, // 36
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_ee12>, // 37
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_ee23>, // 38
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_ee31>, // 39
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_pe11>, // 40
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_pe22>, // 41
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_pe33>, // 42
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_pe12>, // 43
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_pe23>, // 44
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_pe31>, // 45
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_plastic_mises_strain_2d> // 46
		};
		const size_t make_func_num = sizeof(make_funcs) / sizeof(make_funcs[0]);
	}

	inline Hdf5FieldExtraction* make(FieldType type)
	{
		if (type < Hdf5Field_internal::make_func_num)
			return (*Hdf5Field_internal::make_funcs[unsigned int(type)])();
		return nullptr;
	}

	inline const char* get_name(FieldType type)
	{
		if (type < Hdf5Field_internal::make_func_num)
			return Hdf5Field_internal::field_name[unsigned int(type)];
		return nullptr;
	}
};

#endif