#ifndef __Hdf5_Field_h__
#define __Hdf5_Field_h__

#include "Hdf5FieldExtraction.h"
#include "Hdf5FieldExtraction_x.h"
#include "Hdf5FieldExtraction_y.h"
#include "Hdf5FieldExtraction_z.h"
#include "Hdf5FieldExtraction_vol.h"

namespace Hdf5Field
{
	enum FieldType : unsigned int
	{
		x = 0,
		y = 1,
		z = 2,
		vol = 3,
		m = 4,
		density = 5
	};

	namespace Hdf5Field_internal
	{
		const char * field_name[];

		typedef Hdf5FieldExtraction* (*MakeFunc)();

		template <class Field>
		Hdf5FieldExtraction* make_Hdf5FieldExtraction_template() { return new Field; }

		const MakeFunc make_funcs[] = {
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_x>,
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_y>,
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_z>,
			&make_Hdf5FieldExtraction_template<Hdf5FieldExtraction_vol>
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