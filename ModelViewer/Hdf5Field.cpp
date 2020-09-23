#include "ModelViewer_pcp.h"

#include "Hdf5Field.h"

namespace Hdf5Field
{
	namespace Hdf5Field_internal
	{
		const char * field_name[] = {
			"x", // 0
			"y", // 1
			"z", // 2
			"vol", // 3
			"m", // 4
			"density", //5
			"vx", // 6
			"vy", // 7
			"vz", // 8
			"s11", // 9
			"s22", // 10
			"s33", // 11
			"s12", // 12
			"s23", // 13
			"s31", // 14
			"e11", // 15
			"e22", // 16
			"e33", // 17
			"e12", // 18
			"e23", // 19
			"e31", // 20
			"m_s", // 21
			"m_f", // 22
			"density_s", // 23
			"density_f", // 24
			"vol_m", // 25
			"vx_s", // 26
			"vy_s", // 27
			"vz_s", // 28
			"vx_f", // 29
			"vy_f", // 30
			"vz_f", // 31
			"p", // 32
			"Mises Strain" // 33
		};
	}
}