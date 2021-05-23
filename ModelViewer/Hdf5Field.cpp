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
			"volume", // 3
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
			"solid mass", // 21
			"fluid mass", // 22
			"solid density", // 23
			"fluid density", // 24
			"mixture volume", // 25
			"solid vx", // 26
			"solid vy", // 27
			"solid vz", // 28
			"fluid vx", // 29
			"fluid vy", // 30
			"fluid vz", // 31
			"pore pressure", // 32
			"mises strain", // 33
			"ee11", // 34
			"ee22", // 35
			"ee33", // 36
			"ee12", // 37
			"ee23", // 38
			"ee31", // 39
			"pe11",  // 40
			"pe22",  // 41
			"pe33",  // 42
			"pe12",  // 43
			"pe23",  // 44
			"pe31",  // 45
			"plastic mises strain 2d", // 46
			"maximum shear stress", // 47
			"mat_e" // 48
		};
	}
}
