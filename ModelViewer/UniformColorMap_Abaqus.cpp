#include "ModelViewer_pcp.h"

#include "UniformColorMap_Abaqus.h"

UniformColourMap_Abaqus::UniformColourMap_Abaqus()
{
	unsigned char abaqus_color_map[][3] = {
		{ 0,   0,   255 },
		{ 0,   93,  255 },
		{ 0,   185, 255 },
		{ 0,   255, 232 },
		{ 0,   255, 139 },
		{ 0,   255, 46  },
		{ 46,  255, 0   },
		{ 139, 255, 0   },
		{ 232, 255, 0   },
		{ 255, 185, 0   },
		{ 255, 93,  0   },
		{ 255, 0,   0   }
	};
	size_t abaqus_color_num = sizeof(abaqus_color_map) / sizeof(abaqus_color_map[0]);
	set_colour(abaqus_color_map, abaqus_color_num);
}

UniformColourMap_Abaqus::~UniformColourMap_Abaqus()
{

}

