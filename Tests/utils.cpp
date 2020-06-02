#include "Tests_pcp.h"

#include "utils.h"

ValueToColor::Colori ColorScaleExamples::abaqus_color_scale[] = {
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

size_t ColorScaleExamples::abaqus_color_scale_num 
	= sizeof(ColorScaleExamples::abaqus_color_scale)
	/ sizeof(ColorScaleExamples::abaqus_color_scale[0]);