#ifndef __Uniform_Color_Map_Mono_h__
#define __Uniform_Color_Map_Mono_h__

#include "UniformColorMap.h"

class UniformColorMap_Mono : public UniformColorMap
{
public:
	UniformColorMap_Mono();
	~UniformColorMap_Mono();	
	void set_color(float r, float g, float b);
};

#endif