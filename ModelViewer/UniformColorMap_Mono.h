#ifndef __Uniform_Color_Map_Mono_h__
#define __Uniform_Color_Map_Mono_h__

#include "UniformColorMap.h"

// The style of this color map is obtained from abaqus 
class UniformColorMap_Mono : public UniformColorMap
{
public:
	UniformColorMap_Mono();
	~UniformColorMap_Mono();

protected:
	void set_color(float r, float g, float b);
};

#endif