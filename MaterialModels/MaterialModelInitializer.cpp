#include "MaterialModels_pcp.h"

#include "SandHypoplasticityByUmat.h"
#include "MaterialModelInitializer.h"

namespace MatModel
{
	bool MaterialModelInitializer::is_init = false;

	MaterialModelInitializer MaterialModelInitializer::instance;

	int MaterialModelInitializer::init()
	{
		if (is_init) return 0;

		is_init = true;
		MatModel::SandHypoplasticityByUmat::init(
			"../../Vendors/ExtMatModels/SandHypoplasticity.dll");
		return 0;
	}

	MaterialModelInitializer::~MaterialModelInitializer()
	{
		if (is_init)
		{
			SandHypoplasticityByUmat::free();
			is_init = false;
		}
	}

}