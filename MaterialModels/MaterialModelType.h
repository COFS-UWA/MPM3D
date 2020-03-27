#ifndef __Material_Model_Type_H__
#define __Material_Model_Type_H__

namespace MatModel
{
	enum class Type : unsigned short
	{
		InvalidType = 0,
		LinearElasticity = 1,
		ModifiedCamClay = 2
	};
}

#endif