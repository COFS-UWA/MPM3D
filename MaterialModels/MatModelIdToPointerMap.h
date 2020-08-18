#ifndef __Mat_Model_Id_To_Pointer_Map_h__
#define __Mat_Model_Id_To_Pointer_Map_h__

#include <unordered_map>

#include "MatModelContainer.h"

namespace MatModel
{
	class MatModelIdToPointerMap
	{
	protected:
		typedef std::unordered_map<size_t, MaterialModel*> Id2PtMap;
		Id2PtMap map;

	public:
		MatModelIdToPointerMap();
		~MatModelIdToPointerMap();

		MatModelIdToPointerMap(MatModelContainer& mc);

		void init(MatModelContainer &mc);

		inline MaterialModel* get_mm_by_id(size_t id)
		{
			Id2PtMap::iterator iter = map.find(id);
			if (iter == map.end())
				return nullptr;
			return iter->second;
		}
	};

}

#endif