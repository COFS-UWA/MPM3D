#ifndef __Model_Container_h__
#define __Model_Container_h__

#include "ItemBuffer.hpp"

#include "LinearElasticity.h"
#include "ModifiedCamClay.h"

namespace MatModel
{
	template <typename CModel>
	class __Model_Container__
	{
	protected:
		MemoryUtils::ItemBuffer<CModel> buffer;
		ConstitutiveModelList list;

	public:
		CModel *add(size_t num)
		{
			CModel *res = buffer.alloc(num);
			for (size_t m_id = 0; m_id < num; ++m_id)
			{
				new (res + m_id) CModel;
				list.add_cm(res[m_id]);
			}
			return res;
		}
		inline size_t get_num(void) { return list.get_num(); }
		inline CModel *first(void) { return static_cast<CModel *>(list.first()); }
		inline CModel *next(CModel *cm) { return static_cast<CModel *>(list.next(cm)); }
	};

#define __Model_Container_Add_Model__(ModelName) \
protected:                                       \
	__Model_Container__<ModelName> ModelName##Container; \
public:                                          \
	ModelName *add_##ModelName##(size_t num)     \
	{                                            \
		return ModelName##Container.add(num);    \
	}                                            \
	inline size_t get_num_##ModelName##(void)    \
	{                                            \
		return ModelName##Container.get_num();   \
	}                                            \
	inline ModelName *first_##ModelName##(void)  \
	{                                            \
		return ModelName##Container.first();     \
	}                                            \
	inline ModelName *next_##ModelName##(ModelName *cm) \
	{                                            \
		return ModelName##Container.next(cm);    \
	}

	class ModelContainer
	{
		__Model_Container_Add_Model__(LinearElasticity);
		__Model_Container_Add_Model__(ModifiedCamClay);
	};

#undef __Model_Container_Add_Model__
};

#endif