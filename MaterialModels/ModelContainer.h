#ifndef __Model_Container_h__
#define __Model_Container_h__

#include "ItemBuffer.hpp"
#include "LinkList.hpp"

#include "LinearElasticity.h"
#include "ModifiedCamClay.h"

namespace MatModel
{
	template <typename MModel>
	class __Model_Container__
	{
	protected:
		MemoryUtils::ItemBuffer<MModel> buffer;
		LinkList<MModel, offsetof(MModel, pointer_by_container)> list;

	public:
		__Model_Container__() {}
		~__Model_Container__()
		{
			list.reset();
			buffer.clear();
		}

		MModel *add(size_t num)
		{
			MModel *res = buffer.alloc(num);
			MModel *iter = res;
			for (size_t m_id = 0; m_id < num; ++m_id)
			{
				new (iter) MModel;
				list.append(iter);
				++iter;
			}
			return res;
		}

		inline size_t get_num() { return list.get_num(); }

		inline MModel *first() { return list.first(); }
		inline MModel *next(MModel *mm) { return list.next(mm); }
		inline bool is_not_end(MModel *mm) { return list.is_not_end(mm); }
	};

#define __Add_Mat_Model_to_Model_Container__(ModelName) \
protected:                                       \
	__Model_Container__<ModelName> ModelName##Container; \
public:                                          \
	ModelName *add_##ModelName##(size_t num)     \
	{                                            \
		return ModelName##Container.add(num);    \
	}                                            \
	inline size_t get_num_##ModelName##()    \
	{                                            \
		return ModelName##Container.get_num();   \
	}                                            \
	inline ModelName *first_##ModelName##()  \
	{                                            \
		return ModelName##Container.first();     \
	}                                            \
	inline ModelName *next_##ModelName##(ModelName *mm) \
	{                                            \
		return ModelName##Container.next(mm);    \
	}                                            \
	inline bool is_not_end(ModelName *mm)        \
	{                                            \
		return ModelName##Container.is_not_end(mm);  \
	}

	class ModelContainer
	{
		__Add_Mat_Model_to_Model_Container__(LinearElasticity);
		__Add_Mat_Model_to_Model_Container__(ModifiedCamClay);
	};

#undef __Add_Mat_Model_to_Model_Container__
};

#endif