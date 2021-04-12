#ifndef __Mat_Model_Container_h__
#define __Mat_Model_Container_h__

#include "ItemBuffer.hpp"
#include "LinkList.hpp"

#include "LinearElasticity.h"
#include "ModifiedCamClay.h"
#include "UndrainedModifiedCamClay.h"
#include "VonMises.h"
#include "Tresca.h"
#include "SandHypoplasticityByUmat.h"
#include "SandHypoplasticityWrapper.h"

namespace MatModel
{
	class MatModelContainer;

	template <typename MModel>
	class __Mat_Model_Container__
	{
	protected:
		MemoryUtils::ItemBuffer<MModel> buffer;
		LinkList<MModel, offsetof(MModel, pointer_by_container)> list;

	public:
		__Mat_Model_Container__() {}
		~__Mat_Model_Container__()
		{
			list.reset();
			buffer.clear();
		}

		MModel *add(size_t num, MatModelContainer *container)
		{
			MModel *res = buffer.alloc(num);
			MModel *iter = res;
			for (size_t m_id = 0; m_id < num; ++m_id)
			{
				new (iter) MModel;
				iter->id = container->gen_uid();
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
	__Mat_Model_Container__<ModelName> ModelName##Container; \
public:                                          \
	ModelName *add_##ModelName##(size_t num)     \
	{                                            \
		return ModelName##Container.add(num, this);    \
	}                                            \
	inline size_t get_num_##ModelName##()        \
	{                                            \
		return ModelName##Container.get_num();   \
	}                                            \
	inline ModelName *first_##ModelName##()      \
	{                                            \
		return ModelName##Container.first();     \
	}                                            \
	inline ModelName *next_##ModelName##(ModelName *mm) \
	{                                            \
		return ModelName##Container.next(mm);    \
	}                                            \
	inline bool is_not_end_##ModelName##(ModelName *mm) \
	{                                            \
		return ModelName##Container.is_not_end(mm);  \
	}

	class MatModelContainer
	{
	protected:
		template <typename MModel> friend class __Mat_Model_Container__;
		static size_t cur_uid;
		inline size_t gen_uid() noexcept { return cur_uid++; }

		__Add_Mat_Model_to_Model_Container__(LinearElasticity);
		__Add_Mat_Model_to_Model_Container__(ModifiedCamClay);
		__Add_Mat_Model_to_Model_Container__(UndrainedModifiedCamClay);
		__Add_Mat_Model_to_Model_Container__(VonMises);
		__Add_Mat_Model_to_Model_Container__(Tresca);
		__Add_Mat_Model_to_Model_Container__(SandHypoplasticityByUmat);
		__Add_Mat_Model_to_Model_Container__(SandHypoplasticityWrapper);
	};

#undef __Add_Mat_Model_to_Model_Container__
};

#endif