#ifndef __Mat_Model_Container_h__
#define __Mat_Model_Container_h__

#include "CLAItemBuffer.hpp"
#include "LinkList.hpp"

#include "LinearElasticity.h"
#include "ModifiedCamClay.h"
#include "UndrainedModifiedCamClay.h"
#include "VonMises.h"
#include "Tresca.h"
#include "MohrCoulombWrapper.h"
#include "SandHypoplasticityByUmat.h"
#include "SandHypoplasticityWrapper.h"
#include "SandHypoplasticityStbWrapper.h"
#include "NorsandWrapper.h"

namespace MatModel
{
	class MatModelContainer;

	template <typename MModel>
	class __Mat_Model_Container__
	{
	protected:
		MemoryUtils::CLAItemBuffer<MModel> buffer;
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
				iter = buffer.following_item(iter);
			}
			return res;
		}

		inline size_t get_num() { return list.get_num(); }

		inline MModel *first() { return list.first(); }
		inline MModel *next(MModel *mm) { return list.next(mm); }
		inline bool is_not_end(MModel *mm) { return list.is_not_end(mm); }
		inline MModel* following_model(MModel* mm) { return buffer.following_item(mm); }
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
	}                                             \
	inline ModelName* following_##ModelName##(ModelName* mm) \
	{                                            \
		return ModelName##Container.following_model(mm);  \
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
		__Add_Mat_Model_to_Model_Container__(MohrCoulombWrapper);
		__Add_Mat_Model_to_Model_Container__(SandHypoplasticityByUmat);
		__Add_Mat_Model_to_Model_Container__(SandHypoplasticityWrapper);
		__Add_Mat_Model_to_Model_Container__(SandHypoplasticityStbWrapper);
		__Add_Mat_Model_to_Model_Container__(NorsandWrapper);
	};

#undef __Add_Mat_Model_to_Model_Container__
};

#endif