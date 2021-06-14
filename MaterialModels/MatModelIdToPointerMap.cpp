#include "MaterialModels_pcp.h"

#include <exception>

#include "MatModelIdToPointerMap.h"

namespace MatModel
{
	MatModelIdToPointerMap::MatModelIdToPointerMap() {}
	
	MatModelIdToPointerMap::~MatModelIdToPointerMap() {}

	MatModelIdToPointerMap::MatModelIdToPointerMap(MatModelContainer& mc) { init(mc); }

	void MatModelIdToPointerMap::init(MatModelContainer& mc)
	{
#define EXCEPTION_MSG_LEN 200
		char exception_msg[EXCEPTION_MSG_LEN];
		std::pair<Id2PtMap::iterator, bool> res;

		// linear elasticity
		for (MatModel::LinearElasticity* iter = mc.first_LinearElasticity();
			 mc.is_not_end_LinearElasticity(iter);
			 iter = mc.next_LinearElasticity(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(
					exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: Linear elasticity model %zu already exists.",
					iter->get_id()
					);
				throw std::exception(exception_msg);
			}
		}

		// modified cam clay
		for (MatModel::ModifiedCamClay* iter = mc.first_ModifiedCamClay();
			mc.is_not_end_ModifiedCamClay(iter);
			iter = mc.next_ModifiedCamClay(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(
					exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: Modified Cam-clay model %zu already exists.",
					iter->get_id()
					);
				throw std::exception(exception_msg);
			}
		}

		// undrained modified cam clay
		for (MatModel::UndrainedModifiedCamClay* iter = mc.first_UndrainedModifiedCamClay();
			mc.is_not_end_UndrainedModifiedCamClay(iter);
			iter = mc.next_UndrainedModifiedCamClay(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(
					exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: Undrained Modified Cam-clay model %zu already exists.",
					iter->get_id()
					);
				throw std::exception(exception_msg);
			}
		}

		// VonMises
		for (MatModel::VonMises *iter = mc.first_VonMises();
			mc.is_not_end_VonMises(iter);
			iter = mc.next_VonMises(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(
					exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: Von Mises model %zu already exists.",
					iter->get_id()
				);
				throw std::exception(exception_msg);
			}
		}

		// Tresca
		for (MatModel::Tresca *iter = mc.first_Tresca();
			mc.is_not_end_Tresca(iter);
			iter = mc.next_Tresca(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(
					exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: Tresca model %zu already exists.",
					iter->get_id()
					);
				throw std::exception(exception_msg);
			}
		}

		// Mohr Coulomb
		for (MatModel::MohrCoulombWrapper* iter = mc.first_MohrCoulombWrapper();
			mc.is_not_end_MohrCoulombWrapper(iter);
			iter = mc.next_MohrCoulombWrapper(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: MohrCoulomb model %zu already exists.",
					iter->get_id());
				throw std::exception(exception_msg);
			}
		}

		// Sandhypoplasticity
		for (MatModel::SandHypoplasticityWrapper* iter = mc.first_SandHypoplasticityWrapper();
			mc.is_not_end_SandHypoplasticityWrapper(iter);
			iter = mc.next_SandHypoplasticityWrapper(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: SandHypoplasticity model %zu already exists.",
					iter->get_id());
				throw std::exception(exception_msg);
			}
		}
		
		// SandhypoplasticityStb
		for (MatModel::SandHypoplasticityStbWrapper* iter = mc.first_SandHypoplasticityStbWrapper();
			mc.is_not_end_SandHypoplasticityStbWrapper(iter);
			iter = mc.next_SandHypoplasticityStbWrapper(iter))
		{
			res = map.emplace(iter->get_id(), iter);
			if (!res.second)
			{
				snprintf(exception_msg,
					EXCEPTION_MSG_LEN,
					"class MatModelIdToPointerMap error: SandHypoplasticityStb model %zu already exists.",
					iter->get_id());
				throw std::exception(exception_msg);
			}
		}
	}
}