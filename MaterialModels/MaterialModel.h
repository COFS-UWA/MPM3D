#ifndef __Material_Model_h__
#define __Material_Model_h__

#include "LinkList.hpp"

namespace MatModel
{
	class MaterialModel;
	template <typename MModel> class __Mat_Model_Container__;

	typedef int(*CMIntFunc)(MaterialModel *_self, double dstrain[6]);
	// store model to mat_model_mem
	typedef void (*StoreToFunc)(const MaterialModel* _self, char *mat_model_mem);
	// retrieve model from mat_model_mem
	typedef void (*RetrieveFromFunc)(MaterialModel* _self, const char* mat_model_mem);

	class MaterialModel
	{
	protected:
		const char *type;
		size_t id;

		union // stress
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		union // stress increment
		{
			double dstress[6];
			struct { double ds11, ds22, ds33, ds12, ds23, ds31; };
		};
		union // elastic strain increment
		{
			double dstrain_e[6];
			struct { double dee11, dee22, dee33, dee12, dee23, dee31; };
		};
		union // plastic strain increment
		{
			double dstrain_p[6];
			struct { double dep11, dep22, dep33, dep12, dep23, dep31; };
		};
		union // Elastic stiffness
		{
			double De_mat[6][6];
			double De_mat_array[36];
		};
		union // Elasto-plastic stiffness
		{
			double Dep_mat[6][6];
			double Dep_mat_array[36];
		};

		CMIntFunc integration_func;
		StoreToFunc store_to_func;
		RetrieveFromFunc retrieve_from_func;

	public:
		MaterialModel(
			CMIntFunc _inte_func = nullptr,
			const char *_type = "MaterialModel",
			StoreToFunc _sto_func = nullptr,
			RetrieveFromFunc _ret_func = nullptr) :
			integration_func(_inte_func), type(_type),
			store_to_func(_sto_func), retrieve_from_func(_ret_func){}
		~MaterialModel() {}

		inline int integrate(double dstrain[6]) { return (*integration_func)(this, dstrain); }
		inline void store_to(char* mat_model_mem) const noexcept { (*store_to_func)(this, mat_model_mem); }
		inline void retrieve_from(const char* mat_model_mem) noexcept { (*retrieve_from_func)(this, mat_model_mem); }

		inline size_t get_id() const noexcept { return id; }
		inline void set_id(size_t _id) noexcept { id = _id; }

		inline const char *get_type() noexcept { return type; }
		inline const double *get_stress()    noexcept { return stress; }
		inline const double *get_dstress()   noexcept { return dstress; }
		inline const double *get_dstrain_e() noexcept { return dstrain_e; }
		inline const double *get_dstrain_p() noexcept { return dstrain_p; }
		inline const double *get_De_mat()  noexcept { return De_mat_array; }
		inline const double *get_Dep_mat() noexcept { return Dep_mat_array; }

		union
		{
			void *ext_data_pt;
			unsigned long long ext_data_ull;
			long long ext_data_ll;
			double ext_data_d;
		};

	protected:
		template <typename MModel> friend class __Mat_Model_Container__;
		LinkListPointer<MaterialModel> pointer_by_container;
	};
};

#endif