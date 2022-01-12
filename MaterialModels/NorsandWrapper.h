#ifndef __Norsand_Wrapper_h__
#define __Norsand_Wrapper_h__

#include "MaterialModel.h"
#include "norsand.h"

namespace MatModel
{
	int nor_sand_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);
	void nor_sand_wrapper_store_to_function(const MaterialModel* _self, char* mat_model_mem);
	void nor_sand_wrapper_retrieve_from_function(MaterialModel* _self, const char* mat_model_mem);

	class NorsandWrapper : public MaterialModel
	{
		friend int nor_sand_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);
		friend void nor_sand_wrapper_store_to_function(const MaterialModel* _self, char* mat_model_mem);
		friend void nor_sand_wrapper_retrieve_from_function(MaterialModel* _self, const char* mat_model_mem);

	protected:
		NorsandGlobal glb;
		Norsand mat;

	public:
		NorsandWrapper();
		NorsandWrapper(const NorsandWrapper& other);
		~NorsandWrapper() {}
		
		void set_param(const double* ini_stress,
			double e0,
			double phi,
			double gamma, double lambda,
			double N, double chi, double H,
			double Ig, double niu);
		void set_pi(double pi) noexcept { mat.pi = pi; }
		
		double get_phi() const noexcept { return glb.phi; }
		double get_gamma() const noexcept { return glb.gamma; }
		double get_lambda() const noexcept { return glb.lambda; }
		double get_N() const noexcept { return glb.N; }
		double get_chi() const noexcept { return glb.chi; }
		double get_H() const noexcept { return glb.H; }
		double get_Ig() const noexcept { return glb.Ig; }
		double get_niu() const noexcept { return glb.niu; }

		const double* get_stress() const noexcept { return mat.stress; }
		double get_e() const noexcept { return mat.e; }
		double get_pi() const noexcept { return mat.pi; }
	};
}

#endif