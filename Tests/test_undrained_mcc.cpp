#include "Tests_pcp.h"

#include <fstream>

#include "UndrainedModifiedCamClay.h"
#include "test_material_models.h"

void test_undrained_mcc()
{
	enum class AnalysisType : unsigned char
	{
		TriaxialUndrained = 0,
		IsotropicConsolidation = 1
	};

	AnalysisType tp = AnalysisType::TriaxialUndrained;
	//AnalysisType tp = AnalysisType::IsotropicConsolidation;
	double de = -0.5;
	//double de = -0.01;
	double ini_stress[6] = { -20000.0, -12025.0, -12025.0, 0.0, 0.0, 0.0 };
	//double ini_stress[6] = { -20000.0, -20000.0, -20000.0, 0.0, 0.0, 0.0 };
	size_t inc_num = 10000;
	size_t out_num = 100;
	MatModel::UndrainedModifiedCamClay mcc;
	mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.667, ini_stress, 20000.0, 1.0e8);
	std::string res_csv_name = "undrained_mcc_res.csv";
	
	const double (*Dep_mat)[6];
	double dstrain[6];
	const double *stress, *eff_stress;
	double e11, e22, e33;

	de /= double(inc_num);
	e11 = 0.0;
	e22 = 0.0;
	e33 = 0.0;
	size_t out_inv = inc_num / out_num;
	if (out_inv < 1)
		out_inv = 1;
	int res;

	std::fstream res_file;
	res_file.open(res_csv_name.c_str(), std::ios::out | std::ios::binary);
	stress = mcc.get_stress();
	eff_stress = mcc.get_effective_stress();
	res_file << "e11, e22, e33, s11, s22, s33, p, s11', s22', s33', p', q, pore, pc, e, e_NCRC, f, res\n"
			 << e11 << ", " << e22 << ", " << e33 << ", "
			 << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
			 << mcc.get_p() << ", "
			 << eff_stress[0] << ", " << eff_stress[1] << ", " << eff_stress[2] << ", "
			 << mcc.get_effective_p() << ", "
			 << mcc.get_q() << ", " << mcc.get_pore_pressure() << ", "
			 << mcc.get_pc() << ", " << mcc.get_e_by_strain() << ", " << mcc.get_e_by_model() << ", "
			 << mcc.get_f() << ", 0,\n";

	for (size_t inc_id = 0; inc_id < inc_num; ++inc_id)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialUndrained:
			// to be improved
			Dep_mat = reinterpret_cast<const double(*)[6]>(mcc.get_Dep_mat());
			dstrain[0] = de;
			dstrain[1] = -Dep_mat[1][0] / (Dep_mat[1][1] + Dep_mat[1][2]) * de;
			dstrain[2] = -Dep_mat[2][0] / (Dep_mat[2][1] + Dep_mat[2][2]) * de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::IsotropicConsolidation:
			dstrain[0] = de;
			dstrain[1] = de;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		default:
			res_file.close();
			return;
		}

		res = mcc.integrate(dstrain);

		// output result
		e11 += dstrain[0];
		e22 += dstrain[1];
		e33 += dstrain[2];
		stress = mcc.get_stress();
		//if (inc_id % out_inv == 0)
		//{
			res_file << e11 << ", " << e22 << ", " << e33 << ", "
				<< stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
				<< mcc.get_p() << ", "
				<< eff_stress[0] << ", " << eff_stress[1] << ", " << eff_stress[2] << ", "
				<< mcc.get_effective_p() << ", "
				<< mcc.get_q() << ", " << mcc.get_pore_pressure() << ", "
				<< mcc.get_pc() << ", " << mcc.get_e_by_strain() << ", " << mcc.get_e_by_model() << ", "
				<< mcc.get_f() << ", 0,\n";
		//}
	}

	res_file.close();
}
