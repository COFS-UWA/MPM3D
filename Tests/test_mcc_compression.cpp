#include "Tests_pcp.h"

#include <fstream>

#include "ModifiedCamClay.h"
#include "test_material_models.h"

void test_mcc_compression()
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
		OneDConsolidation = 2,
		IsotropicConsolidation = 3
	};

	AnalysisType tp = AnalysisType::IsotropicConsolidation;
	//AnalysisType tp = AnalysisType::TriaxialDrained;
	double de = -0.2;
	size_t inc_num = 10000;
	size_t out_num = 100;
	double ini_stress[6] = { -10000.0, -10000.0, -10000.0, 0.0, 0.0, 0.0 };
	MatModel::ModifiedCamClay mcc;
	mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.667, ini_stress, 20000.0);
	//mcc.set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
	std::string res_csv_name = "mcc_compression_res.csv";
	
	const double(*Dep_mat)[6];
	double dstrain[6];
	const double *stress;
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
	res_file << "e11, e22, e33, s11, s22, s33, p, q, pc, e, e_NCRC, f, res\n"
			 << e11 << ", " << e22 << ", " << e33 << ", "
			 << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
			 << mcc.get_p() << ", " << mcc.get_q() << ", " << mcc.get_pc() << ", "
			 << mcc.get_e_by_strain() << ", " << mcc.get_e_by_model() << ", "
			 << mcc.get_f() << ", 0,\n";

	for (size_t inc_id = 0; inc_id < inc_num; ++inc_id)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			// to be improved
			Dep_mat = reinterpret_cast<const double(*)[6]>(mcc.get_Dep_mat());
			dstrain[0] = de;
			dstrain[1] = -Dep_mat[1][0] / (Dep_mat[1][1] + Dep_mat[1][2]) * de;
			dstrain[2] = -Dep_mat[2][0] / (Dep_mat[2][1] + Dep_mat[2][2]) * de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::TriaxialUndrained:
			dstrain[0] = de;
			dstrain[1] = -0.5 * de;
			dstrain[2] = -0.5 * de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::OneDConsolidation:
			dstrain[0] = de;
			dstrain[1] = 0.0;
			dstrain[2] = 0.0;
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
				<< mcc.get_p() << ", " << mcc.get_q() << ", "
				<< mcc.get_pc() << ", "
				<< mcc.get_e_by_strain() << ", "
				<< mcc.get_e_by_model() << ", "
				<< mcc.get_f() << ", " << res << ",\n";
		//}
	}

	res_file.close();
}
