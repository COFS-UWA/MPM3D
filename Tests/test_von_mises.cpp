#include "Tests_pcp.h"

#include <fstream>
#include "VonMises.h"
#include "test_material_models.h"

void test_von_mises()
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
	};

	AnalysisType tp = AnalysisType::TriaxialDrained;
	//AnalysisType tp = AnalysisType::TriaxialUndrained;
	double de = -0.1;
	size_t inc_num = 10000;

	std::fstream res_file;
	const double(*Dep_mat)[6];
	double dstrain[6];

	double ini_stress[6] = { -1000.0, -1000.0, -1000.0, 0.0, 0.0, 0.0 };
	//double ini_stress[6] = { -1033.34, -983.33, -983.33, 0.0, 0.0, 0.0 };
	MatModel::VonMises vm;
	vm.set_param(1000.0, 0.3, 50.0, ini_stress);
	
	double e11 = 0.0, e22 = 0.0, e33 = 0.0;
	res_file.open("von_mises_res.csv", std::ios::out | std::ios::binary);
	res_file << "e11, e22, e33, s11, s22, s33, s12, s23, s31, "
				"q, f, norm_f, res\n"
				"0.0, 0.0, 0.0, "
			 << ini_stress[0] << ", "
			 << ini_stress[1] << ", "
			 << ini_stress[2] << ", "
			 << ini_stress[3] << ", "
			 << ini_stress[4] << ", "
			 << ini_stress[5] << ", "
			 << vm.get_q() << ", "
			 << vm.get_f() << ", "
			 << vm.get_norm_f()
			 << ", 0\n";

	de /= double(inc_num);
	for (size_t i = 0; i < inc_num; ++i)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			Dep_mat = reinterpret_cast<const double(*)[6]>(vm.get_Dep_mat());
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
		default:
			res_file.close();
			return;
		}

		int res = vm.integrate(dstrain);
		e11 += dstrain[0];
		e22 += dstrain[1];
		e33 += dstrain[2];
		const double* stress = vm.get_stress();
		// output result
		res_file << e11 << ", "
				<< e22 << ", "
				<< e33 << ", "
				<< stress[0] << ", "
				<< stress[1] << ", "
				<< stress[2] << ", "
				<< stress[3] << ", "
				<< stress[4] << ", "
				<< stress[5] << ", "
				 << vm.get_q() << ", "
				 << vm.get_f() << ", "
				 << vm.get_norm_f() << ", "
				 << res << "\n";
	}

	res_file.close();
}
