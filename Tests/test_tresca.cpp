#include "Tests_pcp.h"

#include <fstream>
#include <Eigen/Dense>
#include "Tresca.h"
#include "test_material_models.h"

namespace
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
		Consolidation = 2,
		Specified = 3
	};
	
	void output_var(
		std::ostream &os,
		double cohesion,
		const double out_strain[3],
		const double out_stress[6],
		const double out_pstrain[3]
		)
	{
		Eigen::Matrix3d s_mat;
		s_mat << out_stress[0], out_stress[3], out_stress[5],
				 out_stress[3], out_stress[1], out_stress[4],
				 out_stress[5], out_stress[4], out_stress[2];
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(s_mat);
		const Eigen::Vector3d& pstress = eigen_solver.eigenvalues();
		double smax = pstress[0];
		if (smax < pstress[1])
			smax = pstress[1];
		if (smax < pstress[2])
			smax = pstress[2];
		double smin = pstress[0];
		if (smin > pstress[1])
			smin = pstress[1];
		if (smin > pstress[2])
			smin = pstress[2];
		os << out_strain[0] << ", "
			<< out_strain[1] << ", "
			<< out_strain[2] << ", "
			<< out_stress[0] << ", "
			<< out_stress[1] << ", "
			<< out_stress[2] << ", "
			<< out_stress[3] << ", "
			<< out_stress[4] << ", "
			<< out_stress[5] << ", "
			<< out_pstrain[0] << ", "
			<< out_pstrain[1] << ", "
			<< out_pstrain[2] << ", "
			<< smax - smin - 2.0 * cohesion << "\n";
	}
}

void test_tresca()
{
	double de11, de22, de33;
	AnalysisType tp = AnalysisType::Specified;
	//AnalysisType tp = AnalysisType::TriaxialDrained;
	// tresca 1
	de11 = -0.05;
	de22 = 0.0;
	de33 = 0.0;
	// tresca 2
	//de11 = 0.05;
	//de22 = 0.0;
	//de33 = 0.0;
	// tresca 3
	//de11 = 0.05;
	//de22 = 0.0;
	//de33 = -0.05;
	size_t inc_num = 5000;
	size_t out_num = 100;

	double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	MatModel::Tresca tc;
	tc.set_param(1000.0, 0.1, 1.0, ini_stress);

	std::fstream res_file;
	res_file.open("Tresca_res.csv", std::ios::out | std::ios::binary);
	res_file << "e11, e22, e33, s11, s22, s33, s12, s23, s31,"
				"dep11, dep22, dep33, f\n";

	de11 /= double(inc_num);
	de22 /= double(inc_num);
	de33 /= double(inc_num);
	size_t out_inv = inc_num / out_num;
	double out_strain[3] = { 0.0, 0.0, 0.0 };
	const double(*Dep_mat)[6];
	double dstrain[6];
	for (size_t i = 0; i < inc_num; ++i)
	{
		if (i % out_inv == 0) // output
			output_var(res_file, tc.get_cohesion(), out_strain, tc.get_stress(), tc.get_dstrain_p());

		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			Dep_mat = reinterpret_cast<const double(*)[6]>(tc.get_Dep_mat());
			dstrain[0] = de11;
			dstrain[1] = -Dep_mat[1][0] / (Dep_mat[1][1] + Dep_mat[1][2]) * de11;
			dstrain[2] = -Dep_mat[2][0] / (Dep_mat[2][1] + Dep_mat[2][2]) * de11;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::TriaxialUndrained:
			dstrain[0] = de11;
			dstrain[1] = -0.5 * de11;
			dstrain[2] = -0.5 * de11;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::Consolidation:
			dstrain[0] = de11;
			dstrain[1] = 0.0;
			dstrain[2] = 0.0;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::Specified:
			dstrain[0] = de11;
			dstrain[1] = de22;
			dstrain[2] = de33;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		default:
			res_file.close();
			return;
		}

		int res = tc.integrate(dstrain);

		out_strain[0] += dstrain[0];
		out_strain[1] += dstrain[1];
		out_strain[2] += dstrain[2];
	}

	output_var(res_file, tc.get_cohesion(), out_strain, tc.get_stress(), tc.get_dstrain_p());
	res_file.close();
}
