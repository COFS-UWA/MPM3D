#include "Tests_pcp.h"

#include <fstream>

#include "SandHypoplasticityByUmat.h"
#include "test_simulations.h"

void test_hypoplasticity_match_Wichtmann()
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
		Consolidation = 2
	};

	AnalysisType tp = AnalysisType::TriaxialDrained;
	double de = -0.25;
	constexpr size_t inc_num = 2500;
	const double ini_stress[6] = { -100.0e3, -100.0e3, -100.0e3, 0.0, 0.0, 0.0 };
	constexpr double e0 = 0.817;

	std::fstream res_file;
	const double(*Dep_mat)[6];
	double dstrain[6];
	const double *stress;

	de = de / double(inc_num);
	MatModel::SandHypoplasticityByUmat shp;
	constexpr double R = 1.0e-4;
	const double ig_strain[6] = {
		-R/sqrt(3.0),
		-R/sqrt(3.0),
		-R/sqrt(3.0),
		0.0, 0.0, 0.0
	};
	shp.set_param(ini_stress, e0,
		33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
		2.2, 1.1, 1.0e-4, 0.1, 5.5, ig_strain);
	res_file.open("res_Wichtmann.csv", std::ios::out | std::ios::binary);
	stress = shp.get_stress();
	double e11 = 0.0, e22 = 0.0, e33 = 0.0;
	res_file << " e11, e22, e33, s11, s22, s33, q, p, e\n"
		<< e11 << ", " << e22 << ", " << e33 << ", "
		<< stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
		<< stress[0] - stress[1] << ", "
		<< (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
		<< shp.get_e() << "\n";

	for (size_t i = 0; i < inc_num; ++i)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			Dep_mat = reinterpret_cast<const double(*)[6]>(shp.get_Dep_mat());
			dstrain[0] = -Dep_mat[1][0] / (Dep_mat[1][1] + Dep_mat[1][2]) * de;
			dstrain[1] = -Dep_mat[2][0] / (Dep_mat[2][1] + Dep_mat[2][2]) * de;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::TriaxialUndrained:
			dstrain[0] = -0.5 * de;
			dstrain[1] = -0.5 * de;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		case AnalysisType::Consolidation:
			dstrain[0] = 0.0;
			dstrain[1] = 0.0;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			break;
		default:
			res_file.close();
			return;
		}

		e11 += dstrain[0]; e22 += dstrain[1]; e33 += dstrain[2];
		
		int res = shp.integrate(dstrain);

		// output result
		stress = shp.get_stress();
		res_file << e11 << ", "
				 << e22 << ", "
				 << e33 << ", "
				 << stress[0] << ", "
				 << stress[1] << ", "
				 << stress[2] << ", "
				 <<  stress[0] - stress[1] << ", "
				 << (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
				 << shp.get_e() << "\n";
	}

	res_file.close();
}
