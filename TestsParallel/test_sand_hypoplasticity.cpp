#include "TestsParallel_pcp.h"

#include <fstream>

#include "MatModelUtils.h"
#include "SandHypoplasticityByUmat.h"
#include "test_simulations_omp.h"

template <typename MatModel>
bool integrate_drained_triaxial_test_secant(
	MatModel& model,
	double de, // axial strain
	double sr, // radix stress
	double dstrain[6], // output strain increment
	double tol = 1.0e-3,
	size_t max_iter_num = 100)
{
	const double de_tol = abs(de) * tol;
	const double sr_tol = abs(sr) * tol;
	MatModel tmp_md;
	const double* stress = tmp_md.get_stress();
	const double(*D)[6] = (const double(*)[6])model.get_Dep_mat();
	double er0, er1, er2, sr0, sr1, sr2;
	dstrain[2] = de;
	dstrain[3] = 0.0;
	dstrain[4] = 0.0;
	dstrain[5] = 0.0;
	// 0
	dstrain[0] = 0.0;
	dstrain[1] = 0.0;
	dstrain[2] = de;
	dstrain[3] = 0.0;
	dstrain[4] = 0.0;
	dstrain[5] = 0.0;
	tmp_md = model;
	tmp_md.integrate(dstrain);
	er1 = 0.0;
	sr1 = stress[0];
	// 1
	dstrain[0] = (D[1][2] * D[0][1] - D[1][1] * D[0][2]) * de / (D[0][0] * D[1][1] - D[1][0] * D[0][1]);
	dstrain[1] = (D[0][2] * D[1][0] - D[0][0] * D[1][2]) * de / (D[0][0] * D[1][1] - D[1][0] * D[0][1]);
	dstrain[2] = de;
	dstrain[3] = 0.0;
	dstrain[4] = 0.0;
	dstrain[5] = 0.0;
	tmp_md = model;
	tmp_md.integrate(dstrain);
	er2 = dstrain[0];
	sr2 = stress[0];
	size_t iter_id = 0;
	do
	{
		er0 = er1;
		er1 = er2;
		sr0 = sr1;
		sr1 = sr2;
		// sn+1
		tmp_md = model;
		er2 = er1 + (sr - sr1) * (er1 - er0) / (sr1 - sr0);
		dstrain[0] = er2;
		dstrain[1] = er2;
		dstrain[2] = de;
		dstrain[3] = 0.0;
		dstrain[4] = 0.0;
		dstrain[5] = 0.0;
		tmp_md.integrate(dstrain);
		sr2 = stress[0];
		if (abs(er2 - er1) <= de_tol && abs(sr2 - sr1) <= sr_tol)
		{
			model.integrate(dstrain);
			return true;
		}
	} while ((++iter_id) < max_iter_num);
	return false;
}

namespace
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
		Consolidation = 2
	};
}

void test_sand_hypoplasticity_integration()
{
	double dstrain[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	const double ini_stress[6] = {
		-100.0e3,
		-100.0e3,
		-100.0e3,
		0.0, 0.0, 0.0
	};
	constexpr double R = 1.0e-4;
	const double ig_strain[6] = {
		-R / sqrt(3.0),
		-R / sqrt(3.0),
		-R / sqrt(3.0),
		0.0, 0.0, 0.0
	};

	int res;
	MatModel::SandHypoplasticityByUmat shp;
	shp.set_param(ini_stress, 0.817,
		33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
		2.2, 1.1, R, 0.1, 5.5, ig_strain);

	const double *stress = shp.get_stress();
	std::cout << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
			  << stress[3] << ", " << stress[4] << ", " << stress[5] << "\n";

	dstrain[0] = 5.13e-4;
	dstrain[1] = 5.13e-4;
	//dstrain[0] = 0.000275933;
	//dstrain[1] = 0.000275933;
	dstrain[2] = -2.50e-3;
	res = shp.integrate(dstrain);
	std::cout << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
			  << stress[3] << ", " << stress[4] << ", " << stress[5] << "\n";

	dstrain[0] = 7.67e-4;
	dstrain[1] = 7.67e-4;
	res = shp.integrate(dstrain);
	std::cout << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
		<< stress[3] << ", " << stress[4] << ", " << stress[5] << "\n";

}

void test_sand_hypoplasticity_Herleand_Gudehus_1999()
{

}

void test_sand_hypoplasticity_Wichtmann_2019()
{
	AnalysisType tp = AnalysisType::TriaxialDrained;
	double de = -0.25;
	// need high inc_num to be accurate
	constexpr size_t inc_num = 250000;
	constexpr size_t out_num = 100;
	const double ini_stress[6] = { -100.0e3, -100.0e3, -100.0e3, 0.0, 0.0, 0.0 };
	constexpr double R = 1.0e-4;
	const double ig_strain[6] = { -R/sqrt(3.0), -R/sqrt(3.0), -R/sqrt(3.0), 0.0, 0.0, 0.0 };

	MatModel::SandHypoplasticityByUmat shp;
	shp.set_param(ini_stress, 0.817,
		33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
		2.2, 1.1, R, 0.1, 5.5, ig_strain);
	std::fstream res_file;
	res_file.open("res_Wichtmann.csv", std::ios::out | std::ios::binary);
	const double* stress = shp.get_stress();
	double e11 = 0.0;
	double e22 = 0.0;
	double e33 = 0.0;
	res_file << " e11, e22, e33, s11, s22, s33, q, p, e\n"
		<< e11 << ", "
		<< e22 << ", "
		<< e33 << ", "
		<< stress[0] << ", "
		<< stress[1] << ", "
		<< stress[2] << ", "
		<< stress[0] - stress[2] << ", "
		<< (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
		<< shp.get_e() << "\n";

	de = de / double(inc_num);
	const size_t out_inv = inc_num / out_num;
	double dstrain[6];
	bool bres;
	int res;
	for (size_t i = 0; i < inc_num; ++i)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			bres = MatModel::integrate_drained_triaxial_test(shp, de, dstrain);
			break;
		case AnalysisType::TriaxialUndrained:
			dstrain[0] = -0.5 * de;
			dstrain[1] = -0.5 * de;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			res = shp.integrate(dstrain);
			break;
		case AnalysisType::Consolidation:
			dstrain[0] = 0.0;
			dstrain[1] = 0.0;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			res = shp.integrate(dstrain);
			break;
		default:
			res_file.close();
			return;
		}

		e11 += dstrain[0];
		e22 += dstrain[1];
		e33 += dstrain[2];

		// output result
		if (i % out_inv == out_inv - 1)
		{
			std::cout << "Output at step " << i << ".\n";
			res_file << e11 << ", "
				<< e22 << ", "
				<< e33 << ", "
				<< stress[0] << ", "
				<< stress[1] << ", "
				<< stress[2] << ", "
				<< stress[0] - stress[2] << ", "
				<< (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
				<< shp.get_e() << "\n";
		}
	}

	res_file.close();
}

void test_triaxial_secant()
{
	double de = -1.0e-6;
	const double ini_stress[6] = { -100.0, -100.0, -100.0, 0.0, 0.0, 0.0 };
	constexpr double ini_sr = -100.0;
	constexpr double R = 1.0e-4;
	const double ig_strain[6] = { -R / sqrt(3.0), -R / sqrt(3.0), -R / sqrt(3.0), 0.0, 0.0, 0.0 };
	MatModel::SandHypoplasticityByUmat shp;
	const double* stress = shp.get_stress();
	// intergranular strain
	shp.set_param(ini_stress, 0.817,
		33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
		2.2, 1.1, R, 0.1, 5.5, ig_strain);
	// no intergranular strain
	//shp.set_param(ini_stress, 0.817,
	//	33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
	//	0.0, 1.1, R, 0.1, 5.5, ig_strain);
	double dstrain[6];
	bool bres = integrate_drained_triaxial_test_secant<MatModel::SandHypoplasticityByUmat>(shp, de, ini_sr, dstrain);
	double s11 = stress[0];
	double s22 = stress[1];
	double s33 = stress[2];
	size_t efe = 0;
}

void test_sand_hypoplasticity_triaxial()
{
	AnalysisType tp = AnalysisType::TriaxialDrained;
	//AnalysisType tp = AnalysisType::TriaxialUndrained;
	double de = -0.06;
	// need high inc_num to be accurate
	constexpr size_t inc_num = 60000;
	constexpr size_t out_num = 120;
	//const double ini_stress[6] = { -100.0e3, -100.0e3, -100.0e3, 0.0, 0.0, 0.0 };
	const double ini_stress[6] = { -100.0, -100.0, -100.0, 0.0, 0.0, 0.0 };
	const double ini_sr = -100.0;
	constexpr double R = 1.0e-4;
	const double ig_strain[6] = { -R / sqrt(3.0), -R / sqrt(3.0), -R / sqrt(3.0), 0.0, 0.0, 0.0 };
	MatModel::SandHypoplasticityByUmat shp;
	// intergranular strain
	shp.set_param(ini_stress, 0.817,
		33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
		2.2, 1.1, R, 0.1, 5.5, ig_strain);
	// no intergranular strain
	//shp.set_param(ini_stress, 0.817,
	//	33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
	//	0.0, 1.1, R, 0.1, 5.5, ig_strain);
	std::fstream res_file;
	res_file.open("res_triaxial.csv", std::ios::out | std::ios::binary);
	const double* stress = shp.get_stress();
	double e11 = 0.0;
	double e22 = 0.0;
	double e33 = 0.0;
	res_file << " e11, e22, e33, s11, s22, s33, q, p, e\n"
		<< e11 << ", "
		<< e22 << ", "
		<< e33 << ", "
		<< stress[0] << ", "
		<< stress[1] << ", "
		<< stress[2] << ", "
		<< stress[0] - stress[2] << ", "
		<< (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
		<< shp.get_e() << "\n";

	de = de / double(inc_num);
	const size_t out_inv = inc_num / out_num;
	double dstrain[6];
	bool bres;
	int res;
	for (size_t i = 0; i < inc_num; ++i)
	{
		switch (tp)
		{
		case AnalysisType::TriaxialDrained:
			//bres = MatModel::integrate_drained_triaxial_test(shp, de, dstrain);
			bres = integrate_drained_triaxial_test_secant<MatModel::SandHypoplasticityByUmat>(shp, de, ini_sr, dstrain);
			break;
		case AnalysisType::TriaxialUndrained:
			dstrain[0] = -0.5 * de;
			dstrain[1] = -0.5 * de;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			res = shp.integrate(dstrain);
			break;
		case AnalysisType::Consolidation:
			dstrain[0] = 0.0;
			dstrain[1] = 0.0;
			dstrain[2] = de;
			dstrain[3] = 0.0;
			dstrain[4] = 0.0;
			dstrain[5] = 0.0;
			res = shp.integrate(dstrain);
			break;
		default:
			res_file.close();
			return;
		}

		e11 += dstrain[0];
		e22 += dstrain[1];
		e33 += dstrain[2];

		// output result
		if (i % out_inv == out_inv - 1)
		{
			std::cout << "Output at step " << i + 1 << ".\n";
			res_file << e11 << ", " << e22 << ", " << e33 << ", "
				<< stress[0] << ", "
				<< stress[1] << ", "
				<< stress[2] << ", "
				<< stress[0] - stress[2] << ", "
				<< (stress[0] + stress[1] + stress[2]) / 3.0 << ", "
				<< shp.get_e() << "\n";
		}
	}

	res_file.close();
}
