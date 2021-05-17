#include "TestsParallel_pcp.h"

#include <iostream>
#include <fstream>

#include "MatModelTestUtils.h"
#include "SandHypoplasticityWrapper.h"
#include "test_simulations_omp.h"

template <typename MatModel>
bool integrate_drained_triaxial_test_secant2(
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
	dstrain[0] = 0.3 * de;
	dstrain[1] = 0.3 * de;
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

void test_sand_hypoplasticity_wrapper()
{
	enum class AnalysisType : unsigned char
	{
		TriaxialDrained = 0,
		TriaxialUndrained = 1,
		Consolidation = 2
	};

	//AnalysisType tp = AnalysisType::TriaxialDrained;
	//double de = -0.2;
	//const double ini_stress[6] = { -60.0e3, -60.0e3, -120.0e3, 0.0, 0.0, 0.0 };
	//constexpr size_t inc_num = 20000;
	//constexpr size_t out_num = 200;
	
	AnalysisType tp = AnalysisType::Consolidation;
	double de = -0.03;
	const double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	constexpr size_t inc_num = 3000;
	constexpr size_t out_num = 100;

	const double ini_sr = ini_stress[0]; // radial stress
	MatModel::SandHypoplasticityWrapper shp;
	shp.set_param(
		ini_stress, 0.65,
		30.0, 1354.0e6, 0.34,
		0.49, 0.76, 0.86,
		0.18, 1.27);

	std::fstream res_file;
	res_file.open("res_sand_hypo_wrapper.csv", std::ios::out | std::ios::binary);
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
			bres = integrate_drained_triaxial_test_secant2<MatModel::SandHypoplasticityWrapper>(shp, de, ini_sr, dstrain);
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
