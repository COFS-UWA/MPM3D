#include "TestsParallel_pcp.h"

#include <fstream>

#include "MatModelTestUtils.h"
#include "MohrCoulombWrapper.h"
#include "test_simulations_omp.h"

void test_mohr_coulomb()
{
	const double ini_s[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	MatModel::MohrCoulombWrapper mc;
	mc.set_param(ini_s, 10.0, 0.0, 1.0, 1000.0, 0.0);

	double dstrain[6];

	// elastic
	//dstrain[0] = 0.0;
	//dstrain[1] = 0.0;
	//dstrain[2] = -1.0e-5;
	//dstrain[3] = 0.0;
	//dstrain[4] = 0.0;
	//dstrain[5] = 0.0;

	// area0
	//dstrain[0] = 2.0e-3;
	//dstrain[1] = 0.0;
	//dstrain[2] = -5.0e-3;
	//dstrain[3] = 0.0;
	//dstrain[4] = 0.0;
	//dstrain[5] = 0.0;

	// area1
	//dstrain[0] = 2.0e-3;
	//dstrain[1] = 2.0e-3;
	//dstrain[2] = -5.0e-3;
	//dstrain[3] = 0.0;
	//dstrain[4] = 0.0;
	//dstrain[5] = 0.0;

	// area1
	dstrain[0] = -5.0e-3;
	dstrain[1] = -5.0e-3;
	dstrain[2] = 2.0e-3;
	dstrain[3] = 0.0;
	dstrain[4] = 0.0;
	dstrain[5] = 0.0;

	//// area3
	//dstrain[0] = 3.0e-3;
	//dstrain[1] = 2.0e-3;
	//dstrain[2] = 5.0e-3;
	//dstrain[3] = 0.0;
	//dstrain[4] = 0.0;
	//dstrain[5] = 0.0;

	int32_t res = mc.integrate(dstrain);
	const double* stress = mc.get_stress();

	std::cout << stress[0] << ", "
		<< stress[1] << ", "
		<< stress[2] << ", "
		<< stress[3] << ", "
		<< stress[4] << ", "
		<< stress[5] << "\n";
}
