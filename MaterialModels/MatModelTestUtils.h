#ifndef __Mat_Model_Tests_Utils_h__
#define __Mat_Model_Tests_Utils_h__

enum class AnalysisType : uint16_t
{
	Invalid = 0,
	TriaxialDrained = 1,
	TriaxialUndrained = 2,
	OneDCompression = 3,
	IsoCompression = 4
};

template <typename MatModel, typename FloatType>
bool integrate_drained_triaxial_test_secant(
	MatModel& model,
	FloatType de, // axial strain
	FloatType sr, // radix stress
	FloatType dstrain[6], // output strain increment
	FloatType tol = 1.0e-3,
	size_t max_iter_num = 100)
{
	const FloatType de_tol = FloatType(abs(de)) * tol;
	const FloatType sr_tol = FloatType(abs(sr)) * tol;
	FloatType er0, er1, er2, sr0, sr1, sr2;
	MatModel tmp_md;
	const FloatType* stress = tmp_md.get_stress();
	// 0
	dstrain[0] = 0.0;
	dstrain[1] = 0.0;
	dstrain[2] = de;
	dstrain[3] = 0.0;
	dstrain[4] = 0.0;
	dstrain[5] = 0.0;
	tmp_md = model;
	tmp_md.integrate(dstrain);
	er1 = FloatType(0.0);
	sr1 = stress[0];
	// 1
	//const FloatType(*D)[6] = (const FloatType(*)[6])model.get_Dep_mat();
	//dstrain[0] = (D[1][2] * D[0][1] - D[1][1] * D[0][2]) * de / (D[0][0] * D[1][1] - D[1][0] * D[0][1]);
	//dstrain[1] = (D[0][2] * D[1][0] - D[0][0] * D[1][2]) * de / (D[0][0] * D[1][1] - D[1][0] * D[0][1]);
	dstrain[0] = -FloatType(0.1) * de;
	dstrain[1] = -FloatType(0.1) * de;
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

#endif