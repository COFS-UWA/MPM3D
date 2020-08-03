#include "MaterialModels_pcp.h"

#include "UndrainedModifiedCamClay.h"

namespace MatModel
{
	// return value:
	//     = 0 - elstic
	//     > 0 - plastic
	//     < 0 - convergence failure
	int undrained_modified_cam_clay_integration_function(MaterialModel* _self, double dstrain[6])
	{
		UndrainedModifiedCamClay& self = static_cast<UndrainedModifiedCamClay&>(*_self);
		ModifiedCamClay& mcc = self.mcc;

		int res;
		if ((res = mcc.integrate(dstrain)) < 0)
			return res;
		self.pore_pressure += self.Kw_div_n * -(dstrain[0] + dstrain[1] + dstrain[2]);

		self.Kw_div_n = self.Kw * (1.0 / mcc.get_e_by_model() + 1.0);
		self.cal_total_stress();
		self.cal_Kmat_with_Kw();

		return res;
	}
}