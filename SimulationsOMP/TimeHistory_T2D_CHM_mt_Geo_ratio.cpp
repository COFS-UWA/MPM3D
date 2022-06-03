#include "SimulationsOMP_pcp.h"

#include "Model_T2D_CHM_mt.h"
#include "Step_T2D_CHM_mt_Geo.h"

#include "TimeHistory_T2D_CHM_mt_Geo_ratio.h"

void TimeHistory_T2D_CHM_mt_Geo_ratio::close()
{
	if (is_init)
	{
		res_file.close();
		is_init = false;
	}
}

int time_history_output_func_T2D_CHM_mt_geo_ratio_to_text_res_file(TimeHistory& _self)
{
	TimeHistory_T2D_CHM_mt_Geo_ratio& th = static_cast<TimeHistory_T2D_CHM_mt_Geo_ratio &>(_self);
	Step_T2D_CHM_mt_Geo& step = static_cast<Step_T2D_CHM_mt_Geo&>(th.get_step());
	Model_T2D_CHM_mt& md = static_cast<Model_T2D_CHM_mt&>(step.get_model());

	// substep index, time, f_ub, e_kin, f_ub_ratio, e_kin_ratio;
	th.res_file << step.get_total_substep_index() << ", "
		<< step.get_total_time() << ", "
		<< step.get_f_ub() << ", "
		<< step.get_e_kin() << ", "
		<< step.get_f_ub_ratio() << ", "
		<< step.get_e_kin_ratio() << "\n";

	return 0;
}
