#include "SimulationsOMP_pcp.h"

#include "Model_T3D_CHM_mt.h"
#include "Step_T3D_CHM_ud_mt_subiter.h"

#include "TimeHistory_T3D_CHM_ud_mt_subiter_ratio.h"

void TimeHistory_T3D_CHM_ud_mt_subiter_ratio::close()
{
	if (is_init)
	{
		res_file.close();
		is_init = false;
	}
}

int time_history_output_func_T3D_CHM_ud_mt_subiter_ratio_to_text_res_file(TimeHistory& _self)
{
	TimeHistory_T3D_CHM_ud_mt_subiter_ratio& th = static_cast<TimeHistory_T3D_CHM_ud_mt_subiter_ratio &>(_self);
	Step_T3D_CHM_ud_mt_subiter& step = static_cast<Step_T3D_CHM_ud_mt_subiter&>(th.get_step());
	Model_T3D_CHM_mt& md = static_cast<Model_T3D_CHM_mt&>(step.get_model());

	th.res_file << step.get_total_substep_index() << ", "
		<< step.get_total_time() << ", "
		<< step.get_cur_e_kin() << ", "
		<< step.get_max_e_kin() << ", "
		<< step.get_subiter_index() << "\n";

	return 0;
}
