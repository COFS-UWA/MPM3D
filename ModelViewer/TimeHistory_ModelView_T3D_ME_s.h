#ifndef __Time_History_Model_View_T3D_ME_s_h__
#define __Time_History_Model_View_T3D_ME_s_h__

#include "TimeHistory_ModelView.h"

#include "Model_T3D_ME_s.h"

int time_history_output_func_model_view_t3d_me_s_to_model_view(TimeHistory& _self);

// input data of Model_T3D_ME_s into MPM3DModelView
class TimeHistory_ModelView_T3D_ME_s : public TimeHistory_ModelView
{
protected:
	friend int time_history_output_func_model_view_t3d_me_s_to_model_view(TimeHistory& _self);

	Model_T3D_ME_s *model;
	
public:
	TimeHistory_ModelView_T3D_ME_s(const char* _name = "ModelView") :
		TimeHistory_ModelView(nullptr, 0, nullptr, _name, "TimeHistory_ModelView_T3D_ME_s",
			&time_history_output_func_model_view_t3d_me_s_to_model_view),
		model(nullptr) {}

	TimeHistory_ModelView_T3D_ME_s(Model_T3D_ME_s &md, PrepMPM3DApp &v,
		const char* _name = "ModelView") :
		TimeHistory_ModelView(nullptr, 0, &v.get_model_view(),
			_name, "TimeHistory_ModelView_T3D_ME_s",
			&time_history_output_func_model_view_t3d_me_s_to_model_view),
		model(&md) {}

	~TimeHistory_ModelView_T3D_ME_s() {}

	inline void set_model(Model_T3D_ME_s &md) { model = &md; }

	int initialize_model_view_data() override;
};

#endif