#include "ModelViewer_pcp.h"

#include "TimeHistory_ModelView_T3D_ME_s.h"


int TimeHistory_ModelView_T3D_ME_s::initialize_model_view_data()
{
	int res;
	if (model)
	{
		// background mesh
		QVector3D gray(0.5f, 0.5f, 0.5f);
		if ((res = view->init_bg_mesh<Model_T3D_ME_s>(*model, gray)) != 0)
			return res;

		// particle data
		QVector3D orange(0.8235f, 0.5137f, 0.0314f);
		if ((res = view->init_monocolor_pcl_data<Model_T3D_ME_s::Particle>(model->get_pcls(), model->get_pcl_num(), orange)) != 0)
			return res;
	}

	if (points && point_num != 0)
	{
		QVector3D red(1.0f, 0.0f, 0.0f);
		if ((res = view->init_point_data(points, point_num, 10.0f, red)) != 0)
			return res;
	}

	return 0;
}


int time_history_output_func_model_view_t3d_me_s_to_model_view(TimeHistory& _self)
{
	TimeHistory_ModelView_T3D_ME_s &self = static_cast<TimeHistory_ModelView_T3D_ME_s&>(_self);

	int res;
	res = self.view->update_monocolor_pcl_data<Model_T3D_ME_s::Particle>(self.model->get_pcls());

	// update model (send message to)
	// emit update();

	return 0;
}
