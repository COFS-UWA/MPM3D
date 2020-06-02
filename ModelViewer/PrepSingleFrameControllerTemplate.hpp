#ifndef __Prep_Single_Frame_Controller_Template_h__
#define __Prep_Single_Frame_Controller_Template_h__

#include "PrepSingleFrameControllerBase.h"

template <typename Model>
class PrepSingleFrameControllerTemplate :
	public PrepSingleFrameControllerBase
{
protected:
	Model *model;

public:
	PrepSingleFrameControllerTemplate(MPM3DModelView& v) :
		PrepSingleFrameControllerBase(v), model(nullptr) {}
	~PrepSingleFrameControllerTemplate() { close(); }
	void close() {}

	inline void set_model(Model &md) { model = &md; }

protected:
	int initialize_model_view_data() override
	{
		typedef typename Model::Particle Particle;
		int res;
		if (model)
		{
			// background mesh
			QVector3D gray(0.5f, 0.5f, 0.5f);
			if ((res = view->init_bg_mesh<Model>(*model, gray)) != 0)
				return res;

			// particle data
			//QVector3D orange(0.8235f, 0.5137f, 0.0314f);
			//if ((res = view->init_monocolor_pcl_data<Particle>(model->get_pcls(), model->get_pcl_num(), orange)) != 0)
			//	return res;
			ValueToColor::Colorf orange(0.8235f, 0.5137f, 0.0314f);
			if ((res = view->init_monocolor_pcl_data<Particle>(model->get_pcls(), model->get_pcl_num(), orange)) != 0)
				return res;
		}

		if ((res = PrepSingleFrameControllerBase::initialize_model_view_data()) != 0)
			return res;
		
		return 0;
	}
};

#endif