#ifndef __Qt_Controller_Prep_h__
#define __Qt_Controller_Prep_h__

#include "QtGLView.h"
#include "QtSceneFromModel.h"
#include "QtController.h"

class QtController_Prep : public QtController
{
protected:
	QtGLView *view;
	QtSceneFromModel *scene;

public:
	QtController_Prep();
	QtController_Prep(QtGLView &v, QtSceneFromModel &s);
	~QtController_Prep();

	inline void set_scene(QtSceneFromModel& s) { scene = &s; }

	int initialize(int wd, int ht) override;
	void draw_scene() override;
	void resize_scene(int wd, int ht) override;
};

#endif