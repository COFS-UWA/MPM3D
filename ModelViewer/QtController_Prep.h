#ifndef __Qt_Controller_Prep_h__
#define __Qt_Controller_Prep_h__

#include "QtSceneFromModel_2DMPM.h"
#include "QtGLView.h"
#include "QtController.h"

class QtController_Prep : public QtController
{
protected:
	QtGLView *view;
	QtSceneFromModel_2DMPM *scene;

public:
	QtController_Prep();
	QtController_Prep(QtGLView &v, QtSceneFromModel_2DMPM &s);
	~QtController_Prep();

	inline void set_scene(QtSceneFromModel_2DMPM& s) { scene = &s; }

	int initialize(int wd, int ht) override;
	void draw_scene() override;
	void resize_scene(int wd, int ht) override;
};

#endif