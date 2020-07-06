#ifndef __Qt_Controller_Posp_Static_h__
#define __Qt_Controller_Posp_Static_h__

#include "QtSceneFromHdf5_2DMPM.h"
#include "QtGLView.h"
#include "QtController.h"

class QtController_Posp_Static : public QtController
{
protected:
	QtGLView *view;
	QtSceneFromHdf5_2DMPM *scene;

public:
	QtController_Posp_Static();
	QtController_Posp_Static(QtGLView &v, QtSceneFromHdf5_2DMPM &s);
	~QtController_Posp_Static();

	inline void set_scene(QtSceneFromHdf5_2DMPM& s) { scene = &s; }

	int initialize(int wd, int ht) override;
	void draw_scene() override;
	void resize_scene(int wd, int ht) override;
};

#endif