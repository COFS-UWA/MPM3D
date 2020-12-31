#include "ModelViewer_pcp.h"

#include "QtController_Prep.h"

QtController_Prep::QtController_Prep() :
	view(nullptr), scene(nullptr) {}

QtController_Prep::QtController_Prep(
	QtGLView& v, QtSceneFromModel& s)
{
	set_view(v);
	set_scene(s);
}

QtController_Prep::~QtController_Prep() {}

int QtController_Prep::initialize(int wd, int ht)
{
	return scene->initialize(wd, ht);
}

void QtController_Prep::draw_scene()
{
	scene->draw();
}

void QtController_Prep::resize_scene(int wd, int ht)
{
	scene->resize(wd, ht);
}
