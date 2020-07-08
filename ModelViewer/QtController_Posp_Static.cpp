#include "ModelViewer_pcp.h"

#include "QtController_Posp_Static.h"

QtController_Posp_Static
	::QtController_Posp_Static() :
	scene(nullptr), frame_id(0)
{

}

QtController_Posp_Static::QtController_Posp_Static(
	QtGLView& v,
	QtSceneFromHdf5_2DMPM& s) :
	scene(nullptr), frame_id(0)
{
	set_view(v);
	set_scene(s);
}

QtController_Posp_Static::~QtController_Posp_Static()
{

}

int QtController_Posp_Static::initialize(int wd, int ht)
{
	return scene->init_scene(wd, ht, frame_id);
}

void QtController_Posp_Static::draw_scene()
{
	scene->draw();
}

void QtController_Posp_Static::resize_scene(int wd, int ht)
{
	scene->resize(wd, ht);
}
