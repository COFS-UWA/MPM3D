#include "ModelViewer_pcp.h"

#include "QtController_Posp_Static.h"

QtController_Posp_Static
	::QtController_Posp_Static() :
	view(nullptr), scene(nullptr)
{

}

QtController_Posp_Static::QtController_Posp_Static(
	QtGLView& v,
	QtSceneFromHdf5_2DMPM& s
	)
{
	set_view(v);
	set_scene(s);
}

QtController_Posp_Static::~QtController_Posp_Static()
{

}

int QtController_Posp_Static::initialize(int wd, int ht)
{
	return scene->initialize(wd, ht);
}

void QtController_Posp_Static::draw_scene()
{
	scene->draw();
}

void QtController_Posp_Static::resize_scene(int wd, int ht)
{
	scene->resize(wd, ht);
}
