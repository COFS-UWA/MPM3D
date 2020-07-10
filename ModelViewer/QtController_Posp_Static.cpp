#include "ModelViewer_pcp.h"

#include <iostream>

#include "QtController_Posp_Static.h"

QtController_Posp_Static
	::QtController_Posp_Static() :
	scene(nullptr)
{
	init_self();
}

QtController_Posp_Static::QtController_Posp_Static(
	QtGLView& v,
	QtSceneFromHdf5& s) :
	scene(nullptr)
{
	set_view(v);
	set_scene(s);
	init_self();
}

void QtController_Posp_Static::init_self()
{
	frame_id = 0;
	
	need_output_png = false;
	png_name = "";

	screen_shot_timer.setSingleShot(true);
	connect(&screen_shot_timer, SIGNAL(timeout()), this, SLOT(window_to_png()));
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

	// only take screenshot once
	static size_t draw_num = 0;
	if (draw_num == 0 && need_output_png)
	{
		++draw_num;
		screen_shot_timer.start(0);
	}
}

void QtController_Posp_Static::resize_scene(int wd, int ht)
{
	scene->resize(wd, ht);
}

void QtController_Posp_Static::window_to_png()
{
	if (need_output_png)
	{
		std::string fname = png_name + ".png";
		QPixmap screen_pixels = view->grab(view->geometry());
		screen_pixels.save(fname.c_str(), "png");
		std::cout << "window save to " << fname.c_str() << "\n";
	}
}

void QtController_Posp_Static::set_png_name(const char* fname)
{
	if (!fname || !strcmp(fname, ""))
	{
		need_output_png = false;
		png_name = "";
		return;
	}
	need_output_png = true;
	png_name = fname;
}
