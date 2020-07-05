#include "ModelViewer_pcp.h"

#include "QtController.h"
#include "QtGLView.h"

QtGLView::QtGLView(QWidget* parent) :
	QOpenGLWidget(parent),
	controller(nullptr),
	fully_loaded(false), resize_num(0),
	draw_func(&QtGLView::draw_view_ini),
	resize_func(&QtGLView::resize_view_ini)
{

}

QtGLView::~QtGLView()
{

}

void QtGLView::initializeGL()
{
	initializeOpenGLFunctions();
	controller->initialize(width(), height());
}

void QtGLView::paintGL()
{
	(this->*draw_func)();
}

void QtGLView::resizeGL(int wd, int ht)
{
	(this->*resize_func)(wd, ht);
}

void QtGLView::draw_view_ini()
{
	if (resize_num == 2)
	{
		fully_loaded = true;
		draw_func = &QtGLView::draw_view;
		resize_func = &QtGLView::resize_view;
		update();
	}
}

void QtGLView::resize_view_ini(int wd, int ht)
{
	++resize_num;
}

void QtGLView::draw_view()
{
	controller->draw_scene();
}

void QtGLView::resize_view(int wd, int ht)
{
	controller->resize_scene(wd, ht);
}
