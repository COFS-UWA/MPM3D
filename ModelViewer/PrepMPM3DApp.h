#ifndef __Prep_MPM_3D_App_h__
#define __Prep_MPM_3D_App_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "MPM3DModelView.h"
#include "PrepSingleFrameControllerTemplate.hpp"

class PrepMPM3DApp;

namespace PreMPM3DApp_Internal
{

// main window of PreMPM3DApp
class MainWindow : public QMainWindow
{
	Q_OBJECT
protected:
	friend PrepMPM3DApp;

	MPM3DModelView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();
};

};


class PrepMPM3DApp
{
protected:
	QApplication app;

	PreMPM3DApp_Internal::MainWindow *main_win;
	MPM3DModelView::Controller *view_controller;

public:
	PrepMPM3DApp(int &argc, char **argv);
	~PrepMPM3DApp();

	MPM3DModelView& get_model_view() { return *(main_win->model_view); }
	
	inline void set_display_bg_mesh(bool op = true) { main_win->model_view->set_display_bg_mesh(op); }
	inline void set_display_pcls(bool op = true) { main_win->model_view->set_display_pcls(op); }
	inline void set_display_points(bool op = true) { main_win->model_view->set_display_points(op); }

	inline void set_view_dir(float x, float y, float z) { main_win->model_view->set_view_dir(x, y, z); }
	inline void set_view_dir(float theta, float fai) { main_win->model_view->set_view_dir(theta, fai); }
	inline void set_view_dist_scale(float dist_sc) { main_win->model_view->set_view_dist_scale(dist_sc); }

	template <typename Model>
	inline void set_model(Model &md)
	{
		typedef PrepSingleFrameControllerTemplate<Model> Controller;
		view_controller = new Controller(get_model_view());
		Controller &con = *static_cast<Controller *>(view_controller);
		con.set_model(md);
	}

	inline void set_points(Point3D *points, size_t point_num)
	{
		PrepSingleFrameControllerBase &con
			= *static_cast<PrepSingleFrameControllerBase *>(view_controller);
		if (points && point_num != 0)
			con.set_points(points, point_num);
	}

	int start();
};

#endif