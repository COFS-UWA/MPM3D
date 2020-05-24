#ifndef __Prep_MPM_3D_App_h__
#define __Prep_MPM_3D_App_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "MPM3DModelView.h"

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

public:
	PrepMPM3DApp(int &argc, char* argv[]);
	~PrepMPM3DApp();

	inline void set_display_bg_mesh(bool op = true) { main_win->model_view->set_display_bg_mesh(op); }
	inline void set_display_pcls(bool op = true) { main_win->model_view->set_display_pcls(op); }
	inline void set_display_points(bool op = true) { main_win->model_view->set_display_points(op); }

	inline void set_view_dir(float x, float y, float z) { main_win->model_view->set_view_dir(x, y, z); }
	inline void set_view_dir(float theta, float fai) { main_win->model_view->set_view_dir(theta, fai); }
	inline void set_view_dist_scale(float dist_sc) { main_win->model_view->set_view_dist_scale(dist_sc); }
	void set_view_dir_and_scale(float x, float y, float z, float dist_sc);
	void set_view_dir_and_scale(float theta, float fai, float dist_sc);
	MPM3DModelView &get_model_view() { return *(main_win->model_view); }

	int start();

	// wrapper
	template <typename TMesh>
	inline int init_bg_mesh(TMesh& mesh, QVector3D& _color)
	{
		return main_win->model_view->init_bg_mesh(mesh, _color);
	}

	template <typename Particle>
	inline int init_monocolor_pcl_data(Particle* pcls, size_t pcl_num, QVector3D& _color)
	{
		return main_win->model_view->init_monocolor_pcl_data(pcls, pcl_num, _color);
	}

	template <typename Particle>
	inline int update_monocolor_pcl_data(Particle* pcls)
	{
		return main_win->model_view->update_monocolor_pcl_data(pcls);
	}

	inline int init_point_data(Point3D* pts, size_t pt_num, GLfloat pt_size, QVector3D& _color)
	{
		return main_win->model_view->init_point_data(pts, pt_num, pt_size, _color);
	}
};

#endif