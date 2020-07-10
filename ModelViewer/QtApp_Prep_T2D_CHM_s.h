#ifndef __Qt_App_Prep_T2D_CHM_s_h__
#define __Qt_App_Prep_T2D_CHM_s_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "QtGLView.h"
#include "QtSceneFromModel_T2D_CHM_s.h"
#include "QtController_Prep.h"

class QtApp_Prep_T2D_CHM_s;

namespace QtApp_Prep_T2D_CHM_s_Internal
{

class MainWindow : public QMainWindow
{
	Q_OBJECT
protected:
	friend QtApp_Prep_T2D_CHM_s;

	QtGLView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();
	inline QtGLView& get_view() { return *model_view; }
};

};


class QtApp_Prep_T2D_CHM_s
{
protected:
	QApplication app;

	QtApp_Prep_T2D_CHM_s_Internal::MainWindow window;
	QtSceneFromModel_T2D_CHM_s scene;
	QtController_Prep controller;

public:
	QtApp_Prep_T2D_CHM_s(int &argc, char **argv);
	~QtApp_Prep_T2D_CHM_s();
	
	inline QtGLView& get_view() { return window.get_view(); }
	
	inline void set_win_size(int wd, int ht) { window.resize(wd, ht); }
	inline void set_display_range(double xl, double xu, double yl, double yu)
	{ scene.set_display_range(xl, xu, yl, yu); }

	inline void set_display_bg_mesh(bool op = true)
	{ scene.set_display_bg_mesh(op); }
	inline void set_display_pcls(bool op = true)
	{ scene.set_display_pcls(op); }
	inline void set_display_pts(bool op = true)
	{ scene.set_display_pts(op); }

	inline void set_model(Model_T2D_CHM_s& model) { scene.set_model(model); }
	
	inline int set_pts_from_pcl_id(size_t* ids, size_t id_num, float radius)
	{ return scene.set_pts_from_pcl_id(ids, id_num, radius); }
	inline int set_pts_from_node_id(size_t* ids, size_t id_num, float radius)
	{ return scene.set_pts_from_node_id(ids, id_num, radius); }
	template <typename Point2D>
	inline int set_pts(Point2D *pts, size_t pt_num, float pt_r)
	{ return scene.set_pts(pts, pt_num, pt_r); }

	int start();
};

#endif