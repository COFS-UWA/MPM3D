#ifndef __Qt_App_Prep_T2D_CHM_mt_h__
#define __Qt_App_Prep_T2D_CHM_mt_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "QtGLView.h"
#include "QtSceneFromModel_T2D_CHM_mt.h"
#include "QtController_Prep.h"

class QtApp_Prep_T2D_CHM_mt;

namespace QtApp_Prep_T2D_CHM_mt_Internal
{

class MainWindow : public QMainWindow
{
	Q_OBJECT
protected:
	friend QtApp_Prep_T2D_CHM_mt;

	QtGLView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();
	inline QtGLView& get_view() { return *model_view; }
};

};


class QtApp_Prep_T2D_CHM_mt
{
protected:
	QApplication app;

	QtApp_Prep_T2D_CHM_mt_Internal::MainWindow window;
	QtSceneFromModel_T2D_CHM_mt scene;
	QtController_Prep controller;

public:
	QtApp_Prep_T2D_CHM_mt(int &argc, char **argv);
	~QtApp_Prep_T2D_CHM_mt();
	
	inline QtGLView& get_view() { return window.get_view(); }
	
	inline void set_win_size(int wd, int ht) { window.resize(wd, ht); }

	inline void set_display_whole_model()
	{ scene.set_display_whole_model(); }
	inline void set_display_range(double xl, double xu,
								  double yl, double yu)
	{ scene.set_display_range(xl, xu, yl, yu); }

	inline void set_display_bg_mesh(bool op = true)
	{ scene.set_display_bg_mesh(op); }
	inline void set_display_pcls(bool op = true)
	{ scene.set_display_pcls(op); }
	inline void set_display_rc(bool op = true)
	{ scene.set_display_rc(op); }
	inline void set_display_pts(bool op = true)
	{ scene.set_display_pts(op); }

	inline void set_model(Model_T2D_CHM_mt& model) { scene.set_model(model); }
	
	inline int set_pts_from_pcl_id(size_t* ids, size_t id_num, float radius)
	{ return scene.set_pts_from_pcl_id(ids, id_num, radius); }
	inline int set_pts_from_node_id(size_t* ids, size_t id_num, float radius)
	{ return scene.set_pts_from_node_id(ids, id_num, radius); }
	template <typename Point2D>
	inline int set_pts(Point2D *pts, size_t pt_num, float pt_r)
	{ return scene.set_pts(pts, pt_num, pt_r); }
	inline int set_pts_from_vx_s_bc(float radius) { return scene.set_pts_from_vx_s_bc(radius); }
	inline int set_pts_from_vy_s_bc(float radius) { return scene.set_pts_from_vy_s_bc(radius); }
	inline int set_pts_from_vx_f_bc(float radius) { return scene.set_pts_from_vx_f_bc(radius); }
	inline int set_pts_from_vy_f_bc(float radius) { return scene.set_pts_from_vy_f_bc(radius); }

	int start();
};

#endif