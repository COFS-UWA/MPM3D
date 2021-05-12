#ifndef __Qt_App_Prep_T3D_CHM_mt_Div_h__
#define __Qt_App_Prep_T3D_CHM_mt_Div_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "QtGLView.h"
#include "QtSceneFromModel_T3D_CHM_mt_Div.h"
#include "QtController_Prep.h"

template <class DivisionSet>
class QtApp_Prep_T3D_CHM_mt_Div;

namespace QtApp_Prep_T3D_CHM_mt_Div_Internal
{

class MainWindow : public QMainWindow
{
	Q_OBJECT
protected:
	template <class DivisionSet>
	friend class QtApp_Prep_T3D_CHM_mt_Div;

	QtGLView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();
	inline QtGLView& get_view() { return *model_view; }
};

};

template <class DivisionSet = EmptyDivisionSet>
class QtApp_Prep_T3D_CHM_mt_Div
{
protected:
	QApplication app;

	QtApp_Prep_T3D_CHM_mt_Div_Internal::MainWindow window;
	QtSceneFromModel_T3D_CHM_mt_Div<DivisionSet> scene;
	QtController_Prep controller;

public:
	QtApp_Prep_T3D_CHM_mt_Div(int &argc, char **argv);
	~QtApp_Prep_T3D_CHM_mt_Div();
	
	inline QtGLView& get_view() { return window.get_view(); }
	inline DivisionSet& get_div_set() { return scene.get_div_set(); }

	inline void set_win_size(int wd, int ht) { window.resize(wd, ht); }

	inline void set_view_dir(float x, float y, float z) { scene.set_view_dir(x, y, z); }
	inline void set_view_dir(float theta, float fai) { scene.set_view_dir(theta, fai); }
	inline void set_view_dist_scale(float dist_sc) { scene.set_view_dist_scale(dist_sc); }

	inline void set_light_dir(float x, float y, float z) { scene.set_light_dir(x, y, z); }
	inline void set_light_dir(float theta, float fai) { scene.set_light_dir(theta, fai); }
	inline void set_light_dist_scale(float dist_sc) { scene.set_light_dist_scale(dist_sc); }

	inline void set_bg_color(QVector3D& color) { scene.set_bg_color(color); }

	inline void set_fog_coef(float coef) { scene.set_fog_coef(coef); }
	inline void set_fog_color(QVector3D& color) { scene.set_fog_color(color); }

	inline void set_light_color(QVector3D& color) { scene.set_light_color(color); }
	inline void set_amb_coef(float coef) { scene.set_amb_coef(coef); }
	inline void set_diff_coef(float coef) { scene.set_diff_coef(coef); }
	inline void set_spec_coef(float coef) { scene.set_spec_coef(coef); }
	inline void set_spec_shininess(float shininess) { scene.set_spec_shininess(shininess); }

	inline void set_display_bg_mesh(bool op = true) { scene.set_display_bg_mesh(op); }
	inline void set_display_pcls(bool op = true) { scene.set_display_pcls(op); }
	inline void set_display_pts(bool op = true) { scene.set_display_pts(op); }
	inline void set_model(Model_T3D_CHM_mt& model) { scene.set_model(model); }

	inline int set_pts_from_pcl_id(size_t* ids, size_t id_num, float radius)
	{ return scene.set_pts_from_pcl_id(ids, id_num, radius); }
	inline int set_pts_from_node_id(size_t* ids, size_t id_num, float radius)
	{ return scene.set_pts_from_node_id(ids, id_num, radius); }
	inline int set_pts_from_vx_s_bc(float radius) { return scene.set_pts_from_vx_s_bc(radius); }
	inline int set_pts_from_vy_s_bc(float radius) { return scene.set_pts_from_vy_s_bc(radius); }
	inline int set_pts_from_vz_s_bc(float radius) { return scene.set_pts_from_vz_s_bc(radius); }
	inline int set_pts_from_vx_f_bc(float radius) { return scene.set_pts_from_vx_f_bc(radius); }
	inline int set_pts_from_vy_f_bc(float radius) { return scene.set_pts_from_vy_f_bc(radius); }
	inline int set_pts_from_vz_f_bc(float radius) { return scene.set_pts_from_vz_f_bc(radius); }
	template <typename Point3D>
	inline int set_pts(Point3D* pts, size_t pt_num, float radius)
	{ return scene.set_pts(pts, pt_num, radius); }

	int start();
};

template <class DivisionSet>
QtApp_Prep_T3D_CHM_mt_Div<DivisionSet>::QtApp_Prep_T3D_CHM_mt_Div(int& argc, char** argv) :
	app(argc, argv),
	scene(window.get_view()),
	controller(window.get_view(), scene) {}

template <class DivisionSet>
QtApp_Prep_T3D_CHM_mt_Div<DivisionSet>::~QtApp_Prep_T3D_CHM_mt_Div() {}

template <class DivisionSet>
int QtApp_Prep_T3D_CHM_mt_Div<DivisionSet>::start()
{
	window.show();
	return app.exec();
}

#endif