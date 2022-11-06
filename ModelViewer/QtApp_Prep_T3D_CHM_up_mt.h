#ifndef __Qt_App_Prep_T3D_CHM_up_mt_h__
#define __Qt_App_Prep_T3D_CHM_up_mt_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "QtGLView.h"
#include "QtSceneFromModel_T3D_CHM_up_mt.h"
#include "QtController_Prep.h"

class QtApp_Prep_T3D_CHM_up_mt;

namespace QtApp_Prep_T3D_CHM_up_mt_Internal
{

class MainWindow : public QMainWindow
{
	Q_OBJECT
protected:
	friend QtApp_Prep_T3D_CHM_up_mt;

	QtGLView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();
	inline QtGLView& get_view() { return *model_view; }
};

};


class QtApp_Prep_T3D_CHM_up_mt
{
protected:
	QApplication app;

	QtApp_Prep_T3D_CHM_up_mt_Internal::MainWindow window;
	QtSceneFromModel_T3D_CHM_up_mt scene;
	QtController_Prep controller;

public:
	QtApp_Prep_T3D_CHM_up_mt(int &argc, char **argv);
	~QtApp_Prep_T3D_CHM_up_mt();
	
	inline QtGLView& get_view() { return window.get_view(); }
	
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
	inline void set_model(Model_T3D_CHM_up_mt& model) { scene.set_model(model); }

	inline int set_pts_from_pcl_id(size_t* ids, size_t id_num, float radius)
	{
		return scene.set_pts_from_pcl_id(ids, id_num, radius);
	}
	inline int set_pts_from_node_id(size_t* ids, size_t id_num, float radius)
	{
		return scene.set_pts_from_node_id(ids, id_num, radius);
	}
	inline int set_pts_from_drained_bc(float radius) { return scene.set_pts_from_drained_bc(radius); }
	inline int set_pts_from_vx_bc(float radius) { return scene.set_pts_from_vx_bc(radius); }
	inline int set_pts_from_vy_bc(float radius) { return scene.set_pts_from_vy_bc(radius); }
	inline int set_pts_from_vz_bc(float radius) { return scene.set_pts_from_vz_bc(radius); }
	inline int set_pts_from_vec_bc(float radius) { return scene.set_pts_from_vec_bc(radius); }

	template <typename Point3D>
	inline int set_pts(Point3D* pts, size_t pt_num, float radius)
	{
		return scene.set_pts(pts, pt_num, radius);
	}

	int start();
};

#endif