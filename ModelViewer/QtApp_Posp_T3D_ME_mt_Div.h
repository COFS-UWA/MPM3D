#ifndef __Qt_App_Posp_T3D_ME_mt_Div_h__
#define __Qt_App_Posp_T3D_ME_mt_Div_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "ResultFile_hdf5.h"
#include "QtGLView.h"
#include "QtSceneFromHdf5_T3D_ME_mt_Div.h"
#include "QtController_Posp_Static.h"
#include "QtController_Posp_Animation.h"

template <class DivisionSet>
class QtApp_Posp_T3D_ME_mt_Div;

namespace QtApp_Posp_T3D_ME_mt_Div_Internal
{
class MainWindow : public QMainWindow
{
	Q_OBJECT

protected:
	template <class DivisionSet>
	friend class QtApp_Posp_T3D_ME_mt_Div;

	QtGLView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();

	inline QtGLView& get_view() { return *model_view; }
	
	// fixed window size
	inline void set_fixed_size(int wd, int ht)
	{ setFixedSize(wd, ht); }
};

};

template <class DivisionSet = EmptyDivisionSet>
class QtApp_Posp_T3D_ME_mt_Div
{
public:
	enum Type
	{
		SingleFrame = 0,
		Animation = 1
	};
	
protected:
	QApplication app;

	Type type;

	QtApp_Posp_T3D_ME_mt_Div_Internal::MainWindow window;
	QtSceneFromHdf5_T3D_ME_mt_Div<DivisionSet> scene;
	QtController *pcontroller;

public:
	QtApp_Posp_T3D_ME_mt_Div(int &argc, char **argv, Type tp = SingleFrame);
	~QtApp_Posp_T3D_ME_mt_Div();
	
	inline QtGLView& get_view() { return window.get_view(); }
	inline DivisionSet& get_div_set() { return scene.get_div_set(); }

	inline void set_win_size(int wd, int ht)
	{
		if (type == Animation)
			window.set_fixed_size(wd, ht);
		else
			window.resize(wd, ht);
	}

	inline void set_view_dir(float x, float y, float z) { scene.set_view_dir(x, y, z); }
	inline void set_view_dir(float theta, float fai) { scene.set_view_dir(theta, fai); }
	inline void set_view_dist_scale(float dist_sc) { scene.set_view_dist_scale(dist_sc); }
	inline void move_view_pos(float dx, float dy, float dz) { scene.move_view_pos(dx, dy, dz); }

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

	inline void set_color_map_fld_range(double min, double max)
	{ scene.set_color_map_fld_range(min, max); }
	inline void set_color_map_geometry(float xpos, float ypos, float ht)
	{ scene.set_color_map_geometry(xpos, ypos, ht); }
	inline void set_mono_color_pcl(bool _op = true) noexcept { scene.set_mono_color_pcl(_op); }

	inline void set_png_name(const char* name)
	{
		if (type == SingleFrame)
		{
			QtController_Posp_Static& pc
				= *static_cast<QtController_Posp_Static*>(pcontroller);
			pc.set_png_name(name);
		}
		else if (type == Animation)
		{
			QtController_Posp_Animation& pc
				= *static_cast<QtController_Posp_Animation*>(pcontroller);
			pc.set_png_name(name);
		}
	}

	int start();

	// ================= SingleFrame only =================
	inline int set_res_file(ResultFile_hdf5& rf,
		const char* th_name, size_t fm_id, Hdf5Field::FieldType fld_type)
	{
		if (type != SingleFrame)
			return -1;
		int res;
		if (res = scene.set_res_file(rf, th_name, fld_type))
			return res;
		QtController_Posp_Static& pc
			= *static_cast<QtController_Posp_Static*>(pcontroller);
		if (res = pc.set_frame_id(fm_id))
			return res;
		return 0;
	}

	// ================= Animation only =================
	inline int set_res_file(ResultFile_hdf5& rf,
		const char* th_name, Hdf5Field::FieldType fld_type)
	{
		if (type != Animation)
			return -1;
		int res;
		if (res = scene.set_res_file(rf, th_name, fld_type))
			return res;
		return 0;
	}

	inline void set_ani_time(double ani_t)
	{
		if (type != Animation) return;
		QtController_Posp_Animation& pc
			= *static_cast<QtController_Posp_Animation*>(pcontroller);
		pc.set_ani_time(ani_t);
	}

	inline void set_start_frame(size_t f_id)
	{
		if (type != Animation) return;
		QtController_Posp_Animation& pc
			= *static_cast<QtController_Posp_Animation*>(pcontroller);
		pc.set_start_frame(f_id);
	}

	inline void set_end_frame(size_t f_id)
	{
		if (type != Animation) return;
		QtController_Posp_Animation& pc
			= *static_cast<QtController_Posp_Animation*>(pcontroller);
		pc.set_end_frame(f_id);
	}

	inline void set_gif_name(const char* name)
	{
		if (type != Animation) return;
		QtController_Posp_Animation& pc
			= *static_cast<QtController_Posp_Animation*>(pcontroller);
		pc.set_gif_name(name);
	}
};

template <class DivisionSet>
QtApp_Posp_T3D_ME_mt_Div<DivisionSet>
	::QtApp_Posp_T3D_ME_mt_Div(
		int& argc, char** argv, Type tp) :
	app(argc, argv), type(tp),
	scene(window.get_view()),
	pcontroller(nullptr)
{
	switch (tp)
	{
	case SingleFrame:
		pcontroller = new QtController_Posp_Static(window.get_view(), scene);
		break;
	case Animation:
		pcontroller = new QtController_Posp_Animation(window.get_view(), scene);
		break;
	default:
		break;
	}
}

template <class DivisionSet>
QtApp_Posp_T3D_ME_mt_Div<DivisionSet>
	::~QtApp_Posp_T3D_ME_mt_Div()
{
	if (pcontroller)
	{
		delete pcontroller;
		pcontroller = nullptr;
	}
}

template <class DivisionSet>
int QtApp_Posp_T3D_ME_mt_Div<DivisionSet>::start()
{
	window.show();
	return app.exec();
}

#endif