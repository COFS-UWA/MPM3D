#ifndef __Qt_App_Posp_T2D_CHM_s_h__
#define __Qt_App_Posp_T2D_CHM_s_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "ResultFile_hdf5.h"
#include "QtGLView.h"
#include "QtSceneFromHdf5_T2D_CHM_s.h"
#include "QtController_Posp_Static.h"
#include "QtController_Posp_Animation.h"

class QtApp_Posp_T2D_CHM_s;

namespace QtApp_Posp_T2D_CHM_s_Internal
{

class MainWindow : public QMainWindow
{
	Q_OBJECT
protected:
	friend QtApp_Posp_T2D_CHM_s;

	QtGLView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();

	inline QtGLView& get_view() { return *model_view; }
	// fixed window size
	inline void set_fixed_size(int wd, int ht)
	{
		setFixedSize(wd, ht);
	}
};

};

class QtApp_Posp_T2D_CHM_s
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

	QtApp_Posp_T2D_CHM_s_Internal::MainWindow window;
	QtSceneFromHdf5_T2D_CHM_s scene;
	QtController *pcontroller;

public:
	QtApp_Posp_T2D_CHM_s(int &argc, char **argv, Type tp = SingleFrame);
	~QtApp_Posp_T2D_CHM_s();
	
	inline QtGLView& get_view() { return window.get_view(); }
	
	inline void set_win_size(int wd, int ht)
	{
		if (type == Animation)
			window.set_fixed_size(wd, ht);
		else
			window.resize(wd, ht);
	}
	inline void set_display_range(double xl, double xu, double yl, double yu)
	{ scene.set_display_range(xl, xu, yl, yu); }

	inline void set_display_bg_mesh(bool op = true)
	{ scene.set_display_bg_mesh(op); }
	inline void set_display_pcls(bool op = true)
	{ scene.set_display_pcls(op); }

	inline void set_fld_range(double min, double max)
	{
		scene.set_fld_range(min, max);
	}

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
	inline int set_res_file(ResultFile_hdf5 &rf,
		const char *th_na, size_t f_id, const char *fld_na)
	{
		if (type != SingleFrame)
			return -1;
		int res;
		if (res = scene.set_res_file(rf, th_na, fld_na))
			return res;
		QtController_Posp_Static &pc
			= *static_cast<QtController_Posp_Static *>(pcontroller);
		if (res = pc.set_frame_id(f_id))
			return res;
		return 0;
	}

// ================= Animation only =================
	inline int set_res_file(ResultFile_hdf5& rf,
		const char* th_na, const char* fld_na)
	{
		if (type != Animation)
			return -1;
		int res;
		if (res = scene.set_res_file(rf, th_na, fld_na))
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
			= *static_cast<QtController_Posp_Animation *>(pcontroller);
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

#endif