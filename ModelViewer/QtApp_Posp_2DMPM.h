#ifndef __Qt_App_Posp_2DMPM_h__
#define __Qt_App_Posp_2DMPM_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "ResultFile_hdf5.h"
#include "QtGLView.h"
#include "QtSceneFromHdf5_2DMPM.h"
#include "QtController_Posp_Static.h"

class QtApp_Posp_2DMPM;

namespace QtApp_Posp_2DMPM_Internal
{

class MainWindow : public QMainWindow
{
	Q_OBJECT
protected:
	friend QtApp_Posp_2DMPM;

	QtGLView *model_view;
	QStatusBar *status_bar;

public:
	MainWindow();
	~MainWindow();
	inline QtGLView& get_view() { return *model_view; }
};

};

class QtApp_Posp_2DMPM
{
protected:
	QApplication app;

	QtApp_Posp_2DMPM_Internal::MainWindow window;
	QtSceneFromHdf5_2DMPM scene;
	QtController_Posp_Static controller;

public:
	QtApp_Posp_2DMPM(int &argc, char **argv);
	~QtApp_Posp_2DMPM();
	
	inline QtGLView& get_view() { return window.get_view(); }
	
	inline void set_win_size(int wd, int ht) { window.resize(wd, ht); }

	inline void set_display_bg_mesh(bool op = true)
	{ scene.set_display_bg_mesh(op); }
	inline void set_display_pcls(bool op = true)
	{ scene.set_display_pcls(op); }

	inline void set_res_file(ResultFile_hdf5 &rf,
		const char *th_na, size_t frame_id, const char *fld_na)
	{ scene.set_res_file(rf, th_na, frame_id, fld_na); }
	
	inline void set_fld_range(double min, double max)
	{ scene.set_fld_range(min, max); }

	int start();
};

#endif