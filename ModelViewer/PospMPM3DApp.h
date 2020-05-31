#ifndef __Posp_MPM_3D_App_h__
#define __Posp_MPM_3D_App_h__

#include <QApplication>
#include <QMainWindow>
#include <QStatusBar>

#include "MPM3DModelView.h"
#include "ResultFile_hdf5.h"

class PospMPM3DApp;

namespace PospMPM3DApp_Internal
{
	// main window of PospMPM3DApp
	class MainWindow : public QMainWindow
	{
		Q_OBJECT

	protected:
		friend PospMPM3DApp;

		MPM3DModelView* model_view;
		QStatusBar* status_bar;

	public:
		MainWindow();
		~MainWindow();
	};
};


class PospMPM3DApp
{
public:
	enum DisplayType
	{
		Invalid = -1,
		SingleFrame = 0,
		Animation = 1
	};

protected:
	QApplication app;
	DisplayType type;

	PospMPM3DApp_Internal::MainWindow *main_win;
	MPM3DModelView::Controller *view_controller;

public:
	PospMPM3DApp(int& argc, char* argv[], DisplayType dt);
	~PospMPM3DApp();

	MPM3DModelView& get_model_view() { return *(main_win->model_view); }
	
	inline void set_view_dir(float x, float y, float z) { main_win->model_view->set_view_dir(x, y, z); }
	inline void set_view_dir(float theta, float fai) { main_win->model_view->set_view_dir(theta, fai); }
	inline void set_view_dist_scale(float dist_sc) { main_win->model_view->set_view_dist_scale(dist_sc); }

	int start();
	
	// for single frame display
	int set_res_file(ResultFile_hdf5& rf, const char *th_na,
					 size_t frame_id, const char *field_na);

	// for animation generation
	void set_ani_time(double ani_time);
	int set_res_file(ResultFile_hdf5 &rf, const char *th_na,
					 const char *field_na);
	void set_gif_name(const char *gif_na);

	typedef ValueToColor::Colori Colori;
	inline int init_color_scale(
		double lower, double upper,
		Colori* colors, size_t color_num,
		bool out_of_bound_color = true
		)
	{
		return main_win->model_view->init_color_scale(lower, upper,
					colors, color_num, out_of_bound_color);
	}
};

#endif