#include "ModelViewer_pcp.h"

#include <QDesktopWidget>

#include "PospSingleFrameController.h"
#include "AnimationGenerationController.h"

#include "PospMPM3DApp.h"

namespace PospMPM3DApp_Internal
{
	MainWindow::MainWindow()
	{
		setObjectName(QString::fromUtf8("model_view"));

		QDesktopWidget* desktop = QApplication::desktop();
		QRect screen_rect = desktop->screenGeometry();
		int screen_width = screen_rect.width();
		int screen_height = screen_rect.height();
		int win_size = (screen_width > screen_height ? screen_height : screen_width) * 2 / 3;
		resize(win_size, win_size);

		model_view = new MPM3DModelView(this);
		model_view->setObjectName(QString::fromUtf8("model_view"));
		setCentralWidget(model_view);

		status_bar = new QStatusBar(this);
		status_bar->setObjectName(QString::fromUtf8("status_bar"));
		setStatusBar(status_bar);
	}

	MainWindow::~MainWindow() {}

};

PospMPM3DApp::PospMPM3DApp(int& argc, char* argv[], DisplayType dt) :
	app(argc, argv), main_win(new PospMPM3DApp_Internal::MainWindow),
	view_controller(nullptr)
{
	type = dt;
	switch (dt)
	{
	case SingleFrame:
		view_controller = new PospSingleFrameController(get_model_view());
		break;
	case Animation:
		view_controller = new AnimationGenerationController(get_model_view());
		break;
	default:
		type = Invalid;
	}
}

PospMPM3DApp::~PospMPM3DApp()
{
	delete main_win;

	if (view_controller)
	{
		delete view_controller;
		view_controller = nullptr;
	}
}

int PospMPM3DApp::start()
{
	main_win->show();
	return app.exec();
}

// for single frame display
int PospMPM3DApp::set_res_file(
	ResultFile_hdf5& rf,
	const char* th_na,
	size_t frame_id,
	const char* field_na,
	MPM3DModelView::PclShape shape
	)
{
	if (type != SingleFrame)
		return -1;

	PospSingleFrameController& con =
		*static_cast<PospSingleFrameController *>(view_controller);
	con.set_pcl_type(shape);
	return con.set_res_file(rf, th_na, frame_id, field_na);
}

// for animation generation
void PospMPM3DApp::set_ani_time(double ani_time)
{
	if (type != Animation)
		return ;

	AnimationGenerationController &con = 
		*static_cast<AnimationGenerationController *>(view_controller);
	con.set_ani_time(ani_time);
}

int PospMPM3DApp::set_res_file(
	ResultFile_hdf5& rf, 
	const char* th_na,
	const char* field_na,
	MPM3DModelView::PclShape shape
	)
{
	if (type != Animation)
		return -1;

	AnimationGenerationController& con =
		*static_cast<AnimationGenerationController*>(view_controller);
	con.set_pcl_type(shape);
	return con.set_res_file(rf, th_na, field_na);
}

void PospMPM3DApp::set_gif_name(const char* gif_na)
{
	if (type != Animation)
		return;

	AnimationGenerationController& con =
		*static_cast<AnimationGenerationController*>(view_controller);
	return con.set_gif_name(gif_na);
}
