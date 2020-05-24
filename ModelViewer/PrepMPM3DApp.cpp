#include "ModelViewer_pcp.h"

#include <QDesktopWidget>

#include "PrepMPM3DApp.h"

namespace PreMPM3DApp_Internal
{

MainWindow::MainWindow()
{
	setObjectName(QString::fromUtf8("model_view"));

	QDesktopWidget* desktop = QApplication::desktop();
	QRect screen_rect = desktop->screenGeometry();
	int screen_width = screen_rect.width();
	int screen_height = screen_rect.height();
	int win_size = (screen_width > screen_height ? screen_height : screen_width) * 2/3;
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

PrepMPM3DApp::PrepMPM3DApp(int &argc, char* argv[]) :
	app(argc, argv), main_win(new PreMPM3DApp_Internal::MainWindow) {}

PrepMPM3DApp::~PrepMPM3DApp()
{
	delete main_win;
}

void PrepMPM3DApp::set_view_dir_and_scale(float x, float y, float z, float dist_sc)
{
	main_win->model_view->set_view_dir(x, y, z);
	main_win->model_view->set_view_dist_scale(dist_sc);
}

void PrepMPM3DApp::set_view_dir_and_scale(float theta, float fai, float dist_sc)
{
	main_win->model_view->set_view_dir(theta, fai);
	main_win->model_view->set_view_dist_scale(dist_sc);
}

int PrepMPM3DApp::start()
{
	main_win->show();
	return app.exec();
}
