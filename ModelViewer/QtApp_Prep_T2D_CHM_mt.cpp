#include "ModelViewer_pcp.h"

#include <QDesktopWidget>

#include "QtApp_Prep_T2D_CHM_mt.h"

namespace QtApp_Prep_T2D_CHM_mt_Internal
{

MainWindow::MainWindow()
{
	setObjectName(QString::fromUtf8("model_view"));

	model_view = new QtGLView(this);
	model_view->setObjectName(QString::fromUtf8("model_view"));
	setCentralWidget(model_view);

	status_bar = new QStatusBar(this);
	status_bar->setObjectName(QString::fromUtf8("status_bar"));
	setStatusBar(status_bar);
}

MainWindow::~MainWindow() {}

};


QtApp_Prep_T2D_CHM_mt::QtApp_Prep_T2D_CHM_mt(int &argc, char **argv) :
	app(argc, argv),
	scene(window.get_view()),
	controller(window.get_view(), scene) {}

QtApp_Prep_T2D_CHM_mt::~QtApp_Prep_T2D_CHM_mt() {}

int QtApp_Prep_T2D_CHM_mt::start()
{
	window.show();
	return app.exec();
}
