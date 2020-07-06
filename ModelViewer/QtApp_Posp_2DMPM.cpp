#include "ModelViewer_pcp.h"

#include <QDesktopWidget>

#include "QtApp_Posp_2DMPM.h"

namespace QtApp_Posp_2DMPM_Internal
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

QtApp_Posp_2DMPM::QtApp_Posp_2DMPM(int &argc, char **argv) :
	app(argc, argv),
	scene(window.get_view()),
	controller(window.get_view(), scene) {}

QtApp_Posp_2DMPM::~QtApp_Posp_2DMPM() {}

int QtApp_Posp_2DMPM::start()
{
	window.show();
	return app.exec();
}
