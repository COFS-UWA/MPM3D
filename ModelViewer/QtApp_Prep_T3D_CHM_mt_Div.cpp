#include "ModelViewer_pcp.h"

#include <QDesktopWidget>

#include "QtApp_Prep_T3D_CHM_mt_Div.h"

namespace QtApp_Prep_T3D_CHM_mt_Div_Internal
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
