#include "Simulations_pcp.h"

#include "ModelWindow.h"

ModelWindow::ModelWindow()
{
	setObjectName(QString::fromUtf8("mainWindow"));
	resize(800, 600);

	modelWindow = new GLModelWindow(this);
	modelWindow->setObjectName(QString::fromUtf8("modelWindow"));
	setCentralWidget(modelWindow);

	statusbar = new QStatusBar(this);
	statusbar->setObjectName(QString::fromUtf8("statusbar"));
	setStatusBar(statusbar);
}

ModelWindow::~ModelWindow() {}



