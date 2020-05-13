#include "ResultViewerQt_pcp.h"

#include <iostream>
#include <cmath>

#include <QtWidgets/QApplication>

#include "MainWindow.h"

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);

	MainWindow w;
	w.show();

	emit w.action_open->triggered();

	return a.exec();
}
