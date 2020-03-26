#include "QtPostProcessor_pcp.h"

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

	int status = a.exec();
	//system("pause");
	return status;
}
