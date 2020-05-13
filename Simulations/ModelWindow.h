#ifndef __Model_Window_h__
#define __Model_Window_h__

#include <QApplication>

#include <QMainWindow>
#include <QStatusBar>

#include "GLModelWindow.h"

class ModelWindow : public QMainWindow
{
	Q_OBJECT
protected:
	GLModelWindow *modelWindow;
	QStatusBar *statusbar;

public:
	ModelWindow();
	~ModelWindow();

	void set_model_to_viewer(ModelToViewerBase &m2v) { modelWindow->set_controller(m2v); }
	GLModelWindow &get_gl_win() { return *modelWindow; }
};

#endif