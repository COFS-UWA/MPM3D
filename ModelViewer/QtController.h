#ifndef __Qt_Controller_h__
#define __Qt_Controller_h__

#include <QObject>

#include "QtGLView.h"

class QtController : public QObject
{
	Q_OBJECT
protected:
	QtGLView *view;

public:
	QtController() {}
	virtual ~QtController() {}

	inline void set_view(QtGLView& v)
	{
		view = &v;
		view->set_controller(*this);
	}
	
	virtual int initialize(int wd, int ht) = 0;
	virtual void draw_scene() = 0;
	virtual void resize_scene(int wd, int ht) = 0;
};

#endif