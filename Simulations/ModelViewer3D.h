#ifndef __Model_Viewer_3D_h__
#define __Model_Viewer_3D_h__

#include <QApplication>

#include "ModelWindow.h"
#include "ModelToViewer3D.hpp"

template <typename ModelType>
int view_model_t3d_me(ModelType &model)
{
	int argc = 1;
	char *argv[] = { "ModelViewer3D" };
	QApplication app(argc, argv);

	ModelWindow md_win;
	ModelToViewer3D<ModelType> m2v(model, md_win);
	md_win.set_model_to_viewer(m2v);
	md_win.get_gl_win().set_view_dir(-1.0f, -2.0f, -1.0f);
	md_win.show();
	
	return app.exec();
}
#endif