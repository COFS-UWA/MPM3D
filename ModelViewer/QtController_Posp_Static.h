#ifndef __Qt_Controller_Posp_Static_h__
#define __Qt_Controller_Posp_Static_h__

#include <QTimer>

#include "QtSceneFromHdf5.h"
#include "QtGLView.h"
#include "QtController.h"

class QtController_Posp_Static : public QtController
{
	Q_OBJECT

protected:
	QtSceneFromHdf5 *scene;

	size_t frame_id; // frame to be displayed
	
	void init_self();

public:
	QtController_Posp_Static();
	QtController_Posp_Static(QtGLView &v, QtSceneFromHdf5 &s);
	~QtController_Posp_Static();

	inline void set_scene(QtSceneFromHdf5& s) { scene = &s; }
	inline int set_frame_id(size_t f_id)
	{
		size_t frame_num = scene->get_frame_num();
		if (f_id >= frame_num)
			return -1;
		frame_id = f_id;
		return 0;
	}

	int initialize(int wd, int ht) override;
	void draw_scene() override;
	void resize_scene(int wd, int ht) override;

protected: // png output
	bool need_output_png;
	std::string png_name;
	QTimer screen_shot_timer;

protected slots:
	void window_to_png();

public:
	void set_png_name(const char* fname);
};

#endif