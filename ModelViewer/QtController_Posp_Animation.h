#ifndef __Qt_Controller_Posp_Animation_h__
#define __Qt_Controller_Posp_Animation_h__

#include <chrono>

#include <QTimer>

#include "ItemArray.hpp"
#include "GifCreator.h"

#include "QtSceneFromHdf5_2DMPM.h"
#include "QtGLView.h"
#include "QtController.h"

class QtController_Posp_Animation : public QtController
{
	Q_OBJECT

protected:
	QtGLView *view;
	QtSceneFromHdf5_2DMPM *scene;

	// set by set_ani_time()
	double ani_time; // in ms
	double md_time;

	size_t cur_frame_id;
	size_t end_frame_id; // = last frame_id + 1

	bool has_next_frame;
	// time interval between current and next frame
	double ani_delay;
	double md_delay;

	double ani_div_md;
	// minimum time interval between each frame
	double min_ani_delay;
	double min_md_delay;

	double cur_ani_time;
	double cur_md_time;
	double next_ani_time;
	double next_md_time;

	std::chrono::system_clock::time_point prev_frame_time;
	QTimer ani_timer;

	// gif output parameters
	std::string gif_name;
	bool gif_is_init;
	unsigned short int ani_delay_100;
	GifCreator::GifWriter gif_file;
	size_t view_width, view_height;
	QPixmap screen_pixels;
	QImage screen_img;
	MemoryUtils::ItemArray<uchar> screen_img_rgba;

	void init_self();
	bool find_next_frame();

public:
	QtController_Posp_Animation();
	QtController_Posp_Animation(QtGLView &v, QtSceneFromHdf5_2DMPM &s);
	~QtController_Posp_Animation();

	inline void set_view(QtGLView& v)
	{
		view = &v;
		view->set_controller(*this);
		// fixed window size
		QSizePolicy fixed_size_policy(QSizePolicy::Fixed, QSizePolicy::Fixed);
		view->setSizePolicy(fixed_size_policy);
	}
	inline void set_scene(QtSceneFromHdf5_2DMPM& s) { scene = &s; }

	inline void set_ani_time(double ani_t /* in second */) { ani_time = ani_t * 1000.0; }
	inline void set_start_frame(size_t f_id) { cur_frame_id = f_id; }
	inline void set_end_frame(size_t f_id) { end_frame_id = f_id + 1; }

	int initialize(int wd, int ht) override;
	void draw_scene() override;
	void resize_scene(int wd, int ht) override;

signals:
	void render_next_frame_signal();

protected slots:
	void render_next_frame();
	void draw_cur_frame();
};

#endif