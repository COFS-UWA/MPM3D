#include "ModelViewer_pcp.h"

#include <limits>
#include <iostream>

#include "QtController_Posp_Animation.h"

QtController_Posp_Animation::QtController_Posp_Animation() :
	view(nullptr), scene(nullptr)
{
	init_self();
}

QtController_Posp_Animation::QtController_Posp_Animation(
	QtGLView& v,
	QtSceneFromHdf5_2DMPM& s
	)
{
	set_view(v);
	set_scene(s);
	init_self();
}

void QtController_Posp_Animation::init_self()
{
	min_ani_delay = 20.0 * 0.9999; // 50 fps

	cur_frame_id = 0;
	end_frame_id = std::numeric_limits<size_t>::max();

	ani_timer.setSingleShot(true);
	connect(this, SIGNAL(render_next_frame_signal()),
			this, SLOT(render_next_frame()));
	connect(&ani_timer, SIGNAL(timeout()),
			this, SLOT(draw_cur_frame()));
}

QtController_Posp_Animation::~QtController_Posp_Animation()
{

}

int QtController_Posp_Animation::initialize(int wd, int ht)
{	
	size_t max_output_num = scene->get_frame_num();
	if (end_frame_id > max_output_num)
		end_frame_id = max_output_num;
	if (cur_frame_id > max_output_num - 1)
		cur_frame_id = max_output_num - 1;
	if (cur_frame_id > end_frame_id - 1)
		cur_frame_id = end_frame_id - 1;

	cur_ani_time = 0.0;
	cur_md_time = scene->get_frame_time(cur_frame_id);
	double end_md_time = scene->get_frame_time(end_frame_id-1);
	md_time = end_md_time - cur_md_time;
	ani_div_md = ani_time / md_time;
	min_md_delay = min_ani_delay / ani_div_md;
	has_next_frame = true;

	// draw the first frame
	return scene->init_scene(wd, ht, cur_frame_id);
}

namespace
{
inline void print_frame_info(size_t cur_frame_id,
	double cur_ani_time, double cur_md_time)
{
	std::cout << "Display frame " << cur_frame_id
		<< " at " << cur_ani_time / 1000.0
		<< "s (model time: " << cur_md_time << ")\n";
}
}

void QtController_Posp_Animation::draw_scene()
{
	print_frame_info(cur_frame_id, cur_ani_time, cur_md_time);
	scene->draw();
	if (has_next_frame)
	{
		prev_frame_time = std::chrono::system_clock::now();
		emit render_next_frame_signal();
	}
}

void QtController_Posp_Animation::resize_scene(int wd, int ht)
{
	(void)wd; (void)ht; // fixed size win, do nothing here
}

void QtController_Posp_Animation::render_next_frame()
{
	if((has_next_frame = find_next_frame()) == false)
		return;

	scene->update_scene(cur_frame_id);
	
	std::chrono::system_clock::time_point cur_time = std::chrono::system_clock::now();
	std::chrono::milliseconds elapsed_time
		= std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - prev_frame_time);
	double time_to_next_frame = ani_delay - double(elapsed_time.count());
	if (time_to_next_frame < 0.0)
		time_to_next_frame = 0.0;
	ani_timer.start(time_to_next_frame);
}

void QtController_Posp_Animation::draw_cur_frame() { view->update(); }

bool QtController_Posp_Animation::find_next_frame()
{
	double new_md_time, diff_md_time, diff_ani_time;
	for (size_t frame_id = cur_frame_id + 1; frame_id < end_frame_id; ++frame_id)
	{
		new_md_time = scene->get_frame_time(frame_id);

		diff_md_time = new_md_time - cur_md_time;
		diff_ani_time = diff_md_time * ani_div_md;
		if (diff_ani_time > min_ani_delay)
		{
			ani_delay = diff_ani_time;
			md_delay = diff_md_time;
			cur_frame_id = frame_id;
			cur_ani_time += ani_delay;
			cur_md_time = new_md_time;
			return true;
		}
	}

	return false;
}
