#include "ModelViewer_pcp.h"

#include <limits>
#include <iostream>

#include "file_utils.h"

#include "QtController_Posp_Animation.h"

QtController_Posp_Animation::QtController_Posp_Animation() :
	scene(nullptr), has_frame_undrawed(false)
{
	init_self();
}

QtController_Posp_Animation::QtController_Posp_Animation(
	QtGLView& v,
	QtSceneFromHdf5& s) :
	scene(nullptr), has_frame_undrawed(false)
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

	ani_next_frame_timer.setSingleShot(true);
	connect(&ani_next_frame_timer, SIGNAL(timeout()),
			this, SLOT(render_next_frame()));

	ani_timer.setSingleShot(true);
	connect(&ani_timer, SIGNAL(timeout()),
			view, SLOT(update()));

	// output parameters
	need_output_png = false;
	png_name = "";

	need_output_gif = false;
	gif_name = "";
	gif_file_is_init = false;
}

QtController_Posp_Animation::~QtController_Posp_Animation()
{
	close_gif_file();
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

	render_begin_time = std::chrono::system_clock::now();

	// init png output
	if (need_output_png)
	{
		make_dir(png_name.c_str());
	}

	// draw the first frame
	int res = scene->init_scene(wd, ht, cur_frame_id);
	has_frame_undrawed = true;

	return res;
}

namespace
{
inline void print_frame_info(size_t cur_frame_id,
	double cur_ani_time, double cur_md_time, double render_time)
{
	std::cout << "Frame " << cur_frame_id
		<< " at animation time: " << cur_ani_time / 1000.0
		<< "s (model time: " << cur_md_time << "s), using "
		<< render_time << "s\n";
}
}

void QtController_Posp_Animation::draw_scene()
{
	scene->draw();
	if (has_frame_undrawed)
	{
		has_frame_undrawed = false;
		prev_frame_time = std::chrono::system_clock::now();
		render_time = double(std::chrono::duration_cast<std::chrono::milliseconds>(prev_frame_time - render_begin_time).count())/1000.0;
		print_frame_info(cur_frame_id, cur_ani_time, cur_md_time, render_time);
		if (has_next_frame)
			ani_next_frame_timer.start(0);
	}
}

namespace
{
	void from_bgra_to_rgba(uchar* dst, uchar* src, size_t width, size_t height)
	{
		for (size_t h_id = 0; h_id < height; ++h_id)
			for (size_t w_id = 0; w_id < width; ++w_id)
			{
				dst[0] = src[2];
				dst[1] = src[1];
				dst[2] = src[0];
				dst[3] = src[3];
				dst += 4;
				src += 4;
			}
	}
}

void QtController_Posp_Animation::render_next_frame()
{
	// output prev frame
	char ss_name[256];
	if (need_output_png)
	{
		snprintf(ss_name, 256, "%s\\frame_%zu.png",
				 png_name.c_str(), cur_frame_id);
		view->grabFramebuffer().save(ss_name, "png");
	}

	if (need_output_gif)
	{
		if (!gif_file_is_init)
		{
			gif_width = view->width();
			gif_height = view->height();
			GifCreator::GifBegin(&gif_file,
				gif_name.c_str(),
				gif_width, gif_height, true
			);
			gif_file_is_init = true;
			screen_img_rgba.reserve(gif_width * gif_height * 4);
		}
		
		ani_delay_100 = unsigned short(ani_delay / 10.0);
		if (ani_delay_100 < 2)
			ani_delay_100 = 2;
		screen_img = view->grabFramebuffer();
		from_bgra_to_rgba(
			screen_img_rgba.get_mem(),
			screen_img.bits(),
			gif_width,
			gif_height
			);
		GifCreator::GifWriteFrame(
			&gif_file,
			screen_img_rgba.get_mem(),
			gif_width,
			gif_height,
			ani_delay_100
		);
	}

	// has_next_frame
	if ((has_next_frame = find_next_frame()) == false)
	{
		if (need_output_gif)
		{
			// output the last frame to gif
			GifCreator::GifWriteFrame(
				&gif_file,
				screen_img_rgba.get_mem(),
				gif_width,
				gif_height,
				200 // 2s
				);
			close_gif_file();
		}
		return;
	}

	scene->update_scene(cur_frame_id);
	has_frame_undrawed = true;

	std::chrono::system_clock::time_point cur_time = std::chrono::system_clock::now();
	std::chrono::milliseconds elapsed_time
		= std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - prev_frame_time);
	// ANI_DELAY_ADJUSTMENT is an empirical number to make ani schedule more accurate
#define ANI_DELAY_ADJUSTMENT 10.0
	double time_to_next_frame = ani_delay - ANI_DELAY_ADJUSTMENT - double(elapsed_time.count());
#undef  ANI_DELAY_ADJUSTMENT
	if (time_to_next_frame < 0.0)
		time_to_next_frame = 0.0;
	ani_timer.start(time_to_next_frame);
}

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

void QtController_Posp_Animation::set_png_name(const char* fname)
{
	if (!fname || !strcmp(fname, ""))
	{
		need_output_png = false;
		png_name = "";
		return;
	}
	need_output_png = true;
	png_name = fname;
}

void QtController_Posp_Animation::set_gif_name(const char* fname)
{
	if (!fname || !strcmp(fname, ""))
	{
		need_output_gif = false;
		gif_name = "";
		return;
	}
	need_output_gif = true;
	gif_name = std::string(fname) + ".gif";
}

void QtController_Posp_Animation::close_gif_file()
{
	if (gif_file_is_init)
	{
		GifCreator::GifEnd(&gif_file);
		gif_file_is_init = false;
	}
}

void QtController_Posp_Animation::resize_scene(int wd, int ht)
{
	(void)wd; (void)ht; // fixed size win, do nothing
}
