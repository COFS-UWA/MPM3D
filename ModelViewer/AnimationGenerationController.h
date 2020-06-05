#ifndef __Animation_Generation_Controller_h__
#define __Animation_Generation_Controller_h__

#include <QObject>
#include <QTimer>

#include "ItemArray.hpp"
#include "ResultFile_hdf5.h"
#include "MPM3DModelView.h"
#include "GifCreator.h"

// read result from hdf5 file
// generate and display animation on MPM3DModelView
// write animation to gif file
class AnimationGenerationController : public QObject,
	public MPM3DModelView::Controller
{
	Q_OBJECT
public:
	typedef MPM3DModelView::PclShape PclShape;

protected:
	QTimer ani_timer;

	// set by set_ani_time()
	double ani_time; // in ms
	
	// cal by initialize()
	double md_time;
	size_t frame_num;
	
	bool has_next_frame;
	// cal by find_next_frame()
	size_t next_frame_id;
	// time interval between current and next frame
	double ani_delay;
	double md_delay;

	double ani_div_md;
	// minimum time interval between each frame
	double min_ani_delay;
	double min_md_delay;

	size_t cur_frame_id;
	double cur_ani_time;
	double cur_md_time;
	double next_ani_time;
	double next_md_time;

	std::chrono::system_clock::time_point prev_frame_time;

	bool animation_completed;

	ResultFile_hdf5 *res_file;
	hid_t th_id;
	std::string field_name;
	size_t pcl_num;
	hid_t pcl_dt_id;
	size_t pcl_size;
	size_t pcl_x_off;
	size_t pcl_y_off;
	size_t pcl_z_off;
	size_t pcl_vol_off;
	size_t pcl_fld_off;
	hid_t pcl_fld_type;

	std::string gif_name;
	bool gif_is_init;
	unsigned short int ani_delay_100;
	GifCreator::GifWriter gif_file;
	size_t view_width, view_height;
	QPixmap screen_pixels;
	QImage screen_img;

	PclShape pcl_shape;

public:
	AnimationGenerationController(MPM3DModelView &v);
	~AnimationGenerationController();
	void close();

	inline void set_pcl_type(PclShape shape) { pcl_shape = shape; }
	inline void set_ani_time(double _ani_time) { ani_time = _ani_time * 1000.0; }
	int set_res_file(ResultFile_hdf5& rf, const char* th_na, const char* field_na);
	inline void set_gif_name(const char* gif_na) { gif_name = gif_na; }

protected:
	int initialize_model_view_data() override;
	int before_render() override;
	int after_render() override;

signals:
	void render_finished();

protected slots:
	void prepare_next_frame();

protected: // user defined behaviour
	// update frame_num, cur_md_time and md_time
	// display the the first frame scene
	virtual int initialize();

	// update cur_frame_id, ani_delay and md_delay
	// return true if there is still next frame
	// return false if already reach the last frame
	virtual bool find_next_frame();

	// update each frame
	virtual int render();
};

#endif