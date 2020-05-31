#include "ModelViewer_pcp.h"

#include <iostream>
#include <chrono>

#include "Model_hdf5_utilities.h"
#include "Model_T3D_ME_s_hdf5_utilities.h"
#include "AnimationGenerationController.h"

AnimationGenerationController::AnimationGenerationController(MPM3DModelView& v) :
	MPM3DModelView::Controller(v),
	frame_num(0), md_time(0.0),
	min_ani_delay(20.0 * 0.9999), // maximum frame rate is 50 fps
	res_file(nullptr), th_id(-1), field_name(""), pcl_dt_id(-1),
	gif_name(""), gif_is_init(false),
	animation_completed(false)
{
	connect(this, SIGNAL(render_finished()), this, SLOT(prepare_next_frame()));
	
	ani_timer.setSingleShot(true);
	connect(&ani_timer, SIGNAL(timeout()), view, SLOT(update()));
}

AnimationGenerationController::~AnimationGenerationController() { close(); }

void AnimationGenerationController::close()
{
	// close hdf5 file
	if (th_id >= 0)
	{
		res_file->close_group(th_id);
		th_id = -1;
	}
	if (pcl_dt_id >= 0)
	{
		H5Tclose(pcl_dt_id);
		pcl_dt_id = -1;
	}

	// close gif
	if (gif_is_init)
	{
		GifCreator::GifEnd(&gif_file);
		gif_is_init = false;
	}
}

int AnimationGenerationController::set_res_file(
	ResultFile_hdf5& rf,
	const char *th_na,
	const char *field_na
	)
{
	close();
	res_file = &rf;
	hid_t th_grp_id = rf.get_time_history_grp_id();
	
	// check if the time history exists
	if (!rf.has_group(th_grp_id, th_na))
	{
		close();
		return -1;
	}

	th_id = rf.open_group(th_grp_id, th_na);
	
	if (!rf.has_group(th_id, "frame_0"))
		return -2;

	hid_t frame_grp_id = rf.open_group(th_id, "frame_0");
	hid_t pcl_grp_id = rf.open_group(frame_grp_id, "ParticleData");
	hid_t pcl_dset_id = rf.open_dataset(pcl_grp_id, "field");
	pcl_dt_id = H5Dget_type(pcl_dset_id);
	pcl_size = H5Tget_size(pcl_dt_id);
	int init_mem_num = 0;
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char *mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
		if (strcmp(mem_name, "x") == 0)
		{
			pcl_x_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, "y") == 0)
		{
			pcl_y_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, "z") == 0)
		{
			pcl_z_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, "vol") == 0)
		{
			pcl_vol_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, field_na) == 0)
		{
			field_name = field_na;
			pcl_fld_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			pcl_fld_type = H5Tget_member_class(pcl_dt_id, mem_id);
			++init_mem_num;
		}
	}

	rf.close_dataset(pcl_dset_id);
	rf.close_group(pcl_grp_id);
	rf.close_group(frame_grp_id);
	if (init_mem_num == 5)
		return 0;

	// no such data member
	return -3;
}

int AnimationGenerationController::initialize_model_view_data()
{
	cur_frame_id = 0;
	cur_ani_time = 0.0;
	cur_md_time = 0.0;

	int res = initialize();
	if (res < 0)
		return res;

	if (frame_num == 0 || md_time <= 0.0)
		return -1;

	ani_div_md = ani_time / md_time;
	min_md_delay = min_ani_delay / ani_div_md;

	has_next_frame = find_next_frame();

	return 0;
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

int AnimationGenerationController::before_render()
{
	if (!animation_completed)
		print_frame_info(cur_frame_id, cur_ani_time, cur_md_time);
	return 0;
}

int AnimationGenerationController::after_render()
{
	// prepare next frame
	if (has_next_frame)
		prev_frame_time = std::chrono::system_clock::now();
	if (!animation_completed);
		emit render_finished();
	return 0;
}

void AnimationGenerationController::prepare_next_frame()
{	
	if (animation_completed)
		return;

	// init gif output
	if (!gif_is_init && gif_name.length())
	{
		view_width = view->width();
		view_height = view->height();
		GifCreator::GifBegin(&gif_file, gif_name.c_str(), view_width, view_height, true);
		gif_is_init = true;
	}

	// output frame to gif
	if (gif_is_init)
	{
		ani_delay_100 = unsigned short int(ani_delay / 10.0);
		if (ani_delay_100 < 2)
			ani_delay_100 = 2;
		screen_pixels = view->grab(view->geometry());
		//char ss_name[50];
		//snprintf(ss_name, 50, "frame_%zu.png", cur_frame_id);
		//screen_pixels.save(ss_name, "png");
		//std::cout << "output frame " << cur_frame_id << "\n";
		screen_img = screen_pixels.toImage();
		GifCreator::GifWriteFrame(
			&gif_file,
			screen_img.bits(), //(const uint8_t *)pixels_rgb_data,
			view_width,
			view_height,
			ani_delay_100
			);
		if (!has_next_frame) // the last frame
		{
			GifCreator::GifWriteFrame(
				&gif_file,
				screen_img.bits(),
				view_width,
				view_height,
				2
				);
		}
	}
	
	// finalize gif output
	if (!has_next_frame)
	{
		if (gif_is_init)
		{
			GifCreator::GifEnd(&gif_file);
			gif_is_init = false;
			animation_completed = true;
		}
		return;
	}

	cur_frame_id = next_frame_id;
	cur_ani_time += ani_delay;
	cur_md_time += md_delay;

	render();

	has_next_frame = find_next_frame();

	std::chrono::system_clock::time_point cur_time = std::chrono::system_clock::now();
	std::chrono::milliseconds elapsed_time
		= std::chrono::duration_cast<std::chrono::milliseconds>(cur_time - prev_frame_time);
	double time_to_next_frame = ani_delay - double(elapsed_time.count());
	if (time_to_next_frame < 0.0)
		time_to_next_frame = 0.0;
	ani_timer.start(time_to_next_frame); // display this frame
}


// user defined part
using namespace Model_hdf5_utilities;
using namespace Model_T3D_ME_s_hdf5_utilities;

int AnimationGenerationController::initialize()
{
	ResultFile_hdf5& rf = *res_file;
	
	// get frame_num
	rf.read_attribute(th_id, "output_num", frame_num);
	if (frame_num == 0)
		return -1;

	// get cur_md_time
	hid_t first_frame_id = rf.open_group(th_id, "frame_0");
	rf.read_attribute(first_frame_id, "total_time", cur_md_time);
	rf.close_group(first_frame_id);

	// get md_time
	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_num-1);
	hid_t last_frame_id = rf.open_group(th_id, frame_name);
	rf.read_attribute(last_frame_id, "total_time", md_time);
	rf.close_group(last_frame_id);
	md_time -= cur_md_time;

	// write model data to view
	hid_t md_data_grp_id = rf.get_model_data_grp_id();
	hid_t bg_mesh_id = rf.open_group(md_data_grp_id, "BackgroundMesh");
	size_t node_num, elem_num;
	// read nodes
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	hid_t node_dt_id = get_nd_3d_dt_id();
	Node3DData *nodes_data = new Node3DData[node_num];
	rf.read_dataset(bg_mesh_id, "NodeCoordinate", node_num, nodes_data, node_dt_id);
	H5Tclose(node_dt_id);
	// read elements
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	hid_t elem_dt_id = get_ed_3d_dt_id();
	Elem3DData *elems_data = new Elem3DData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementTopology", elem_num, elems_data, elem_dt_id);
	H5Tclose(elem_dt_id);
	// write nodess and elems to window
	view->init_bg_mesh<Node3DData, Elem3DData>(
		nodes_data, node_num,
		elems_data, elem_num,
		QVector3D(0.5f, 0.5f, 0.5f)
		);
	delete[] nodes_data;
	delete[] elems_data;
	rf.close_group(bg_mesh_id);

	// write particle data to view
	hid_t frame_id = rf.open_group(th_id, "frame_0");
	hid_t pcl_data_id = rf.open_group(frame_id, "ParticleData");
	// pcl num
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	// pcl data
	char* pcls_data = new char[pcl_size * pcl_num];
	rf.read_dataset(pcl_data_id, "field", pcl_num, (void *)pcls_data, pcl_dt_id);
	rf.close_group(pcl_data_id);
	rf.close_group(frame_id);
	// write pcl data to view
	// support field data of different type
	//if (data_type == H5T_NATIVE_DOUBLE)
	view->init_multicolor_pcl_data<double>(
		pcls_data, pcl_size, pcl_num,
		pcl_x_off, pcl_y_off, pcl_z_off,
		pcl_vol_off, 0.125,
		pcl_fld_off, view->get_color_scale());
	delete[] pcls_data;

	return 0;
}

bool AnimationGenerationController::find_next_frame()
{
	size_t frame_id;
	double new_frame_md_time, diff_md_time, diff_ani_time;
	char frame_name[50];
	hid_t frame_h5_id;
	ResultFile_hdf5 &rf = *res_file;
	
	for (frame_id = cur_frame_id + 1; frame_id < frame_num; ++frame_id)
	{
		snprintf(frame_name, 50, "frame_%zu", frame_id);
		frame_h5_id = rf.open_group(th_id, frame_name);
		rf.read_attribute(frame_h5_id, "total_time", new_frame_md_time);
		rf.close_group(frame_h5_id);

		diff_md_time = new_frame_md_time - cur_md_time;
		diff_ani_time = diff_md_time * ani_div_md;
		if (diff_ani_time > min_ani_delay)
		{
			ani_delay = diff_ani_time;
			md_delay = diff_md_time;
			next_frame_id = frame_id;
			return true;
		}
	}

	// not more frame
	ani_delay = 2500.0; // 1.0s
	md_delay = ani_delay / ani_delay_100;
	next_frame_id = frame_num; // invalid value
	return false;
}

int AnimationGenerationController::render()
{
	ResultFile_hdf5 &rf = *res_file;
	
	// write particle data to view
	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", cur_frame_id);
	hid_t frame_id = rf.open_group(th_id, frame_name);
	hid_t pcl_data_id = rf.open_group(frame_id, "ParticleData");
	// pcl num
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	// pcl data
	char* pcls_data = new char[pcl_size * pcl_num];
	rf.read_dataset(pcl_data_id, "field", pcl_num, pcls_data, pcl_dt_id);
	rf.close_group(pcl_data_id);
	rf.close_group(frame_id);
	// write pcl data to view, support field data of different type
	//if (data_type == H5T_NATIVE_DOUBLE)
	view->update_multicolor_pcl_data<double>(pcls_data);
	delete[] pcls_data;

	return 0;
}
