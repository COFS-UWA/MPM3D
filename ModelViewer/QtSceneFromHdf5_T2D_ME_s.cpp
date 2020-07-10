#include "ModelViewer_pcp.h"

#include "Model_hdf5_utilities.h"
#include "QtSceneFromHdf5_T2D_ME_s.h"

QtSceneFromHdf5_T2D_ME_s::QtSceneFromHdf5_T2D_ME_s(
	QOpenGLFunctions_3_3_Core &_gl) :
	QtSceneFromHdf5(_gl), res_file(nullptr),
	frame_grp_id(-1), th_id(-1), pcl_dt_id(),
	display_bg_mesh(true), display_pcls(true),
	bg_mesh_obj(_gl), pcls_obj(_gl),
	color_map_texture(0),
	display_whole_model(true), padding_ratio(0.05f),
	bg_color(0.2f, 0.3f, 0.3f)
{

}

QtSceneFromHdf5_T2D_ME_s::~QtSceneFromHdf5_T2D_ME_s()
{
	close_file();
	clear();
}

void QtSceneFromHdf5_T2D_ME_s::close_file()
{
	if (!res_file)
		return;
	ResultFile_hdf5& rf = *res_file;
	if (frame_grp_id >= 0)
	{
		rf.close_group(frame_grp_id);
		frame_grp_id = -1;
	}
	if (th_id >= 0)
	{
		rf.close_group(th_id);
		th_id = -1;
	}
	if (pcl_dt_id >= 0)
	{
		H5Tclose(pcl_dt_id);
		pcl_dt_id = -1;
	}
	res_file = nullptr;
}

void QtSceneFromHdf5_T2D_ME_s::clear()
{
	if (color_map_texture)
	{
		gl.glDeleteTextures(1, &color_map_texture);
		color_map_texture = 0;
	}
}

void QtSceneFromHdf5_T2D_ME_s::set_viewport(
	int wd, int ht, GLfloat xlen, GLfloat ylen)
{
	int wd2, ht2, padding;
	ht2 = wd * (ylen / xlen);
	if (ht2 <= ht)
	{
		padding = (ht - ht2) / 2;
		vp_x_pos = 0;
		vp_y_pos = padding;
		vp_x_size = wd;
		vp_y_size = ht2;
	}
	else
	{
		wd2 = ht * (xlen / ylen);
		padding = (wd - wd2) / 2;
		vp_x_pos = padding;
		vp_y_pos = 0;
		vp_x_size = wd2;
		vp_y_size = ht;
	}
}

void QtSceneFromHdf5_T2D_ME_s::draw()
{
	gl.glViewport(vp_x_pos, vp_y_pos, vp_x_size, vp_y_size);

	gl.glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	gl.glClear(GL_COLOR_BUFFER_BIT);

	shader_plain2D.bind();

	if (display_bg_mesh)
		bg_mesh_obj.draw(shader_plain2D);

	shader_circles.bind();

	if (display_pcls)
		pcls_obj.draw(shader_circles);
}

void QtSceneFromHdf5_T2D_ME_s::resize(int wd, int ht)
{
	set_viewport(wd, ht, xu - xl, yu - yl);
}

// ================== animation ==================
int QtSceneFromHdf5_T2D_ME_s::set_res_file(
	ResultFile_hdf5& rf,
	const char* th_na,
	const char* field_na
	)
{
	close_file();
	res_file = &rf;

	hid_t th_grp_id = rf.get_time_history_grp_id();
	// check if the time history exists
	if (!rf.has_group(th_grp_id, th_na))
	{
		close_file();
		return -1;
	}
	th_id = rf.open_group(th_grp_id, th_na);

	if (!rf.has_group(th_id, "frame_0"))
	{
		close_file();
		return -2;
	}
	frame_grp_id = rf.open_group(th_id, "frame_0");

	hid_t pcl_grp_id = rf.open_group(frame_grp_id, "ParticleData");
	hid_t pcl_dset_id = rf.open_dataset(pcl_grp_id, "field");

	// get field data type
	pcl_dt_id = H5Dget_type(pcl_dset_id);
	pcl_size = H5Tget_size(pcl_dt_id);
	int init_mem_num = 0;
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char* mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
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
		else if (strcmp(mem_name, "vol") == 0)
		{
			pcl_vol_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}

		if (strcmp(mem_name, field_na) == 0)
		{
			field_name = field_na;
			pcl_fld_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			pcl_fld_type = H5Tget_member_class(pcl_dt_id, mem_id);
			++init_mem_num;
		}
	}

	rf.close_dataset(pcl_dset_id);
	rf.close_group(pcl_grp_id);
	if (init_mem_num == 4)
		return 0;

	// not enough members in data type
	close_file();
	return -3;
}

size_t QtSceneFromHdf5_T2D_ME_s::get_frame_num()
{
	ResultFile_hdf5& rf = *res_file;
	size_t frame_num;
	rf.read_attribute(th_id, "output_num", frame_num);
	return frame_num;
}

double QtSceneFromHdf5_T2D_ME_s::get_frame_time(size_t frame_id)
{
	ResultFile_hdf5& rf = *res_file;
	double frame_time;
	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	rf.read_attribute(frame_grp_id, "total_time", frame_time);
	rf.close_group(frame_grp_id);
	return frame_time;
}

int QtSceneFromHdf5_T2D_ME_s::init_scene(int wd, int ht, size_t frame_id)
{
	using namespace Model_hdf5_utilities;

	if (!res_file || !res_file->is_open())
		return -1;

	int res;
	ResultFile_hdf5& rf = *res_file;

	// write model data to view
	hid_t md_data_grp_id = rf.get_model_data_grp_id();
	hid_t bg_mesh_id = rf.open_group(md_data_grp_id, "BackgroundMesh");
	// read nodes
	size_t node_num;
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	hid_t node_dt_id = get_nd_2d_dt_id();
	Node2DData* nodes_data = new Node2DData[node_num];
	rf.read_dataset(bg_mesh_id, "NodeCoordinate", node_num, nodes_data, node_dt_id);
	H5Tclose(node_dt_id);
	// read elements
	size_t elem_num;
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	hid_t elem_dt_id = get_ed_2d_dt_id();
	Elem2DData* elems_data = new Elem2DData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementTopology", elem_num, elems_data, elem_dt_id);
	H5Tclose(elem_dt_id);
	// init bg_mesh buffer
	QVector3D gray(0.5f, 0.5f, 0.5f);
	res = bg_mesh_obj.init_from_elements(
		nodes_data,
		node_num,
		elems_data,
		elem_num,
		gray
		);
	// get bounding box
	if (display_whole_model)
	{
		xl = nodes_data[0].x;
		xu = xl;
		yl = nodes_data[0].y;
		yu = yl;
		for (size_t n_id = 1; n_id < node_num; ++n_id)
		{
			Node2DData& n = nodes_data[n_id];
			if (xl > n.x)
				xl = n.x;
			if (xu < n.x)
				xu = n.x;
			if (yl > n.y)
				yl = n.y;
			if (yu < n.y)
				yu = n.y;
		}
		GLfloat xlen = xu - xl;
		GLfloat ylen = yu - yl;
		GLfloat padding = (xlen > ylen ? xlen : ylen) * padding_ratio;
		xl -= padding;
		xu += padding;
		yl -= padding;
		yu += padding;
	}
	delete[] nodes_data;
	delete[] elems_data;
	rf.close_group(bg_mesh_id);
	if (res)
		return res;

	// init particle data buffer
	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	hid_t pcl_data_id = rf.open_group(frame_grp_id, "ParticleData");
	// pcl num
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	// pcl data
	pcls_data_mem.reserve(pcl_size * pcl_num);
	char* pcls_data = pcls_data_mem.get_mem();
	rf.read_dataset(pcl_data_id, "field", pcl_num, (void *)pcls_data, pcl_dt_id);
	rf.close_group(pcl_data_id);
	rf.close_group(frame_grp_id);
	// init pcl gl buffer
	res = pcls_obj.init<double>(
		pcls_data, pcl_size, pcl_num,
		pcl_x_off, pcl_y_off,
		pcl_vol_off, 0.5, pcl_fld_off
		);
	if (res)
		return res;

	// color map texture
	size_t color_map_texture_size;
	unsigned char* color_map_texture_data
		= color_map.gen_1Dtexture(20, color_map_texture_size);
	if (!color_map_texture_data || !color_map_texture_size)
		return -2;
	gl.glGenTextures(1, &color_map_texture);
	if (color_map_texture == 0)
		return -1;
	gl.glBindTexture(GL_TEXTURE_1D, color_map_texture);
	gl.glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	gl.glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	gl.glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	gl.glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB,
		color_map_texture_size, 0, GL_RGB,
		GL_UNSIGNED_BYTE, color_map_texture_data);
	gl.glGenerateMipmap(GL_TEXTURE_1D);
	delete[] color_map_texture_data;

	gl.glActiveTexture(GL_TEXTURE0);
	gl.glBindTexture(GL_TEXTURE_1D, color_map_texture);

	// viewport
	set_viewport(wd, ht, xu - xl, yu - yl);

	// init shaders
	// shader_plain2D
	shader_plain2D.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_plain2D.vert"
	);
	shader_plain2D.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_plain2D.frag"
	);
	shader_plain2D.link();

	// shader_circles
	shader_circles.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_circles.vert"
		);
	shader_circles.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_circles.frag"
		);
	shader_circles.link();

	// view matrix
	view_mat.setToIdentity();
	view_mat.ortho(xl, xu, yl, yu, -1.0f, 1.0f);

	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", view_mat);

	shader_circles.bind();
	shader_circles.setUniformValue("view_mat", view_mat);
	// color map
	GLfloat color_value_lower = color_map.get_lower_bound();
	shader_circles.setUniformValue("color_value_lower", color_value_lower);
	GLfloat color_value_range = color_map.get_range();
	shader_circles.setUniformValue("color_value_range", color_value_range);
	shader_circles.setUniformValue("color_map", 0);

	return 0;
}

void QtSceneFromHdf5_T2D_ME_s::update_scene(size_t frame_id)
{
	ResultFile_hdf5 &rf = *res_file;
	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	hid_t pcl_data_id = rf.open_group(frame_grp_id, "ParticleData");
	// pcl num
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	// pcl data
	pcls_data_mem.reserve(pcl_size * pcl_num);
	char* pcls_data = pcls_data_mem.get_mem();
	rf.read_dataset(pcl_data_id, "field", pcl_num, (void*)pcls_data, pcl_dt_id);
	rf.close_group(pcl_data_id);
	rf.close_group(frame_grp_id);
	// update pcl gl buffer
	pcls_obj.update<double>(
		pcls_data, pcl_size, pcl_num,
		pcl_x_off, pcl_y_off,
		pcl_vol_off, 0.5, pcl_fld_off
		);
}
