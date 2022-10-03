#include "ModelViewer_pcp.h"

#include "RigidObject/RigidObject_hdf5_utilities.h"
#include "Model_T2D_ME_mt_hdf5_utilities.h"
#include "QtSceneFromHdf5_T2D_ME_mt.h"

QtSceneFromHdf5_T2D_ME_mt::QtSceneFromHdf5_T2D_ME_mt(
	QOpenGLFunctions_3_3_Core &_gl) :
	QtSceneFromHdf5(_gl),
	res_file(nullptr), th_id(-1),
	x_fld(data_loader), y_fld(data_loader),
	vol_fld(data_loader), pfld(nullptr),
	display_bg_mesh(true), bg_mesh_obj(_gl),
	display_pcls(true), pcls_obj(_gl),
	display_rc(true), has_rc_obj(false), rc_obj(_gl),
	display_rr(true), has_rr_obj(false), rr_obj(_gl),
	display_rb(true), has_rb_obj(false), rb_obj(_gl),
	has_color_map(false), color_map_obj(_gl), color_map_texture(0),
	display_whole_model(true), padding_ratio(0.05f),
	bg_color(0.2f, 0.3f, 0.3f),
	rb_color(0.5176f, 0.4392f, 1.0f), // light_slate_blue
	mesh_color(0.5f, 0.5f, 0.5f), // gray
	pcl_is_mono_color(false) {}

QtSceneFromHdf5_T2D_ME_mt::~QtSceneFromHdf5_T2D_ME_mt()
{
	clear();
	close_file();
}

void QtSceneFromHdf5_T2D_ME_mt::clear()
{
	if (color_map_texture)
	{
		gl.glDeleteTextures(1, &color_map_texture);
		color_map_texture = 0;
	}
}

void QtSceneFromHdf5_T2D_ME_mt::close_file()
{
	if (!res_file)
		return;

	if (pfld)
	{
		delete pfld;
		pfld = nullptr;
	}
	data_loader.close_res_file();
	res_file = nullptr;
}

void QtSceneFromHdf5_T2D_ME_mt::set_viewport(
	int wd, int ht, GLfloat xlen, GLfloat ylen)
{
	win_wd = wd;
	win_ht = ht;

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

void QtSceneFromHdf5_T2D_ME_mt::draw()
{
	gl.glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	gl.glClear(GL_COLOR_BUFFER_BIT);

	// =================== mpm model ===================
	gl.glViewport(vp_x_pos, vp_y_pos, vp_x_size, vp_y_size);

	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", view_mat);

	if (display_bg_mesh)
		bg_mesh_obj.draw(shader_plain2D);

	if (has_rc_obj && display_rc)
		rc_obj.draw(shader_plain2D);

	if (has_rr_obj && display_rr)
		rr_obj.draw(shader_plain2D);

	if (has_rb_obj && display_rb)
	{
		shader_rigid_mesh_2D.bind();
		shader_rigid_mesh_2D.setUniformValue("view_mat", view_mat);
		rb_obj.draw(shader_rigid_mesh_2D);
	}

	shader_circles.bind();

	if (display_pcls)
		pcls_obj.draw(shader_circles);

	// =================== color map ===================
	gl.glViewport(0, 0, win_wd, win_ht);

	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", hud_view_mat);

	if (has_color_map)
		color_map_obj.draw(shader_plain2D, shader_char);
}

void QtSceneFromHdf5_T2D_ME_mt::resize(int wd, int ht)
{
	set_viewport(wd, ht, bbox.xu - bbox.xl, bbox.yu - bbox.yl);

	hud_view_mat.setToIdentity();
	hud_view_mat.ortho(0.0f, GLfloat(wd)/GLfloat(ht), 0.0f, 1.0f, -1.0f, 1.0f);
	shader_char.bind();
	shader_char.setUniformValue("view_mat", hud_view_mat);
}

// ================== animation ==================
int QtSceneFromHdf5_T2D_ME_mt::set_res_file(
	ResultFile_hdf5& rf,
	const char *th_name,
	Hdf5Field::FieldType fld_type)
{
	char exception_msg[256];
	close_file();
	res_file = &rf;
	data_loader.set_time_history(rf, th_name);
	th_id = data_loader.get_time_history_id();

	if (!x_fld.validate_data_type())
		throw std::exception("x field not found in hdf5 field data.");
	if (!y_fld.validate_data_type())
		throw std::exception("y field not found in hdf5 field data.");
	if (!vol_fld.validate_data_type())
		throw std::exception("vol field not found in hdf5 field data.");

	pfld = Hdf5Field::make(fld_type);
	if (!pfld)
	{
		snprintf(exception_msg, sizeof(exception_msg),
			"Field type %u is not supported.", unsigned int(fld_type));
		throw std::exception(exception_msg);
	}
	pfld->set_data_loader(data_loader);
	if (!(pfld->validate_data_type()))
	{
		snprintf(exception_msg, sizeof(exception_msg),
			"%s field not found in hdf5 field data.",
			Hdf5Field::get_name(fld_type));
		throw std::exception(exception_msg);
	}
	field_name = Hdf5Field::get_name(fld_type);

	return 0;
}

size_t QtSceneFromHdf5_T2D_ME_mt::get_frame_num()
{
	ResultFile_hdf5& rf = *res_file;
	size_t frame_num;
	rf.read_attribute(th_id, "output_num", frame_num);
	return frame_num;
}

double QtSceneFromHdf5_T2D_ME_mt::get_frame_time(size_t frame_id)
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

int QtSceneFromHdf5_T2D_ME_mt::init_scene(int wd, int ht, size_t frame_id)
{
	using namespace Model_hdf5_utilities;

	if (!res_file || !res_file->is_open())
		return -1;

	int res;
	ResultFile_hdf5 &rf = *res_file;

	// write model data to view
	hid_t md_data_grp_id = rf.get_model_data_grp_id();
	hid_t bg_mesh_id = rf.open_group(md_data_grp_id, "BackgroundMesh");
	// read nodes
	size_t node_num;
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	using Model_T2D_ME_mt_hdf5_utilities::get_node_dt_id;
	using Model_T2D_ME_mt_hdf5_utilities::NodeData;
	hid_t node_dt_id = Model_T2D_ME_mt_hdf5_utilities::get_node_dt_id();
	NodeData* nodes_data = new NodeData[node_num];
	rf.read_dataset(bg_mesh_id, "NodeData", node_num, nodes_data, node_dt_id);
	H5Tclose(node_dt_id);
	// read elements
	size_t elem_num;
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	using Model_T2D_ME_mt_hdf5_utilities::get_element_dt_id;
	using Model_T2D_ME_mt_hdf5_utilities::ElementData;
	hid_t elem_dt_id = get_element_dt_id();
	ElementData* elems_data = new ElementData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementData", elem_num, elems_data, elem_dt_id);
	H5Tclose(elem_dt_id);
	// init bg_mesh buffer
	res = bg_mesh_obj.init_from_elements(
		nodes_data,
		node_num,
		elems_data,
		elem_num,
		mesh_color);
	// get bounding box
	if (display_whole_model)
	{
		bbox.xl = nodes_data[0].x;
		bbox.xu = bbox.xl;
		bbox.yl = nodes_data[0].y;
		bbox.yu = bbox.yl;
		for (size_t n_id = 1; n_id < node_num; ++n_id)
		{
			NodeData &n = nodes_data[n_id];
			if (bbox.xl > n.x)
				bbox.xl = n.x;
			if (bbox.xu < n.x)
				bbox.xu = n.x;
			if (bbox.yl > n.y)
				bbox.yl = n.y;
			if (bbox.yu < n.y)
				bbox.yu = n.y;
		}
	}
	delete[] nodes_data;
	delete[] elems_data;
	rf.close_group(bg_mesh_id);
	if (res)
		return res;

	// init particle data
	need_mat_model_data = pfld->need_mat_model_data();
	res = data_loader.load_frame_data(
		frame_id,
		need_mat_model_data);
	if (res) return res;

	size_t pcl_num = data_loader.get_pcl_num();
	pcl_fld_mem.reserve(pcl_num * 4);
	// x coord
	double *pcl_x_data = pcl_fld_mem.get_mem();
	x_fld.extract_pcl_fld_data(pcl_x_data);
	// y coord
	double *pcl_y_data = pcl_x_data + pcl_num;
	y_fld.extract_pcl_fld_data(pcl_y_data);
	// vol
	double *pcl_vol_data = pcl_y_data + pcl_num;
	vol_fld.extract_pcl_fld_data(pcl_vol_data);
	// pcl fld
	double *pcl_fld_data = pcl_vol_data + pcl_num;
	pfld->extract_pcl_fld_data(pcl_fld_data);
	pcls_obj.init(
		pcl_num,
		pcl_x_data,
		pcl_y_data,
		pcl_vol_data,
		pcl_fld_data,
		0.5);
	
	// color map texture
	size_t color_map_texture_size;
	unsigned char* color_map_texture_data;
	if (pcl_is_mono_color)
		color_map_texture_data = color_map1.gen_1Dtexture(1, color_map_texture_size);
	else
		color_map_texture_data = color_map.gen_1Dtexture(20, color_map_texture_size);
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

	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	Rect rb_bbox;
	// rigid rect
	if (rf.has_group(frame_grp_id, "RigidRect"))
	{
		has_rr_obj = true;
		hid_t rb_grp_id = rf.open_group(frame_grp_id, "RigidRect");
		double rr_hx, rr_hy, rr_x, rr_y, rr_ang;
		rf.read_attribute(rb_grp_id, "hx", rr_hx);
		rf.read_attribute(rb_grp_id, "hy", rr_hy);
		rf.read_attribute(rb_grp_id, "x", rr_x);
		rf.read_attribute(rb_grp_id, "y", rr_y);
		rf.read_attribute(rb_grp_id, "angle", rr_ang);
		rr_obj.init(rr_x, rr_y, rr_ang, rr_hx, rr_hy, rb_color, 3.0f);
		rf.close_group(rb_grp_id);
		double cos_ang = cos(rr_ang);
		double sin_ang = sin(rr_ang);
		rr_hx *= 0.5;
		rr_hy *= 0.5;
		double xr = abs(cos_ang * rr_hx) + abs(sin_ang * rr_hy);
		double yr = abs(sin_ang * rr_hx) + abs(cos_ang * rr_hy);
		rb_bbox.xl = rr_x - xr;
		rb_bbox.xu = rr_x + xr;
		rb_bbox.yl = rr_y - yr;
		rb_bbox.yu = rr_y + yr;
		if (display_whole_model)
			bbox.envelop(rb_bbox);
	}
	rf.close_group(frame_grp_id);
	// rigid body by t2d mesh
	QVector3D navajowhite(1.0f, 0.871f, 0.678f);
	if (rf.has_group(md_data_grp_id, "RigidObjectByT2DMesh"))
	{
		has_rb_obj = true;
		hid_t rb_grp_id = rf.open_group(md_data_grp_id, "RigidObjectByT2DMesh");
		size_t edge_num;
		rf.read_attribute(rb_grp_id, "edge_num", edge_num);
		double rb_x, rb_y, rb_ang;
		rf.read_attribute(rb_grp_id, "x", rb_x);
		rf.read_attribute(rb_grp_id, "y", rb_y);
		rf.read_attribute(rb_grp_id, "ang", rb_ang);
		PointToLineDistance *pt_ln_dist = new PointToLineDistance[edge_num];
		hid_t pt_ln_dt = RigidObject_hdf5_utilities::get_pt_to_ln_dist_dt_id();
		rf.read_dataset(
			rb_grp_id,
			"pt_ln_dist",
			edge_num,
			pt_ln_dist,
			pt_ln_dt);
		H5Tclose(pt_ln_dt);
		rb_obj.init_edges(
			pt_ln_dist,
			edge_num,
			rb_x, rb_y, rb_ang,
			navajowhite);
		delete[] pt_ln_dist;
		double grid_xl, grid_yl, grid_hx, grid_hy;
		size_t grid_x_num, grid_y_num;
		rf.read_attribute(rb_grp_id, "grid_xl", grid_xl);
		rf.read_attribute(rb_grp_id, "grid_yl", grid_yl);
		rf.read_attribute(rb_grp_id, "grid_hx", grid_hx);
		rf.read_attribute(rb_grp_id, "grid_hy", grid_hy);
		rf.read_attribute(rb_grp_id, "grid_x_num", grid_x_num);
		rf.read_attribute(rb_grp_id, "grid_y_num", grid_y_num);
		Rect rmesh_bbox(grid_xl + rb_x, grid_xl + double(grid_x_num) * grid_hx + rb_x,
						grid_yl + rb_y, grid_yl + double(grid_y_num) * grid_hy + rb_y);
		if (display_whole_model)
			bbox.envelop(rmesh_bbox);
	}

	// color map display
	if (has_color_map)
	{
		color_map_obj.init(
			cm_xpos, cm_ypos, cm_ht,
			color_map,
			field_name.c_str(),
			"%8.3e",
			"../../Asset/times_new_roman.ttf");
	}

	// viewport
	double xlen = bbox.xu - bbox.xl;
	double ylen = bbox.yu - bbox.yl;
	double padding = (xlen > ylen ? xlen : ylen) * padding_ratio;
	bbox.xl -= padding;
	bbox.xu += padding;
	bbox.yl -= padding;
	bbox.yu += padding;
	set_viewport(wd, ht, GLfloat(bbox.xu - bbox.xl), GLfloat(bbox.yu - bbox.yl));

	// init shaders
	// shader_plain2D
	shader_plain2D.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_plain2D.vert");
	shader_plain2D.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_plain2D.frag");
	shader_plain2D.link();

	// shader_rigid_mesh_2D
	shader_rigid_mesh_2D.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_rigid_mesh2D.vert");
	shader_rigid_mesh_2D.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_rigid_mesh2D.frag");
	shader_rigid_mesh_2D.link();

	// shader_circles
	shader_circles.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_circles.vert");
	shader_circles.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_circles.frag");
	shader_circles.link();

	// shader for displaying character
	shader_char.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_char.vert");
	shader_char.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_char.frag");
	shader_char.link();

	// view matrix
	// mpm model view matrix
	view_mat.setToIdentity();
	view_mat.ortho(
		GLfloat(bbox.xl),
		GLfloat(bbox.xu),
		GLfloat(bbox.yl),
		GLfloat(bbox.yu),
		-1.0f, 1.0f);
	// hud view matrix
	hud_view_mat.setToIdentity();
	hud_view_mat.ortho(0.0f, GLfloat(wd) / GLfloat(ht), 0.0f, 1.0f, -1.0f, 1.0f);

	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", view_mat);

	shader_rigid_mesh_2D.bind();
	shader_rigid_mesh_2D.setUniformValue("view_mat", view_mat);

	shader_circles.bind();
	shader_circles.setUniformValue("view_mat", view_mat);
	// color map data used by shader_circles
	GLfloat color_value_lower = color_map.get_lower_bound();
	shader_circles.setUniformValue("color_value_lower", color_value_lower);
	GLfloat color_value_range = color_map.get_range();
	shader_circles.setUniformValue("color_value_range", color_value_range);
	shader_circles.setUniformValue("color_map", 0);

	shader_char.bind();
	shader_char.setUniformValue("view_mat", hud_view_mat);

	return 0;
}

void QtSceneFromHdf5_T2D_ME_mt::update_scene(size_t frame_id)
{
	data_loader.load_frame_data(frame_id, need_mat_model_data);

	size_t pcl_num = data_loader.get_pcl_num();
	pcl_fld_mem.reserve(pcl_num * 4);
	// x coord
	double *pcl_x_data = pcl_fld_mem.get_mem();
	x_fld.extract_pcl_fld_data(pcl_x_data);
	// y coord
	double *pcl_y_data = pcl_x_data + pcl_num;
	y_fld.extract_pcl_fld_data(pcl_y_data);
	// vol
	double *pcl_vol_data = pcl_y_data + pcl_num;
	vol_fld.extract_pcl_fld_data(pcl_vol_data);
	// pcl fld
	double *pcl_fld_data = pcl_vol_data + pcl_num;
	pfld->extract_pcl_fld_data(pcl_fld_data);
	pcls_obj.update(
		pcl_num,
		pcl_x_data,
		pcl_y_data,
		pcl_vol_data,
		pcl_fld_data,
		0.5);
	
	ResultFile_hdf5& rf = *res_file;
	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	// rigid rect
	if (rf.has_group(frame_grp_id, "RigidRect"))
	{
		hid_t rb_grp_id = rf.open_group(frame_grp_id, "RigidRect");
		double rr_x, rr_y, rr_ang;
		rf.read_attribute(rb_grp_id, "x", rr_x);
		rf.read_attribute(rb_grp_id, "y", rr_y);
		rf.read_attribute(rb_grp_id, "angle", rr_ang);
		rr_obj.update(rr_x, rr_y, rr_ang);
		rf.close_group(rb_grp_id);
	}
	// rigid t2d mesh
	if (rf.has_group(frame_grp_id, "RigidObjectByT2DMesh"))
	{
		hid_t rb_grp_id = rf.open_group(frame_grp_id, "RigidObjectByT2DMesh");
		double rb_x, rb_y, rb_ang;
		rf.read_attribute(rb_grp_id, "x", rb_x);
		rf.read_attribute(rb_grp_id, "y", rb_y);
		rf.read_attribute(rb_grp_id, "ang", rb_ang);
		rb_obj.update(rb_x, rb_y, rb_ang);
		rf.close_group(rb_grp_id);
	}
	//
	rf.close_group(frame_grp_id);
}
