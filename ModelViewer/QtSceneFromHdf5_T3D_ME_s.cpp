#include "ModelViewer_pcp.h"

#include "Model_hdf5_utilities.h"
#include "QtSceneFromHdf5_T3D_ME_s.h"

QtSceneFromHdf5_T3D_ME_s::QtSceneFromHdf5_T3D_ME_s(
	QOpenGLFunctions_3_3_Core &_gl) :
	QtSceneFromHdf5(_gl), res_file(nullptr),
	frame_grp_id(-1), th_id(-1), pcl_dt_id(),
	display_bg_mesh(true), display_pcls(true),
	bg_mesh_obj(_gl), pcls_obj(_gl),
	has_color_map(false), color_map_obj(_gl),
	color_map_texture(0),
	bg_color(0.2f, 0.3f, 0.3f)
{

}

QtSceneFromHdf5_T3D_ME_s::~QtSceneFromHdf5_T3D_ME_s()
{
	close_file();
	clear();
}

void QtSceneFromHdf5_T3D_ME_s::close_file()
{
	if (!res_file)
		return;
	ResultFile_hdf5 &rf = *res_file;
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

void QtSceneFromHdf5_T3D_ME_s::clear()
{
	if (color_map_texture)
	{
		gl.glDeleteTextures(1, &color_map_texture);
		color_map_texture = 0;
	}
}

void QtSceneFromHdf5_T3D_ME_s::draw()
{
	gl.glViewport(0, 0, win_wd, win_ht);

	gl.glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	gl.glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (display_bg_mesh)
	{
		shader_plain3D.bind();
		bg_mesh_obj.draw(shader_plain3D);
	}

	if (display_pcls)
	{
		shader_balls.bind();
		pcls_obj.draw(shader_balls);
	}

	if (has_color_map)
		color_map_obj.draw(shader_plain2D, shader_char);
}

void QtSceneFromHdf5_T3D_ME_s::resize(int wd, int ht)
{
	win_wd = wd;
	win_ht = ht;

	gl.glViewport(0, 0, win_wd, win_ht);

	update_proj_mat();
	shader_plain3D.bind();
	shader_plain3D.setUniformValue("proj_mat", proj_mat);
	shader_balls.bind();
	shader_balls.setUniformValue("proj_mat", proj_mat);

	hud_view_mat.setToIdentity();
	hud_view_mat.ortho(0.0f, GLfloat(wd)/GLfloat(ht), 0.0f, 1.0f, -1.0f, 1.0f);
	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", hud_view_mat);
	shader_char.bind();
	shader_char.setUniformValue("view_mat", hud_view_mat);
}

// ================== animation ==================
int QtSceneFromHdf5_T3D_ME_s::set_res_file(
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

size_t QtSceneFromHdf5_T3D_ME_s::get_frame_num()
{
	ResultFile_hdf5 &rf = *res_file;
	size_t frame_num;
	rf.read_attribute(th_id, "output_num", frame_num);
	return frame_num;
}

double QtSceneFromHdf5_T3D_ME_s::get_frame_time(size_t frame_id)
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

int QtSceneFromHdf5_T3D_ME_s::init_scene(int wd, int ht, size_t frame_id)
{
	using namespace Model_hdf5_utilities;

	if (!res_file || !res_file->is_open())
		return -1;

	win_wd = wd;
	win_ht = ht;

	int res;
	ResultFile_hdf5& rf = *res_file;

	// write model data to view
	hid_t md_data_grp_id = rf.get_model_data_grp_id();
	hid_t bg_mesh_id = rf.open_group(md_data_grp_id, "BackgroundMesh");
	// read nodes
	size_t node_num;
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	hid_t node_dt_id = get_nd_3d_dt_id();
	Node3DData* nodes_data = new Node3DData[node_num];
	rf.read_dataset(bg_mesh_id, "NodeCoordinate", node_num, nodes_data, node_dt_id);
	H5Tclose(node_dt_id);
	// read elements
	size_t elem_num;
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	hid_t elem_dt_id = get_ed_3d_dt_id();
	Elem3DData* elems_data = new Elem3DData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementTopology", elem_num, elems_data, elem_dt_id);
	H5Tclose(elem_dt_id);
	// get bounding box
	Cube mh_bbox;
	mh_bbox.xl = nodes_data[0].x;
	mh_bbox.xu = mh_bbox.xl;
	mh_bbox.yl = nodes_data[0].y;
	mh_bbox.yu = mh_bbox.yl;
	mh_bbox.zl = nodes_data[0].z;
	mh_bbox.zu = mh_bbox.zl;
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		Node3DData& n = nodes_data[n_id];
		if (mh_bbox.xl > n.x)
			mh_bbox.xl = n.x;
		if (mh_bbox.xu < n.x)
			mh_bbox.xu = n.x;
		if (mh_bbox.yl > n.y)
			mh_bbox.yl = n.y;
		if (mh_bbox.yu < n.y)
			mh_bbox.yu = n.y;
		if (mh_bbox.zl > n.z)
			mh_bbox.zl = n.z;
		if (mh_bbox.zu < n.z)
			mh_bbox.zu = n.z;
	}
	md_centre.setX(float(mh_bbox.xl + mh_bbox.xu) * 0.5f);
	md_centre.setY(float(mh_bbox.yl + mh_bbox.yu) * 0.5f);
	md_centre.setZ(float(mh_bbox.zl + mh_bbox.zu) * 0.5f);
	float dx = mh_bbox.xu - mh_bbox.xl;
	float dy = mh_bbox.yu - mh_bbox.yl;
	float dz = mh_bbox.zu - mh_bbox.zl;
	md_radius = sqrt(dx * dx + dy * dy + dz * dz) * 0.5f;
	// init bg_mesh buffer
	QVector3D gray(0.5f, 0.5f, 0.5f);
	res = bg_mesh_obj.init_from_elements(
		nodes_data,
		node_num,
		elems_data,
		elem_num,
		gray
		);
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
	unsigned char* color_map_texture_data = color_map.gen_1Dtexture(20, color_map_texture_size);
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

	// complete reading hdf5
	rf.close_group(frame_grp_id);

	// color map display
	if (has_color_map)
	{
		color_map_obj.init(
			cm_xpos, cm_ypos, cm_ht,
			color_map,
			field_name.c_str(),
			"%8.3e",
			"../../Asset/times_new_roman.ttf"
			);
	}

	// init shaders
	// shader_plain3D
	shader_plain3D.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_plain3D.vert"
	);
	shader_plain3D.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_plain3D.frag"
	);
	shader_plain3D.link();
	// shader_balls
	shader_balls.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_balls.vert"
	);
	shader_balls.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_balls.frag"
	);
	shader_balls.link();
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
	// shader_char
	shader_char.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_char.vert"
	);
	shader_char.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_char.frag"
	);
	shader_char.link();

	update_view_mat();
	update_proj_mat();
	update_light_pos();
	hud_view_mat.setToIdentity();
	hud_view_mat.ortho(0.0f, GLfloat(wd) / GLfloat(ht), 0.0f, 1.0f, -1.0f, 1.0f);

	// plain3D shader
	shader_plain3D.bind();
	shader_plain3D.setUniformValue("view_mat", view_mat);
	shader_plain3D.setUniformValue("proj_mat", proj_mat);

	// ball shader
	shader_balls.bind();
	shader_balls.setUniformValue("view_mat", view_mat);
	shader_balls.setUniformValue("proj_mat", proj_mat);

	shader_balls.setUniformValue("view_pos", view_pos);

	// fog effect
	shader_balls.setUniformValue("fog_coef", fog_coef);
	shader_balls.setUniformValue("fog_color", fog_color);

	// phong model parameters
	shader_balls.setUniformValue("light_pos", light_pos);
	shader_balls.setUniformValue("light_color", light_color);
	shader_balls.setUniformValue("amb_coef", amb_coef);
	shader_balls.setUniformValue("diff_coef", diff_coef);
	shader_balls.setUniformValue("spec_coef", spec_coef);
	shader_balls.setUniformValue("spec_shininess", spec_shininess);

	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", hud_view_mat);

	shader_char.bind();
	shader_char.setUniformValue("view_mat", hud_view_mat);

	return 0;
}

void QtSceneFromHdf5_T3D_ME_s::update_scene(size_t frame_id)
{
	ResultFile_hdf5 &rf = *res_file;
	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	
	// pcl data
	hid_t pcl_data_id = rf.open_group(frame_grp_id, "ParticleData");
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	pcls_data_mem.reserve(pcl_size * pcl_num);
	char* pcls_data = pcls_data_mem.get_mem();
	rf.read_dataset(pcl_data_id, "field", pcl_num, (void*)pcls_data, pcl_dt_id);
	rf.close_group(pcl_data_id);
	// update pcl gl buffer
	pcls_obj.update<double>(
		pcls_data, pcl_size, pcl_num,
		pcl_x_off, pcl_y_off,
		pcl_vol_off, 0.5, pcl_fld_off
		);

	rf.close_group(frame_grp_id);
}
