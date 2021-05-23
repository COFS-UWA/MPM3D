#include "ModelViewer_pcp.h"

#include "RigidObject/RigidObject_hdf5_utilities.h"
#include "Model_T3D_ME_mt_hdf5_utilities.h"
#include "QtSceneFromHdf5_T3D_ME_mt.h"

QtSceneFromHdf5_T3D_ME_mt::QtSceneFromHdf5_T3D_ME_mt(
	QOpenGLFunctions_3_3_Core &_gl) :
	QtSceneFromHdf5(_gl),
	// view direction
	fov_angle(45.0f), view_dist_scale(1.0f),
	view_dir(1.0f, 0.0f, 0.0f), up_dir(0.0f, 0.0f, 1.0f),
	bg_color(0.2f, 0.3f, 0.3f),
	// fog effect
	fog_coef(0.0f), fog_color(bg_color),
	// phong model
	light_color(1.0f, 1.0f, 1.0f),
	amb_coef(0.0f), diff_coef(1.0f),
	spec_coef(0.0f), spec_shininess(30.0f),
	light_dir(view_dir), light_dist_scale(2.0f),
	// model
	res_file(nullptr), th_id(-1),
	x_fld(data_loader), y_fld(data_loader), z_fld(data_loader),
	vol_fld(data_loader), pfld(nullptr),
	display_bg_mesh(true), bg_mesh_obj(_gl),
	display_pcls(true), pcls_obj(_gl),
	has_color_map(false), color_map_obj(_gl), color_map_texture(0),
	need_mat_model_data(false),
	display_rcy(true), has_rcy(false), rcy_obj(_gl),
	display_rco(true), has_rco(false), rco_obj(_gl),
	display_rcu(true), has_rcu(false), rcu_obj(_gl),
	display_rmesh(true), has_rmesh(false), rmesh_obj(_gl),
	pcl_is_mono_color(false)
{
	amb_coef = 0.3f;
	diff_coef = 1.0f;
	spec_coef = 0.5f;
	spec_shininess = 32.0f;
}

QtSceneFromHdf5_T3D_ME_mt::~QtSceneFromHdf5_T3D_ME_mt()
{
	clear();
	close_file();
}

void QtSceneFromHdf5_T3D_ME_mt::clear()
{
	if (color_map_texture)
	{
		gl.glDeleteTextures(1, &color_map_texture);
		color_map_texture = 0;
	}
}

void QtSceneFromHdf5_T3D_ME_mt::close_file()
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

void QtSceneFromHdf5_T3D_ME_mt::update_view_mat()
{
	float dist_from_obj = md_radius * view_dist_scale / sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	view_dir.normalize();
	view_pos = md_centre - dist_from_obj * view_dir;
	view_mat.setToIdentity();
	view_mat.lookAt(view_pos, md_centre, up_dir);
}

void QtSceneFromHdf5_T3D_ME_mt::update_proj_mat()
{
	float aspect_ratio = float(win_wd) / float(win_ht);
	proj_mat.setToIdentity();
	if (aspect_ratio >= 1.0f)
	{
		proj_mat.perspective(fov_angle, aspect_ratio, 0.01f, 10000.0f);
	}
	else
	{
		float fov_angle2 = atan(tan(fov_angle / 180.0 * 3.14159265359) / aspect_ratio) / 3.14159265359 * 180.0;
		proj_mat.perspective(fov_angle2, aspect_ratio, 0.01f, 10000.0f);
	}
}

void QtSceneFromHdf5_T3D_ME_mt::update_light_pos()
{
	const float dist_from_obj = md_radius * light_dist_scale
		/ (float)sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	light_dir.normalize();
	light_pos = md_centre - dist_from_obj * light_dir;
}

// ================== animation ==================
int QtSceneFromHdf5_T3D_ME_mt::set_res_file(
	ResultFile_hdf5& rf,
	const char* th_name,
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
	if (!z_fld.validate_data_type())
		throw std::exception("z field not found in hdf5 field data.");
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

size_t QtSceneFromHdf5_T3D_ME_mt::get_frame_num()
{
	ResultFile_hdf5& rf = *res_file;
	size_t frame_num;
	rf.read_attribute(th_id, "output_num", frame_num);
	return frame_num;
}

double QtSceneFromHdf5_T3D_ME_mt::get_frame_time(size_t frame_id)
{
	ResultFile_hdf5& rf = *res_file;
	double frame_time;
	char frame_name[50];
	snprintf(frame_name, sizeof(frame_name), "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	rf.read_attribute(frame_grp_id, "total_time", frame_time);
	rf.close_group(frame_grp_id);
	return frame_time;
}

int QtSceneFromHdf5_T3D_ME_mt::init_scene(int wd, int ht, size_t frame_id)
{
	using namespace Model_hdf5_utilities;

	if (!res_file || !res_file->is_open())
		return -1;

	win_wd = wd; win_ht = ht;

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
	rf.read_dataset(bg_mesh_id, "NodeData", node_num, nodes_data, node_dt_id);
	H5Tclose(node_dt_id);
	// read elements
	size_t elem_num;
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	hid_t elem_dt_id = get_ed_3d_dt_id();
	Elem3DData* elems_data = new Elem3DData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementData", elem_num, elems_data, elem_dt_id);
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

	// init particle data
	need_mat_model_data = pfld->need_mat_model_data();
	res = data_loader.load_frame_data(
		frame_id,
		need_mat_model_data);
	if (res)
		return res;

	size_t pcl_num = data_loader.get_pcl_num();
	if (pcl_num)
	{
		pcl_fld_mem.reserve(pcl_num * 5);
		// x coord
		double* pcl_x_data = pcl_fld_mem.get_mem();
		x_fld.extract_pcl_fld_data(pcl_x_data);
		// y coord
		double* pcl_y_data = pcl_x_data + pcl_num;
		y_fld.extract_pcl_fld_data(pcl_y_data);
		// z coord
		double* pcl_z_data = pcl_y_data + pcl_num;
		z_fld.extract_pcl_fld_data(pcl_z_data);
		// vol
		double* pcl_vol_data = pcl_z_data + pcl_num;
		vol_fld.extract_pcl_fld_data(pcl_vol_data);
		// pcl fld
		double* pcl_fld_data = pcl_vol_data + pcl_num;
		pfld->extract_pcl_fld_data(pcl_fld_data);
		pcls_obj.init(
			pcl_num,
			pcl_x_data,
			pcl_y_data,
			pcl_z_data,
			pcl_vol_data,
			pcl_fld_data,
			0.5
		);
	}

	init_rigid_objects_buffer(mh_bbox, rf);
	if (!init_color_map_texture())
		return -2;
	init_shaders(mh_bbox, wd, ht);
	return 0;
}

void QtSceneFromHdf5_T3D_ME_mt::update_scene(size_t frame_id)
{
	ResultFile_hdf5 &rf = *res_file;

	data_loader.load_frame_data(frame_id, need_mat_model_data);
	size_t pcl_num = data_loader.get_pcl_num();
	if (pcl_num)
	{
		pcl_fld_mem.reserve(pcl_num * 5);
		// x coord
		double* pcl_x_data = pcl_fld_mem.get_mem();
		x_fld.extract_pcl_fld_data(pcl_x_data);
		// y coord
		double* pcl_y_data = pcl_x_data + pcl_num;
		y_fld.extract_pcl_fld_data(pcl_y_data);
		// z coord
		double* pcl_z_data = pcl_y_data + pcl_num;
		z_fld.extract_pcl_fld_data(pcl_z_data);
		// vol
		double* pcl_vol_data = pcl_z_data + pcl_num;
		vol_fld.extract_pcl_fld_data(pcl_vol_data);
		// pcl fld
		double* pcl_fld_data = pcl_vol_data + pcl_num;
		pfld->extract_pcl_fld_data(pcl_fld_data);
		pcls_obj.update(
			pcl_num,
			pcl_x_data,
			pcl_y_data,
			pcl_z_data,
			pcl_vol_data,
			pcl_fld_data,
			0.5
			);
	}

	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);	
	update_rigid_objects_buffer(frame_grp_id, rf);
	rf.close_group(frame_grp_id);
}

void QtSceneFromHdf5_T3D_ME_mt::draw()
{
	gl.glEnable(GL_DEPTH_TEST);

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
	
	if (display_rcy && has_rcy)
		rcy_obj.draw(shader_rigid_mesh);

	if (display_rco && has_rco)
		rco_obj.draw(shader_rigid_mesh);

	if (display_rcu && has_rcu)
		rcu_obj.draw(shader_rigid_mesh);

	if (display_rmesh && has_rmesh)
		rmesh_obj.draw(shader_rigid_mesh);

	// all 3d objects should be drawn before this
	gl.glDisable(GL_DEPTH_TEST);

	if (has_color_map)
		color_map_obj.draw(shader_plain2D, shader_char);
}

void QtSceneFromHdf5_T3D_ME_mt::resize(int wd, int ht)
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
	hud_view_mat.ortho(0.0f, GLfloat(wd) / GLfloat(ht), 0.0f, 1.0f, -1.0f, 1.0f);
	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", hud_view_mat);
	shader_char.bind();
	shader_char.setUniformValue("view_mat", hud_view_mat);
}

void QtSceneFromHdf5_T3D_ME_mt::init_rigid_objects_buffer(
	Cube &mh_bbox,
	ResultFile_hdf5& rf
	)
{
	QVector3D navajowhite(1.0f, 0.871f, 0.678f);
	hid_t md_data_grp_id = rf.get_model_data_grp_id();
	
	// rigid cylinder
	if (rf.has_group(md_data_grp_id, "RigidCylinder"))
	{
		hid_t rcy_id = rf.open_group(md_data_grp_id, "RigidCylinder");
		double rcy_x, rcy_y, rcy_z, rcy_h, rcy_r;
		rf.read_attribute(rcy_id, "x", rcy_x);
		rf.read_attribute(rcy_id, "y", rcy_y);
		rf.read_attribute(rcy_id, "z", rcy_z);
		rf.read_attribute(rcy_id, "h", rcy_h);
		rf.read_attribute(rcy_id, "r", rcy_r);
		rcy_obj.init(rcy_x, rcy_y, rcy_z, rcy_h, rcy_r, navajowhite);
		Cube rcy_bbox(rcy_x - rcy_r, rcy_x + rcy_r,
					  rcy_y - rcy_r, rcy_y + rcy_r,
					  rcy_z - rcy_h * 0.5,
					  rcy_z + rcy_h * 0.5);
		mh_bbox.envelop(rcy_bbox);
		rf.close_group(rcy_id);
		has_rcy = true;
	}

	// rigid cone
	if (rf.has_group(md_data_grp_id, "RigidCone"))
	{
		hid_t rco_id = rf.open_group(md_data_grp_id, "RigidCone");
		double rco_x, rco_y, rco_z;
		double rco_h_tip, rco_h_shaft, rco_r;
		rf.read_attribute(rco_id, "x", rco_x);
		rf.read_attribute(rco_id, "y", rco_y);
		rf.read_attribute(rco_id, "z", rco_z);
		rf.read_attribute(rco_id, "h_tip", rco_h_tip);
		rf.read_attribute(rco_id, "h_shaft", rco_h_shaft);
		rf.read_attribute(rco_id, "r", rco_r);
		rco_obj.init(rco_x, rco_y, rco_z,
			rco_h_tip, rco_h_shaft, rco_r,
			navajowhite);
		Cube rco_bbox(rco_x - rco_r, rco_x + rco_r,
			rco_y - rco_r, rco_y + rco_r,
			rco_z - rco_h_tip, rco_z + rco_h_shaft);
		mh_bbox.envelop(rco_bbox);
		rf.close_group(rco_id);
		has_rco = true;
	}

	// rigid cube
	if (rf.has_group(md_data_grp_id, "RigidCube"))
	{
		hid_t rcu_id = rf.open_group(md_data_grp_id, "RigidCube");
		double rcu_x, rcu_y, rcu_z;
		double rcu_hx, rcu_hy, rcu_hz;
		rf.read_attribute(rcu_id, "x", rcu_x);
		rf.read_attribute(rcu_id, "y", rcu_y);
		rf.read_attribute(rcu_id, "z", rcu_z);
		rf.read_attribute(rcu_id, "hx", rcu_hx);
		rf.read_attribute(rcu_id, "hy", rcu_hy);
		rf.read_attribute(rcu_id, "hz", rcu_hz);
		rcu_obj.init(rcu_x, rcu_y, rcu_z,
			rcu_hx, rcu_hy, rcu_hz, navajowhite);
		Cube rcu_bbox(rcu_x - 0.5 * rcu_hx, rcu_x + 0.5 * rcu_hy,
			rcu_y - 0.5 * rcu_hy, rcu_y + 0.5 * rcu_hy,
			rcu_z - 0.5 * rcu_hz, rcu_z + 0.5 * rcu_hz);
		mh_bbox.envelop(rcu_bbox);
		rf.close_group(rcu_id);
		has_rcu = true;
	}

	// rigid t3d mesh
	if (rf.has_group(md_data_grp_id, "RigidObjectByT3DMesh"))
	{
		hid_t rmesh_id = rf.open_group(md_data_grp_id, "RigidObjectByT3DMesh");
		size_t face_num;
		rf.read_attribute(rmesh_id, "face_num", face_num);
		if (face_num)
		{
			Point3D rm_pos;
			rf.read_attribute(rmesh_id, "x", rm_pos.x);
			rf.read_attribute(rmesh_id, "y", rm_pos.y);
			rf.read_attribute(rmesh_id, "z", rm_pos.z);
			Vector3D rm_ang;
			rf.read_attribute(rmesh_id, "x_ang", rm_ang.x);
			rf.read_attribute(rmesh_id, "y_ang", rm_ang.y);
			rf.read_attribute(rmesh_id, "z_ang", rm_ang.z);
			PointToTriangleDistance* pt_tri_dist
				= new PointToTriangleDistance[face_num];
			hid_t pt_tri_dt = RigidObject_hdf5_utilities::get_pt_to_tri_dist_dt_id();
			rf.read_dataset(
				rmesh_id,
				"pt_tri_dist",
				face_num,
				pt_tri_dist,
				pt_tri_dt
				);
			H5Tclose(pt_tri_dt);
			rmesh_obj.init_faces(
				pt_tri_dist,
				face_num,
				rm_pos,
				rm_ang,
				navajowhite
				);
			delete[] pt_tri_dist;
			double grid_xl, grid_yl, grid_zl;
			double grid_hx, grid_hy, grid_hz;
			size_t grid_x_num, grid_y_num, grid_z_num;
			rf.read_attribute(rmesh_id, "grid_xl", grid_xl);
			rf.read_attribute(rmesh_id, "grid_yl", grid_yl);
			rf.read_attribute(rmesh_id, "grid_zl", grid_zl);
			rf.read_attribute(rmesh_id, "grid_hx", grid_hx);
			rf.read_attribute(rmesh_id, "grid_hy", grid_hy);
			rf.read_attribute(rmesh_id, "grid_hz", grid_hz);
			rf.read_attribute(rmesh_id, "grid_x_num", grid_x_num);
			rf.read_attribute(rmesh_id, "grid_y_num", grid_y_num);
			rf.read_attribute(rmesh_id, "grid_z_num", grid_z_num);
			Cube rmesh_bbox(grid_xl + rm_pos.x, grid_xl + double(grid_x_num) * grid_hx + rm_pos.x,
							grid_yl + rm_pos.y, grid_yl + double(grid_y_num) * grid_hy + rm_pos.y,
							grid_zl + rm_pos.z, grid_zl + double(grid_z_num) * grid_hz + rm_pos.z);
			mh_bbox.envelop(rmesh_bbox);
			rf.close_group(rmesh_id);
			has_rmesh = true;
		}
	}
}

bool QtSceneFromHdf5_T3D_ME_mt::init_color_map_texture()
{
	// color map texture
	size_t color_map_texture_size;
	unsigned char* color_map_texture_data;
	if (pcl_is_mono_color)
		color_map_texture_data = color_map1.gen_1Dtexture(1, color_map_texture_size);
	else
		color_map_texture_data = color_map.gen_1Dtexture(20, color_map_texture_size);
	if (!color_map_texture_data || !color_map_texture_size)
		return false;
	gl.glGenTextures(1, &color_map_texture);
	if (color_map_texture == 0)
		return false;
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

	// color_map legend display
	if (has_color_map)
	{
		color_map_obj.init(
			cm_xpos, cm_ypos, cm_ht,
			color_map,
			field_name.c_str(),
			"%11.6e",
			"../../Asset/times_new_roman.ttf"
			);
	}

	return true;
}

void QtSceneFromHdf5_T3D_ME_mt::init_shaders(
	Cube& mh_bbox,
	int wd,
	int ht
	)
{
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

	// shader_rigid_mesh
	shader_rigid_mesh.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_rigid_mesh.vert"
	);
	shader_rigid_mesh.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_rigid_mesh.frag"
	);
	shader_rigid_mesh.link();

	md_centre.setX(float(mh_bbox.xl + mh_bbox.xu) * 0.5f);
	md_centre.setY(float(mh_bbox.yl + mh_bbox.yu) * 0.5f);
	md_centre.setZ(float(mh_bbox.zl + mh_bbox.zu) * 0.5f);
	float dx = mh_bbox.xu - mh_bbox.xl;
	float dy = mh_bbox.yu - mh_bbox.yl;
	float dz = mh_bbox.zu - mh_bbox.zl;
	md_radius = sqrt(dx * dx + dy * dy + dz * dz) * 0.5f;

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

	GLfloat color_value_lower = color_map.get_lower_bound();
	shader_balls.setUniformValue("color_value_lower", color_value_lower);
	GLfloat color_value_range = color_map.get_range();
	shader_balls.setUniformValue("color_value_range", color_value_range);
	shader_balls.setUniformValue("color_map", 0);

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

	// color map
	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", hud_view_mat);

	shader_char.bind();
	shader_char.setUniformValue("view_mat", hud_view_mat);

	// rigid mesh shader
	shader_rigid_mesh.bind();
	shader_rigid_mesh.setUniformValue("view_mat", view_mat);
	shader_rigid_mesh.setUniformValue("proj_mat", proj_mat);

	shader_rigid_mesh.setUniformValue("view_pos", view_pos);

	// fog effect
	shader_rigid_mesh.setUniformValue("fog_coef", fog_coef);
	shader_rigid_mesh.setUniformValue("fog_color", fog_color);

	// phong model parameters
	shader_rigid_mesh.setUniformValue("light_pos", light_pos);
	shader_rigid_mesh.setUniformValue("light_color", light_color);
	shader_rigid_mesh.setUniformValue("amb_coef", amb_coef);
	shader_rigid_mesh.setUniformValue("diff_coef", diff_coef);
	shader_rigid_mesh.setUniformValue("spec_coef", spec_coef);
	shader_rigid_mesh.setUniformValue("spec_shininess", spec_shininess);
}

void QtSceneFromHdf5_T3D_ME_mt::update_rigid_objects_buffer(
	hid_t frame_grp_id,
	ResultFile_hdf5& rf
	)
{
	if (rf.has_group(frame_grp_id, "RigidCylinder"))
	{
		hid_t rcy_id = rf.open_group(frame_grp_id, "RigidCylinder");
		double rcy_x, rcy_y, rcy_z;
		rf.read_attribute(rcy_id, "x", rcy_x);
		rf.read_attribute(rcy_id, "y", rcy_y);
		rf.read_attribute(rcy_id, "z", rcy_z);
		rcy_obj.update(rcy_x, rcy_y, rcy_z);
		rf.close_group(rcy_id);
	}

	if (rf.has_group(frame_grp_id, "RigidCone"))
	{
		hid_t rco_id = rf.open_group(frame_grp_id, "RigidCone");
		double rco_x, rco_y, rco_z;
		rf.read_attribute(rco_id, "x", rco_x);
		rf.read_attribute(rco_id, "y", rco_y);
		rf.read_attribute(rco_id, "z", rco_z);
		rco_obj.update(rco_x, rco_y, rco_z);
		rf.close_group(rco_id);
	}

	if (rf.has_group(frame_grp_id, "RigidCube"))
	{
		hid_t rcu_id = rf.open_group(frame_grp_id, "RigidCube");
		double rcu_x, rcu_y, rcu_z;
		rf.read_attribute(rcu_id, "x", rcu_x);
		rf.read_attribute(rcu_id, "y", rcu_y);
		rf.read_attribute(rcu_id, "z", rcu_z);
		rcu_obj.update(rcu_x, rcu_y, rcu_z);
		rf.close_group(rcu_id);
	}

	if (rf.has_group(frame_grp_id, "RigidObjectByT3DMesh"))
	{
		hid_t rmesh_id = rf.open_group(frame_grp_id, "RigidObjectByT3DMesh");
		Point3D rm_cen;
		Vector3D rm_ang;
		rf.read_attribute(rmesh_id, "x", rm_cen.x);
		rf.read_attribute(rmesh_id, "y", rm_cen.y);
		rf.read_attribute(rmesh_id, "z", rm_cen.z);
		rf.read_attribute(rmesh_id, "x_ang", rm_ang.x);
		rf.read_attribute(rmesh_id, "y_ang", rm_ang.y);
		rf.read_attribute(rmesh_id, "z_ang", rm_ang.z);
		rmesh_obj.update(rm_cen, rm_ang);
		rf.close_group(rmesh_id);
	}
}
