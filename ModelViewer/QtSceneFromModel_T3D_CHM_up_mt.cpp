#include "ModelViewer_pcp.h"

#include "QtSceneFromModel_T3D_CHM_up_mt.h"

QtSceneFromModel_T3D_CHM_up_mt::QtSceneFromModel_T3D_CHM_up_mt(
	QOpenGLFunctions_3_3_Core &_gl) :
	QtSceneFromModel(_gl),
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
	// model data
	model(nullptr), pt_num(0), pts(nullptr),
	display_bg_mesh(true), bg_mesh_obj(_gl),
	display_pcls(true), pcls_obj(_gl),
	display_pts(true), pts_obj(_gl),
	display_rmesh_obj(true), has_rmesh_obj(false), rmesh_obj(_gl) {}

QtSceneFromModel_T3D_CHM_up_mt::~QtSceneFromModel_T3D_CHM_up_mt() {}

void QtSceneFromModel_T3D_CHM_up_mt::update_view_mat()
{
	float dist_from_obj = md_radius * view_dist_scale
		/ sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	view_dir.normalize();
	view_pos = md_centre - dist_from_obj * view_dir;
	view_mat.setToIdentity();
	view_mat.lookAt(view_pos, md_centre, up_dir);
}

void QtSceneFromModel_T3D_CHM_up_mt::update_proj_mat()
{
	float aspect_ratio = float(width) / float(height);
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

void QtSceneFromModel_T3D_CHM_up_mt::update_light_pos()
{
	float dist_from_obj = md_radius * light_dist_scale
		/ sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	light_dir.normalize();
	light_pos = md_centre - dist_from_obj * light_dir;
}

int QtSceneFromModel_T3D_CHM_up_mt::set_pts_from_pcl_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_radius = radius;
	pt_num = id_num;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	const Position *pcl_pos = model->get_pcl_pos();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		const Position &pcl = pcl_pos[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
		pd.z = GLfloat(pcl.z);
	}
	return 0;
}

int QtSceneFromModel_T3D_CHM_up_mt::set_pts_from_node_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_radius = radius;
	pt_num = id_num;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	const Position* node_pos = model->get_node_pos();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		const Position &n = node_pos[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(n.x);
		pd.y = GLfloat(n.y);
		pd.z = GLfloat(n.z);
	}
	return 0;
}

int QtSceneFromModel_T3D_CHM_up_mt::set_pts_from_drained_bc(float radius)
{
	Model_T3D_CHM_up_mt& md = *model;
	size_t node_num = md.get_node_num();
	const Model_T3D_CHM_up_mt::NodeHasVBC* nhv = md.get_node_has_vbc();
	const Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].is_drained)
		{
			const Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pd_tmp.z = n.z;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T3D_CHM_up_mt::set_pts_from_vx_bc(float radius)
{
	Model_T3D_CHM_up_mt& md = *model;
	size_t node_num = md.get_node_num();
	const Model_T3D_CHM_up_mt::NodeHasVBC* nhv = md.get_node_has_vbc();
	const Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].has_vx_bc)
		{
			const Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pd_tmp.z = n.z;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T3D_CHM_up_mt::set_pts_from_vy_bc(float radius)
{
	Model_T3D_CHM_up_mt& md = *model;
	size_t node_num = md.get_node_num();
	const Model_T3D_CHM_up_mt::NodeHasVBC* nhv = md.get_node_has_vbc();
	const Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].has_vy_bc)
		{
			const Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pd_tmp.z = n.z;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T3D_CHM_up_mt::set_pts_from_vz_bc(float radius)
{
	Model_T3D_CHM_up_mt& md = *model;
	size_t node_num = md.get_node_num();
	const Model_T3D_CHM_up_mt::NodeHasVBC* nhv = md.get_node_has_vbc();
	const Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].has_vz_bc)
		{
			const Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pd_tmp.z = n.z;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T3D_CHM_up_mt::set_pts_from_vec_bc(float radius)
{
	Model_T3D_CHM_up_mt& md = *model;
	size_t node_num = md.get_node_num();
	const Model_T3D_CHM_up_mt::NodeVBCVec* nvec = md.get_node_vbc_vec_s();
	const Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nvec[n_id].x != 0.0 || nvec[n_id].y != 0.0 ||
			nvec[n_id].z != 0.0)
		{
			const Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pd_tmp.z = n.z;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T3D_CHM_up_mt::initialize(int wd, int ht)
{
	gl.glEnable(GL_DEPTH_TEST);
	width = wd; height = ht;

	init_shaders();

	// init bg_mesh
	if (display_bg_mesh)
	{
		QVector3D gray(0.5f, 0.5f, 0.5f);
		bg_mesh_obj.init_from_elements(
			model->get_node_pos(),
			model->get_node_num(),
			model->get_elem_node_index(),
			model->get_elem_num(),
			gray);
	}

	// init pcls
	QVector3D moccasin(1.0f, 0.8941f, 0.7098f);
	auto* pcl_index = model->get_pcl_index0();
	auto* pcl_pos = model->get_pcl_pos();
	const size_t pcl_num = model->get_pcl_num();
	Model_T3D_CHM_up_mt::Position* pcl_pos1
		= (Model_T3D_CHM_up_mt::Position*)::operator new(sizeof(Model_T3D_CHM_up_mt::Position) * pcl_num);
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		const size_t ori_p_id = pcl_index[p_id];
		pcl_pos1[p_id].x = pcl_pos[ori_p_id].x;
		pcl_pos1[p_id].y = pcl_pos[ori_p_id].y;
		pcl_pos1[p_id].z = pcl_pos[ori_p_id].z;
	}
	pcls_obj.init(
		pcl_pos1,
		model->get_pcl_vol(),
		pcl_num, moccasin, 0.5f);
	::operator delete ((void*)pcl_pos1);

	// init pts
	init_pts_buffer();
	// init rigid object
	init_rigid_objects_buffer();
	return 0;
}

void QtSceneFromModel_T3D_CHM_up_mt::draw()
{
	gl.glViewport(0, 0, width, height);

	gl.glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	gl.glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (display_bg_mesh)
	{
		shader_plain3D.bind();
		bg_mesh_obj.draw(shader_plain3D);
	}

	shader_balls.bind();

	if (display_pcls)
		pcls_obj.draw(shader_balls);
	
	if (display_pts)
		pts_obj.draw(shader_balls);

	shader_rigid_mesh.bind();

	if (display_rmesh_obj && has_rmesh_obj)
		rmesh_obj.draw(shader_rigid_mesh);
}

void QtSceneFromModel_T3D_CHM_up_mt::resize(int wd, int ht)
{
	width = wd;
	height = ht;

	gl.glViewport(0, 0, width, height);

	update_proj_mat();

	shader_plain3D.bind();
	shader_plain3D.setUniformValue("proj_mat", proj_mat);

	shader_balls.bind();
	shader_balls.setUniformValue("proj_mat", proj_mat);

	shader_rigid_mesh.bind();
	shader_rigid_mesh.setUniformValue("proj_mat", proj_mat);
}

void QtSceneFromModel_T3D_CHM_up_mt::init_shaders()
{
	shader_plain3D.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_plain3D.vert"
	);
	shader_plain3D.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_plain3D.frag"
	);
	shader_plain3D.link();

	shader_balls.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_balls.vert"
	);
	shader_balls.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_balls.frag"
	);
	shader_balls.link();

	shader_rigid_mesh.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_rigid_mesh.vert"
	);
	shader_rigid_mesh.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_rigid_mesh.frag"
	);
	shader_rigid_mesh.link();

	// bounding box
	Cube mh_bbox = model->get_mesh_bbox();
	if (model->has_t3d_rigid_mesh())
	{
		RigidObjectByT3DMesh &rmesh = model->get_t3d_rigid_mesh();
		Cube rmesh_bbox;
		rmesh.get_bbox(rmesh_bbox);
		mh_bbox.envelop(rmesh_bbox);
	}
	md_centre.setX(float(mh_bbox.xl + mh_bbox.xu) * 0.5f);
	md_centre.setY(float(mh_bbox.yl + mh_bbox.yu) * 0.5f);
	md_centre.setZ(float(mh_bbox.zl + mh_bbox.zu) * 0.5f);
	const float dx = mh_bbox.xu - mh_bbox.xl;
	const float dy = mh_bbox.yu - mh_bbox.yl;
	const float dz = mh_bbox.zu - mh_bbox.zl;
	md_radius = sqrt(dx * dx + dy * dy + dz * dz) * 0.5f;

	update_view_mat();
	update_proj_mat();
	update_light_pos();

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

	// phong model
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

void QtSceneFromModel_T3D_CHM_up_mt::init_pts_buffer()
{
	QVector3D red(1.0f, 0.0f, 0.0f);
	if (pts && pt_num)
		pts_obj.init(pts, pt_num, pt_radius, red);
}

void QtSceneFromModel_T3D_CHM_up_mt::init_rigid_objects_buffer()
{
	QVector3D navajowhite(1.0f, 0.871f, 0.678f);

	if (model->has_t3d_rigid_mesh())
	{
		const RigidObjectByT3DMesh& rmesh = model->get_t3d_rigid_mesh();
		const Point3D &cen = rmesh.get_pos();
		rmesh_obj.init_faces(rmesh,	navajowhite);
		has_rmesh_obj = true;
	}
}
