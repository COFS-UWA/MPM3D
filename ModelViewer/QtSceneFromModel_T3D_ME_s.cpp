#include "ModelViewer_pcp.h"

#include "QtSceneFromModel_T3D_ME_s.h"

QtSceneFromModel_T3D_ME_s::QtSceneFromModel_T3D_ME_s(
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
	display_bg_mesh(true), display_pcls(true), display_pts(true),
	bg_mesh_obj(_gl), pcls_obj(_gl), pts_obj(_gl) {}

QtSceneFromModel_T3D_ME_s::~QtSceneFromModel_T3D_ME_s() {}

void QtSceneFromModel_T3D_ME_s::update_view_mat()
{
	float dist_from_obj = md_radius * view_dist_scale / sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	view_dir.normalize();
	view_pos = md_centre - dist_from_obj * view_dir;
	view_mat.setToIdentity();
	view_mat.lookAt(view_pos, md_centre, up_dir);
}

void QtSceneFromModel_T3D_ME_s::update_proj_mat()
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

void QtSceneFromModel_T3D_ME_s::update_light_pos()
{
	float dist_from_obj = md_radius * light_dist_scale
			/ sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	light_dir.normalize();
	light_pos = md_centre - dist_from_obj * light_dir;
}

int QtSceneFromModel_T3D_ME_s::set_pts_from_pcl_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_radius = radius;
	pt_num = id_num;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	Model_T3D_ME_s::Particle* pcls = model->get_pcls();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		Model_T3D_ME_s::Particle& pcl = pcls[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
		pd.z = GLfloat(pcl.z);
	}
	return 0;
}

int QtSceneFromModel_T3D_ME_s::set_pts_from_node_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_radius = radius;
	pt_num = id_num;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	Model_T3D_ME_s::Node* nodes = model->get_nodes();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		Model_T3D_ME_s::Node& n = nodes[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(n.x);
		pd.y = GLfloat(n.y);
		pd.z = GLfloat(n.z);
	}
	return 0;
}

int QtSceneFromModel_T3D_ME_s::initialize(int wd, int ht)
{
	gl.glEnable(GL_DEPTH_TEST);

	width = wd;
	height = ht;

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

	// bounding circle
	Cube mh_bbox = model->get_bounding_box();
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

	// init bg_mesh
	QVector3D gray(0.5f, 0.5f, 0.5f);
	bg_mesh_obj.init_from_edges(
		model->get_nodes(),
		model->get_node_num(),
		model->get_edges(),
		model->get_edge_num(),
		gray
		);

	// init pcls
	QVector3D moccasin(1.0f, 0.8941f, 0.7098f);
	pcls_obj.init(
		model->get_pcls(),
		model->get_pcl_num(),
		moccasin,
		0.5f
		);

	// init pts
	QVector3D red(1.0f, 0.0f, 0.0f);
	if (pts && pt_num)
		pts_obj.init(pts, pt_num, pt_radius, red);

	return 0;
}

void QtSceneFromModel_T3D_ME_s::draw()
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
}

void QtSceneFromModel_T3D_ME_s::resize(int wd, int ht)
{
	width = wd;
	height = ht;

	gl.glViewport(0, 0, width, height);

	update_proj_mat();

	shader_plain3D.bind();
	shader_plain3D.setUniformValue("proj_mat", proj_mat);

	shader_balls.bind();
	shader_balls.setUniformValue("proj_mat", proj_mat);
}
