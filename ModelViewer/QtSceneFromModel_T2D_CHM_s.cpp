#include "ModelViewer_pcp.h"

#include "QtSceneFromModel_T2D_CHM_s.h"

QtSceneFromModel_T2D_CHM_s::QtSceneFromModel_T2D_CHM_s(
	QOpenGLFunctions_3_3_Core &_gl) :
	QtSceneFromModel(_gl),
	model(nullptr),
	pt_num(0), pts(nullptr),
	display_bg_mesh(true),
	display_pcls(true),
	display_rigid_circle(true),
	display_pts(true),
	bg_mesh_obj(_gl),
	pcls_obj(_gl),
	rc_obj(_gl),
	pts_obj(_gl),
	display_whole_model(true),
	padding_ratio(0.05f),
	bg_color(0.2f, 0.3f, 0.3f)
{

}

QtSceneFromModel_T2D_CHM_s::~QtSceneFromModel_T2D_CHM_s()
{

}

void QtSceneFromModel_T2D_CHM_s::set_viewport(
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

int QtSceneFromModel_T2D_CHM_s::initialize(int wd, int ht)
{
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

	// get bounding box
	if (display_whole_model)
	{
		Rect bbox = model->get_bounding_box();
		GLfloat xlen = GLfloat(bbox.xu - bbox.xl);
		GLfloat ylen = GLfloat(bbox.yu - bbox.yl);
		GLfloat padding = (xlen > ylen ? xlen : ylen) * padding_ratio;
		xl = GLfloat(bbox.xl) - padding;
		xu = GLfloat(bbox.xu) + padding;
		yl = GLfloat(bbox.yl) - padding;
		yu = GLfloat(bbox.yu) + padding;
	}

	// viewport
	set_viewport(wd, ht, xu - xl, yu - yl);

	// view matrix
	view_mat.setToIdentity();
	view_mat.ortho(xl, xu, yl, yu, -1.0f, 1.0f);
	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", view_mat);
	shader_circles.bind();
	shader_circles.setUniformValue("view_mat", view_mat);

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

	// init rigid circle
	QVector3D light_slate_blue(0.5176f, 0.4392, 1.0f);
	if (model->rigid_circle_is_valid())
	{
		RigidCircle& rc = model->get_rigid_circle();
		rc_obj.init(
			rc.get_x(),
			rc.get_y(),
			rc.get_radius(),
			light_slate_blue,
			3.0f
		);
	}

	// init pts
	QVector3D red(1.0f, 0.0f, 0.0f);
	if (pts && pt_num)
		pts_obj.init(pts, pt_num, pt_radius, red);

	return 0;
}

void QtSceneFromModel_T2D_CHM_s::draw()
{
	gl.glViewport(vp_x_pos, vp_y_pos, vp_x_size, vp_y_size);

	gl.glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	gl.glClear(GL_COLOR_BUFFER_BIT);

	shader_plain2D.bind();

	if (display_bg_mesh)
		bg_mesh_obj.draw(shader_plain2D);

	if (display_rigid_circle && model->rigid_circle_is_valid())
		rc_obj.draw(shader_plain2D);
	
	shader_circles.bind();

	if (display_pcls)
		pcls_obj.draw(shader_circles);

	if (display_pts && pts && pt_num)
		pts_obj.draw(shader_circles);
}

void QtSceneFromModel_T2D_CHM_s::resize(int wd, int ht)
{
	set_viewport(wd, ht, xu - xl, yu - yl);
}

int QtSceneFromModel_T2D_CHM_s::set_pts_from_pcl_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_num = id_num;
	pt_radius = radius;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	Model_T2D_CHM_s::Particle* pcls = model->get_pcls();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle& pcl = pcls[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
	}
	return 0;
}

int QtSceneFromModel_T2D_CHM_s::set_pts_from_node_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_num = id_num;
	pt_radius = radius;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	Model_T2D_CHM_s::Node* nodes = model->get_nodes();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		Model_T2D_CHM_s::Node& n = nodes[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(n.x);
		pd.y = GLfloat(n.y);
	}
	return 0;
}
