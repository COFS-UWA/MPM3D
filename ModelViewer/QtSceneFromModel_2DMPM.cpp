#include "ModelViewer_pcp.h"

#include "QtSceneFromModel_2DMPM.h"

QtSceneFromModel_2DMPM::QtSceneFromModel_2DMPM(
	QOpenGLFunctions_3_3_Core &_gl) :
	gl(_gl), model(nullptr), pt_num(0), pts(nullptr),
	display_bg_mesh(true), display_pcls(true), display_pts(true),
	bg_mesh_obj(_gl), pcls_obj(_gl), pts_obj(_gl),
	padding_ratio(0.05f), bg_color(0.2f, 0.3f, 0.3f)
{

}

QtSceneFromModel_2DMPM::~QtSceneFromModel_2DMPM()
{

}

void QtSceneFromModel_2DMPM::set_viewport(
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

int QtSceneFromModel_2DMPM::initialize(int wd, int ht)
{
	// init shaders
	// shader_plain2D
	shader_plain_2D.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_plain2D.vert"
		);
	shader_plain_2D.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_plain2D.frag"
		);
	shader_plain_2D.link();

	// shader_circles
	shader_circle_insts.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"../../Asset/shader_circles.vert"
		);
	shader_circle_insts.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"../../Asset/shader_circles.frag"
		);
	shader_circle_insts.link();

	// view matrix
	Rect bbox = model->get_bounding_box();
	GLfloat xlen = GLfloat(bbox.xu - bbox.xl);
	GLfloat ylen = GLfloat(bbox.yu - bbox.yl);
	GLfloat padding = (xlen > ylen ? xlen : ylen) * padding_ratio;
	xl = GLfloat(bbox.xl) - padding;
	xu = GLfloat(bbox.xu) + padding;
	yl = GLfloat(bbox.yl) - padding;
	yu = GLfloat(bbox.yu) + padding;

	set_viewport(wd, ht, xu - xl, yu - yl);

	view_mat.setToIdentity();
	view_mat.ortho(xl, xu, yl, yu, -1.0f, 1.0f);
	shader_plain_2D.bind();
	shader_plain_2D.setUniformValue("view_mat", view_mat);
	shader_circle_insts.bind();
	shader_circle_insts.setUniformValue("view_mat", view_mat);

	// init bg_mesh
	QVector3D gray(0.5f, 0.5f, 0.5f);
	bg_mesh_obj.init(
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

void QtSceneFromModel_2DMPM::draw()
{
	gl.glViewport(vp_x_pos, vp_y_pos, vp_x_size, vp_y_size);

	gl.glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	gl.glClear(GL_COLOR_BUFFER_BIT);

	shader_plain_2D.bind();

	if (display_bg_mesh)
		bg_mesh_obj.draw(shader_plain_2D);

	shader_circle_insts.bind();

	if (display_pcls)
		pcls_obj.draw(shader_circle_insts);

	if (display_pts && pts && pt_num)
		pts_obj.draw(shader_circle_insts);
}

void QtSceneFromModel_2DMPM::resize(int wd, int ht)
{
	set_viewport(wd, ht, xu - xl, yu - yl);
}

int QtSceneFromModel_2DMPM::set_pts_from_pcl_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_num = id_num;
	pt_radius = radius;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	Model_T2D_ME_s::Particle* pcls = model->get_pcls();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		Model_T2D_ME_s::Particle& pcl = pcls[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
	}
	return 0;
}

int QtSceneFromModel_2DMPM::set_pts_from_node_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_num = id_num;
	pt_radius = radius;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	Model_T2D_ME_s::Node* nodes = model->get_nodes();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		Model_T2D_ME_s::Node& n = nodes[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(n.x);
		pd.y = GLfloat(n.y);
	}
	return 0;
}
