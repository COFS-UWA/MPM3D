#include "ModelViewer_pcp.h"

#include "QtSceneFromModel_T2D_CHM_mt.h"

QtSceneFromModel_T2D_CHM_mt::QtSceneFromModel_T2D_CHM_mt(
	QOpenGLFunctions_3_3_Core &_gl) :
	QtSceneFromModel(_gl),
	model(nullptr),
	pt_num(0), pts(nullptr),
	display_bg_mesh(true), bg_mesh_obj(_gl),
	display_pcls(true), pcls_obj(_gl),
	display_pts(true), pts_obj(_gl),
	display_rigid_circle(true), rc_obj(_gl),
	display_whole_model(true), padding_ratio(0.05f),
	bg_color(0.2f, 0.3f, 0.3f) {}

QtSceneFromModel_T2D_CHM_mt::~QtSceneFromModel_T2D_CHM_mt() {}

void QtSceneFromModel_T2D_CHM_mt::set_viewport(
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

int QtSceneFromModel_T2D_CHM_mt::initialize(int wd, int ht)
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
		Rect bbox = model->get_mesh_bbox();
		if (model->has_rigid_circle())
		{
			Rect rc_bbox;
			model->get_rigid_circle().get_bbox(rc_bbox);
			bbox.envelop(rc_bbox);
		}
		GLfloat xlen = GLfloat(bbox.xu - bbox.xl);
		GLfloat ylen = GLfloat(bbox.yu - bbox.yl);
		GLfloat padding = (xlen > ylen ? xlen : ylen) * padding_ratio;
		display_bbox.xl = GLfloat(bbox.xl) - padding;
		display_bbox.xu = GLfloat(bbox.xu) + padding;
		display_bbox.yl = GLfloat(bbox.yl) - padding;
		display_bbox.yu = GLfloat(bbox.yu) + padding;
	}

	// viewport
	set_viewport(wd, ht, 
		display_bbox.xu - display_bbox.xl,
		display_bbox.yu - display_bbox.yl
		);

	// view matrix
	view_mat.setToIdentity();
	view_mat.ortho(
		display_bbox.xl,
		display_bbox.xu,
		display_bbox.yl,
		display_bbox.yu,
		-1.0f, 1.0f
		);
	shader_plain2D.bind();
	shader_plain2D.setUniformValue("view_mat", view_mat);
	shader_circles.bind();
	shader_circles.setUniformValue("view_mat", view_mat);

	// init bg_mesh
	QVector3D gray(0.5f, 0.5f, 0.5f);
	bg_mesh_obj.init_from_elements(
		model->get_node_pos(),
		model->get_node_num(),
		model->get_elem_node_index(),
		model->get_elem_num(),
		gray);

	// init pcls
	QVector3D moccasin(1.0f, 0.8941f, 0.7098f);
	pcls_obj.init(
		model->get_pcl_pos(),
		model->get_pcl_vol(),
		model->get_pcl_num(),
		moccasin,
		0.5f
		);

	// init rigid circle
	QVector3D light_slate_blue(0.5176f, 0.4392, 1.0f);
	if (model->has_rigid_circle())
	{
		auto &rc = model->get_rigid_circle();
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

void QtSceneFromModel_T2D_CHM_mt::draw()
{
	gl.glViewport(vp_x_pos, vp_y_pos, vp_x_size, vp_y_size);

	gl.glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	gl.glClear(GL_COLOR_BUFFER_BIT);

	shader_plain2D.bind();

	if (display_bg_mesh)
		bg_mesh_obj.draw(shader_plain2D);

	if (model->has_rigid_circle() && display_rigid_circle)
		rc_obj.draw(shader_plain2D);

	shader_circles.bind();

	if (display_pcls)
		pcls_obj.draw(shader_circles);

	if (display_pts && pts && pt_num)
		pts_obj.draw(shader_circles);
}

void QtSceneFromModel_T2D_CHM_mt::resize(int wd, int ht)
{
	set_viewport(wd, ht,
		display_bbox.xu - display_bbox.xl,
		display_bbox.yu - display_bbox.yl
		);
}

int QtSceneFromModel_T2D_CHM_mt::set_pts_from_pcl_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_num = id_num;
	pt_radius = radius;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	const Model_T2D_CHM_mt::Position *pcl_pos = model->get_pcl_pos();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		const Model_T2D_CHM_mt::Position &pp = pcl_pos[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(pp.x);
		pd.y = GLfloat(pp.y);
	}
	return 0;
}

int QtSceneFromModel_T2D_CHM_mt::set_pts_from_node_id(
	size_t* ids, size_t id_num, float radius)
{
	if (!model || !ids || !id_num || radius <= 0.0f)
		return -1;

	pt_num = id_num;
	pt_radius = radius;
	pts_mem.reserve(id_num);
	pts = pts_mem.get_mem();
	const Model_T2D_CHM_mt::Position *node_pos = model->get_node_pos();
	for (size_t p_id = 0; p_id < id_num; ++p_id)
	{
		const Model_T2D_CHM_mt::Position &np = node_pos[ids[p_id]];
		PointData& pd = pts[p_id];
		pd.x = GLfloat(np.x);
		pd.y = GLfloat(np.y);
	}
	return 0;
}

int QtSceneFromModel_T2D_CHM_mt::set_pts_from_vx_s_bc(float radius)
{
	Model_T2D_CHM_mt& md = *model;
	size_t node_num = md.get_node_num();
	const Model_T2D_CHM_mt::NodeHasVBC* nhv = md.get_has_vbc_s();
	const Model_T2D_CHM_mt::Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].has_vx_bc)
		{
			const Model_T2D_CHM_mt::Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T2D_CHM_mt::set_pts_from_vy_s_bc(float radius)
{
	Model_T2D_CHM_mt &md = *model;
	size_t node_num = md.get_node_num();
	const Model_T2D_CHM_mt::NodeHasVBC* nhv = md.get_has_vbc_s();
	const Model_T2D_CHM_mt::Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].has_vy_bc)
		{
			const Model_T2D_CHM_mt::Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T2D_CHM_mt::set_pts_from_vx_f_bc(float radius)
{
	Model_T2D_CHM_mt& md = *model;
	const size_t node_num = md.get_node_num();
	const Model_T2D_CHM_mt::NodeHasVBC* nhv = md.get_has_vbc_f();
	const Model_T2D_CHM_mt::Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].has_vx_bc)
		{
			const Model_T2D_CHM_mt::Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}

int QtSceneFromModel_T2D_CHM_mt::set_pts_from_vy_f_bc(float radius)
{
	Model_T2D_CHM_mt& md = *model;
	const size_t node_num = md.get_node_num();
	const Model_T2D_CHM_mt::NodeHasVBC* nhv = md.get_has_vbc_f();
	const Model_T2D_CHM_mt::Position* node_pos = md.get_node_pos();
	pt_radius = radius;
	pt_num = 0;
	pts_mem.reserve(100);
	PointData pd_tmp;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		if (nhv[n_id].has_vy_bc)
		{
			const Model_T2D_CHM_mt::Position& n = node_pos[n_id];
			pd_tmp.x = n.x;
			pd_tmp.y = n.y;
			pts_mem.add(pd_tmp);
			++pt_num;
		}
	}
	pts = pts_mem.get_mem();
	return 0;
}
