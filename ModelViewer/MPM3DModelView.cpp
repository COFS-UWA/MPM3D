#include "ModelViewer_pcp.h"

#include "MPM3DModelView.h"

MPM3DModelView::MPM3DModelView(QWidget *parent) :
	QOpenGLWidget(parent),
	fov_angle(45.0f), view_dist_scale(1.0f), 
	view_dir(1.0f, 0.0f, 0.0f), up_dir(0.0f, 0.0f, 1.0f), 
	bg_mesh_buf(*this), need_to_paint_bg_mesh(true),
	pcl_buf(*this), cpcl_buf(*this),
	is_monocolor_pcl(true), need_to_paint_pcl_buf(true),
	point_buf(*this), need_to_paint_point_buf(true),
	win_is_fully_loaded(false)
{
	connect(&init_timer, SIGNAL(timeout()), this, SLOT(init_timer_func()));
	init_timer.setSingleShot(true);
}

MPM3DModelView::~MPM3DModelView() { win_is_fully_loaded = false; }

void MPM3DModelView::initializeGL()
{
	// this function must and can only called here
	initializeOpenGLFunctions();

	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, width(), height());

	shader_unicolor.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"..\\..\\Asset\\shader_unicolor.vs"
		);
	shader_unicolor.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"..\\..\\Asset\\shader_unicolor.fs"
		);
	shader_unicolor.link();
	
	shader_multicolor.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"..\\..\\Asset\\shader_multicolor.vs"
		);
	shader_multicolor.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"..\\..\\Asset\\shader_multicolor.fs"
		);
	shader_multicolor.link();

	// need to init bounding circle
	controller->initialize_model_view_data();

	update_view_mat();
	update_proj_mat();

	shader_unicolor.bind();
	shader_unicolor.setUniformValue("view_mat", view_mat);
	shader_unicolor.setUniformValue("proj_mat", proj_mat);

	shader_multicolor.bind();
	shader_multicolor.setUniformValue("view_mat", view_mat);
	shader_multicolor.setUniformValue("proj_mat", proj_mat);
}

void MPM3DModelView::paintGL()
{
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	
	if (!is_fully_loaded())
		return;

	controller->before_render();

	shader_unicolor.bind();

	if (need_to_paint_bg_mesh)
		bg_mesh_buf.draw(shader_unicolor);

	if (need_to_paint_pcl_buf)
	{
		if (is_monocolor_pcl)
		{
			pcl_buf.draw(shader_unicolor);
		}
		else
		{
			shader_multicolor.bind();
			cpcl_buf.draw();
		}
	}
	
	if (need_to_paint_point_buf)
	{
		shader_unicolor.bind();
		point_buf.draw(shader_unicolor);
	}

	controller->after_render();
}

void MPM3DModelView::resizeGL(int width, int height)
{
	static size_t call_resize_time = 0;
	// The window is fully loaded after resized twice
	if (call_resize_time == 1)
		init_timer.start(100);
	++call_resize_time;

	glViewport(0, 0, width, height);

	update_proj_mat();
	shader_unicolor.bind();
	shader_unicolor.setUniformValue("proj_mat", proj_mat);
	shader_multicolor.bind();
	shader_multicolor.setUniformValue("proj_mat", proj_mat);
}

void MPM3DModelView::update_view_mat()
{
	float dist_from_obj;
	dist_from_obj = md_radius * view_dist_scale / sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	view_dir.normalize();
	QVector3D cam_cen = md_centre - dist_from_obj * view_dir;
	view_mat.setToIdentity();
	view_mat.lookAt(cam_cen, md_centre, up_dir);
}

void MPM3DModelView::update_proj_mat()
{
	float aspect_ratio = float(width()) / float(height());
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

void MPM3DModelView::init_timer_func()
{
	win_is_fully_loaded = true;
	update();
}
