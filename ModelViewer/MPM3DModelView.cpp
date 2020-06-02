#include "ModelViewer_pcp.h"

#include "MPM3DModelView.h"

MPM3DModelView::MPM3DModelView(QWidget *parent) :
	QOpenGLWidget(parent),
	fov_angle(45.0f), view_dist_scale(1.0f), 
	view_dir(1.0f, 0.0f, 0.0f), up_dir(0.0f, 0.0f, 1.0f), 
	bg_color(0.2f, 0.3f, 0.3f),
	// fog effect
	fog_coef(0.0f), fog_color(bg_color),
	// phong model
	light_color(1.0f, 1.0f, 1.0f),
	amb_coef(0.0f), diff_coef(1.0f),
	spec_coef(0.0f), spec_shininess(20.0f),
	light_dir(view_dir), light_dist_scale(2.0f),
	// model data gl buffer
	need_to_paint_bg_mesh(true), bg_mesh_buf(*this),
	need_to_paint_pcl_buf(true), phong_pcl_buf(*this),
	need_to_paint_point_buf(true), point_buf(*this),
	// whether window is fully loaded
	win_is_fully_loaded(false)
{
	connect(&init_timer, SIGNAL(timeout()), this, SLOT(init_timer_func()));
	init_timer.setSingleShot(true);
}

MPM3DModelView::~MPM3DModelView() {}

void MPM3DModelView::initializeGL()
{
	// this function must and can only called here
	initializeOpenGLFunctions();

	glEnable(GL_DEPTH_TEST);
	glViewport(0, 0, width(), height());
	
	shader_unicolor.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"..\\..\\Asset\\shader_unicolor.vert"
		);
	shader_unicolor.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"..\\..\\Asset\\shader_unicolor.frag"
		);
	shader_unicolor.link();
	
	shader_phong.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"..\\..\\Asset\\shader_phong.vert"
	);
	shader_phong.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"..\\..\\Asset\\shader_phong.frag"
	);
	shader_phong.link();
	
	// Need to init md_center and md_radius
	controller->initialize_model_view_data();

	// Need md_center and md_radius cal by initialize_model_view_data().
	update_view_mat();
	update_proj_mat();

	shader_unicolor.bind();
	shader_unicolor.setUniformValue("view_mat", view_mat);
	shader_unicolor.setUniformValue("proj_mat", proj_mat);

	shader_phong.bind();
	shader_phong.setUniformValue("view_mat", view_mat);
	shader_phong.setUniformValue("proj_mat", proj_mat);

	shader_phong.setUniformValue("view_pos", view_pos);

	// fog effect
	shader_phong.setUniformValue("fog_coef", fog_coef);
	shader_phong.setUniformValue("fog_color", fog_color);
	
	// phong model parameters
	update_light_pos();
	shader_phong.setUniformValue("light_pos", light_pos);
	shader_phong.setUniformValue("light_color", light_color);
	shader_phong.setUniformValue("amb_coef", amb_coef);
	shader_phong.setUniformValue("diff_coef", diff_coef);
	shader_phong.setUniformValue("spec_coef", spec_coef);
	shader_phong.setUniformValue("spec_shininess", spec_shininess);
}

void MPM3DModelView::paintGL()
{
	if (!is_fully_loaded())
		return;

	controller->before_render();

	glClearColor(bg_color.x(), bg_color.y(), bg_color.z(), 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	
	// draw bg mesh
	if (need_to_paint_bg_mesh)
	{
		shader_unicolor.bind();
		bg_mesh_buf.draw(shader_unicolor);
	}

	if (need_to_paint_pcl_buf)
	{
		shader_phong.bind();
		phong_pcl_buf.draw();
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

	shader_phong.bind();
	shader_phong.setUniformValue("proj_mat", proj_mat);
}

void MPM3DModelView::update_view_mat()
{
	float dist_from_obj;
	dist_from_obj = md_radius * view_dist_scale
		/ sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	view_dir.normalize();
	view_pos = md_centre - dist_from_obj * view_dir;
	view_mat.setToIdentity();
	view_mat.lookAt(view_pos, md_centre, up_dir);
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

void MPM3DModelView::update_light_pos()
{
	float dist_from_obj;
	dist_from_obj = md_radius * light_dist_scale / sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	light_dir.normalize();
	light_pos = md_centre - dist_from_obj * light_dir;
}

void MPM3DModelView::init_timer_func()
{
	win_is_fully_loaded = true;
	update();
}
