#include "Simulations_pcp.h"

#include "GLModelWindow.h"

GLModelWindow::GLModelWindow(QWidget *parent) :
	QOpenGLWidget(parent),
	vao_mh(0), vbo_mh(0), ebo_mh(0),
	fov_angle(45.0f), view_dir(0.0f, 1.0f, 0.0f), up_dir(0.0f, 0.0f, 1.0f),
	controller(nullptr)
{

}

GLModelWindow::~GLModelWindow()
{
	if (ebo_mh)
	{
		glDeleteBuffers(1, &ebo_mh);
		ebo_mh = 0;
	}
	if (vbo_mh)
	{
		glDeleteBuffers(1, &vbo_mh);
		vbo_mh = 0;
	}
	if (vao_mh)
	{
		glDeleteVertexArrays(1, &vao_mh);
		vao_mh = 0;
	}
}

void GLModelWindow::set_fov(float _fov) { fov_angle = _fov; }

void GLModelWindow::set_view_dir(float dir_x, float dir_y, float dir_z)
{
	view_dir.setX(dir_x);
	view_dir.setY(dir_y);
	view_dir.setZ(dir_z);
	view_dir.normalize();
}

void GLModelWindow::set_up_dir(float dir_x, float dir_y, float dir_z)
{
	up_dir.setX(dir_x);
	up_dir.setY(dir_y);
	up_dir.setZ(dir_z);
}

void GLModelWindow::initializeGL()
{
	initializeOpenGLFunctions();

	shader.addShaderFromSourceFile(
		QOpenGLShader::Vertex,
		"..\\..\\Asset\\shader_unicolor.vs"
	);
	shader.addShaderFromSourceFile(
		QOpenGLShader::Fragment,
		"..\\..\\Asset\\shader_unicolor.fs"
	);
	shader.link();
	
	glViewport(0, 0, width(), height());
	glEnable(GL_DEPTH_TEST);
	
	if (controller) controller->init_win_data();
}

void GLModelWindow::paintGL()
{
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	update_cam_matrix();
	QVector3D orange(1.0f, 0.647f, 0.310f);
	shader.bind();
	shader.setUniformValue("view_mat", view_mat);
	shader.setUniformValue("proj_mat", proj_mat);
	shader.setUniformValue("color", orange);

	if (vao_mh)
	{
		glLineWidth(1.0f);
		glBindVertexArray(vao_mh);
		glDrawElements(GL_LINES, line_num * 2, GL_UNSIGNED_INT, (GLvoid *)0);
	}

	glFlush();
}

void GLModelWindow::resizeGL(int w, int h)
{
	glViewport(0, 0, w, h);
	//update();
}

void GLModelWindow::update_cam_matrix()
{
	// view matrix
	float dist_from_obj = mh_radius / sin(fov_angle * 0.5 / 180.0 * 3.14159265359);
	QVector3D cam_cen = mh_centre - dist_from_obj * view_dir;
	view_mat.setToIdentity();
	view_mat.lookAt(cam_cen, mh_centre, up_dir);
	// projection matrix
	float aspect_ratio = float(width()) / float(height());
	proj_mat.setToIdentity();
	if (aspect_ratio >= 1.0f)
	{
		proj_mat.perspective(fov_angle, aspect_ratio, 0.01f, 10000.0f);
	}
	else
	{
		float fov_angle2 = atan(tan(fov_angle/180.0*3.14159265359) / aspect_ratio) / 3.14159265359 * 180.0;
		proj_mat.perspective(fov_angle2, aspect_ratio, 0.01f, 10000.0f);
	}
}
