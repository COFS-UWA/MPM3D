#include "ResultViewerQt_pcp.h"

#include <iostream>

#define _USE_MATH_DEFINES // use PI
#include <math.h>

#include "GLWindow.h"

GLWindow::GLWindow(QWidget *parent) :
	QOpenGLWidget(parent),
	// model
	shader(this),
	vao(0), vbo(0), ebo(0),
	// move camera
	is_roaming(false),
	cam_x(0.0f), cam_y(0.0f), cam_z(0.0f),
	cam_yaw(float(M_PI_2)), cam_pitch(0.0f), cam_roll(0.0f),
	need_update_loc_dir(true),
	loc_x_dir(0.0f, 1.0f, 0.0f), loc_y_dir(-1.0f, 0.0f, 0.0f), loc_z_dir(0.0f, 0.0f, 1.0f),
	ang_sensitivity(1.0f),
	velocity_index(0), velocity0(5.0f), velocity(velocity0),
	// render frame
	should_render_frame(true), min_frame_inv(1000.0f/30.0f), adj_frame_inv(0.0f)
{
	// OpenGL version
	QSurfaceFormat format;
	format.setRenderableType(QSurfaceFormat::OpenGL);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setVersion(3, 3);
	setFormat(format);

	frame_timer.setSingleShot(true);
	connect(&frame_timer, SIGNAL(timeout()), this, SLOT(on_frame_timer_timeout()));
}

GLWindow::~GLWindow()
{
	glDeleteBuffers(1, &vao);
	vao = 0;
	glDeleteBuffers(1, &vbo);
	vbo = 0;
	glDeleteBuffers(1, &ebo);
	ebo = 0;
}

void GLWindow::initializeGL()
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

	// projection matrix
	float aspect_ratio = float(width()) / float(height());
	proj_mat.setToIdentity();
	proj_mat.perspective(45.0f, aspect_ratio, 0.01f, 10000.0f);
	shader.setUniformValue("proj_mat", proj_mat);

	glEnable(GL_DEPTH_TEST);
	
	// Model data
	GLfloat vertices[] =
	{
		-1.0f, -1.0f, -1.0f,
		 1.0f, -1.0f, -1.0f,
		 1.0f,  1.0f, -1.0f,
		-1.0f,  1.0f, -1.0f,
		-1.0f, -1.0f,  1.0f,
		 1.0f, -1.0f,  1.0f,
		 1.0f,  1.0f,  1.0f,
		-1.0f,  1.0f,  1.0f
	};
	
	// edges
	GLuint indices[] =
	{
		0, 1,
		1, 2,
		2, 3,
		3, 0,
		4, 5,
		5, 6,
		6, 7,
		7, 4,
		0, 4,
		1, 5,
		4, 1,
		1, 6,
		2, 6,
		3, 7,
		3, 6,
		4, 3
	};

	glGenVertexArrays(1, &vao);
	glGenBuffers(1, &vbo);
	glGenBuffers(1, &ebo);

	glBindVertexArray(vao); 
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

	// coordinates
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

	cam_time.start();
	frame_time.start();
	last_frame_time = frame_time.elapsed();
}

void GLWindow::paintGL()
{
	if (!should_render_frame)
		return;
	should_render_frame = false;
	float start_rendering_time = frame_time.elapsed();
	
	// view matrix
	shader.bind();
	QVector3D cam_pos(cam_x, cam_y, cam_z);
	update_loc_dir();
	view_mat.setToIdentity();
	view_mat.lookAt(cam_pos, cam_pos + loc_x_dir, loc_z_dir);
	shader.setUniformValue("view_mat", view_mat);

	// color
	QVector3D oriange(1.0f, 1.0f, 1.0f);
	shader.setUniformValue("color", oriange);

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	glLineWidth(1.0f);
	glBindVertexArray(vao);
	glDrawElements(GL_LINES, 32, GL_UNSIGNED_INT, 0);

	glFlush();

	if (frame_timer.isActive())
		return;

	float remaining_time_to_next_frame = min_frame_inv - (frame_time.elapsed() + adj_frame_inv - last_frame_time);
	if (remaining_time_to_next_frame > 0)
	{
		frame_timer.start(remaining_time_to_next_frame);
	}
	else // timeout, render next frame immediately
	{
		should_render_frame = true;
		update();
	}

	frame_time.restart();
	last_frame_time = frame_time.elapsed();
}

void GLWindow::on_frame_timer_timeout()
{
	should_render_frame = true;
	update();
}

void GLWindow::resizeGL(int w, int h)
{
	glViewport(0, 0, w, h);

	// projection matrix
	float aspect_ratio = float(w) / float(h);
	proj_mat.setToIdentity();
	proj_mat.perspective(45.0f, aspect_ratio, 0.01f, 10000.0f);
	shader.bind();
	shader.setUniformValue("proj_mat", proj_mat);
}

void GLWindow::mousePressEvent(QMouseEvent *event)
{
	if (event->button() == Qt::LeftButton)
	{
		is_roaming = true;
		setFocus();
		setMouseTracking(true);

		cursor.setShape(Qt::BlankCursor);
		setCursor(cursor);

		// set curosr to windows centre
		QPoint g_center = mapToGlobal(QPoint(width()/2, height()/2));
		QCursor::setPos(g_center);
	}
	else if (event->button() == Qt::RightButton)
	{
		if (is_roaming)
		{
			is_roaming = false;
			setMouseTracking(false);

			cursor.setShape(Qt::ArrowCursor);
			setCursor(cursor);
		}
	}

	event->accept();
}

void GLWindow::mouseReleaseEvent(QMouseEvent *event)
{
	Q_UNUSED(event);
}

// trim angle into [-pi, pi]
inline void trim_angle_to_pi(float &ang)
{
	if (ang > M_PI)
		ang -= int(( ang + M_PI) / (M_PI*2.0)) * (M_PI*2.0);
	else if (ang < -M_PI)
		ang += int((-ang + M_PI) / (M_PI*2.0)) * (M_PI*2.0);
}

void GLWindow::mouseMoveEvent(QMouseEvent *event)
{
	if (is_roaming)
	{
		float drx =  float(event->x()) / float(width()) - 0.5;
		float dry = -float(event->y()) / float(height()) + 0.5;
		
		QPoint g_center = mapToGlobal(QPoint(width()/2, height()/2));
		QCursor::setPos(g_center);

		cam_yaw -= drx * M_PI * ang_sensitivity;
		trim_angle_to_pi(cam_yaw);

		cam_pitch += dry * M_PI * ang_sensitivity;
		if (cam_pitch > 1.569051f) // 89.9 degree
			cam_pitch = 1.569051f;
		if (cam_pitch < -1.569051f)
			cam_pitch = -1.569051f;

		need_update_loc_dir = true;
	}

	event->accept();
}

void GLWindow::wheelEvent(QWheelEvent *event)
{
	// move camera backward and forward
	update_loc_dir();
	float disp = velocity * 0.002f * float(event->angleDelta().y());
	cam_x += disp * loc_x_dir.x();
	cam_y += disp * loc_x_dir.y();
	cam_z += disp * loc_x_dir.z();

	event->accept();
}

void GLWindow::keyPressEvent(QKeyEvent *event)
{
	const float one_step_time_len = 0.05f;
	float cur_time;
	float disp;

	switch (event->key())
	{
		// velocity adjustment
	case Qt::Key_Minus:
		--velocity_index;
		velocity = velocity0 * exp(float(velocity_index) * 0.05);
		break;
	case Qt::Key_Plus:
		++velocity_index;
		velocity = velocity0 * exp(float(velocity_index) * 0.05);
		break;
		// move the camera
	case Qt::Key_Left:
		update_loc_dir();
		if (event->isAutoRepeat()) // long press
		{
			cur_time = cam_time.elapsed();
			disp = float(cur_time - time_start_moving_left) * 0.001f * velocity;
			time_start_moving_left = cur_time;
		}
		else // one time press
		{
			time_start_moving_left = cam_time.elapsed();
			disp = one_step_time_len * velocity;
		}
		cam_x += disp * loc_y_dir.x();
		cam_y += disp * loc_y_dir.y();
		cam_z += disp * loc_y_dir.z();
		break;
	case Qt::Key_Right:
		update_loc_dir();
		if (event->isAutoRepeat()) // long press
		{
			cur_time = cam_time.elapsed();
			disp = float(cur_time - time_start_moving_left) * 0.001f * velocity;
			time_start_moving_left = cur_time;
		}
		else // one time press
		{
			time_start_moving_left = cam_time.elapsed();
			disp = one_step_time_len * velocity;
		}
		cam_x -= disp * loc_y_dir.x();
		cam_y -= disp * loc_y_dir.y();
		cam_z -= disp * loc_y_dir.z();
		break;
	case Qt::Key_Down:
		update_loc_dir();
		if (event->isAutoRepeat()) // long press
		{
			cur_time = cam_time.elapsed();
			disp = float(cur_time - time_start_moving_left) * 0.001f * velocity;
			time_start_moving_left = cur_time;
		}
		else // one time press
		{
			time_start_moving_left = cam_time.elapsed();
			disp = one_step_time_len * velocity;
		}
		cam_x -= disp * loc_x_dir.x();
		cam_y -= disp * loc_x_dir.y();
		cam_z -= disp * loc_x_dir.z();
		break;
	case Qt::Key_Up:
		update_loc_dir();
		if (event->isAutoRepeat()) // long press
		{
			cur_time = cam_time.elapsed();
			disp = float(cur_time - time_start_moving_left) * 0.001f * velocity;
			time_start_moving_left = cur_time;
		}
		else // one time press
		{
			time_start_moving_left = cam_time.elapsed();
			disp = one_step_time_len * velocity;
		}
		cam_x += disp * loc_x_dir.x();
		cam_y += disp * loc_x_dir.y();
		cam_z += disp * loc_x_dir.z();
		break;
	default:
		break;
	}

	event->accept();
}

void GLWindow::keyReleaseEvent(QKeyEvent *event)
{
	//event->accept();
}
