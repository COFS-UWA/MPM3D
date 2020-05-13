#ifndef __GL_Model_Window_h__
#define __GL_Model_Window_h__

#include <QOpenGLWidget>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "Geometry.h"
#include "ModelToViewerBase.h"

template <typename ModelType> class ModelToViewer3D;

class GLModelWindow : public QOpenGLWidget,
	public QOpenGLFunctions_3_3_Core
{
	Q_OBJECT
	template <typename ModelType> friend class ModelToViewer3D;
protected:
	// shader data
	QOpenGLShaderProgram shader;
	QMatrix4x4 view_mat;
	QMatrix4x4 proj_mat;

	// background mesh data
	GLuint vao_mh, vbo_mh, ebo_mh;
	size_t line_num;

	// camera info
	// fov
	GLfloat fov_angle;
	// camera direction
	QVector3D view_dir;
	QVector3D up_dir;
	// mesh position
	QVector3D mh_centre;
	float mh_radius;

	// deduced view_mat and proj_mat
	void update_cam_matrix();

	ModelToViewerBase *controller;

public:
	explicit GLModelWindow(QWidget *parent = Q_NULLPTR);
	~GLModelWindow();

	void initializeGL();
	void paintGL();
	void resizeGL(int w, int h);
	
	void set_fov(float _fov);
	void set_view_dir(float dir_x, float dir_y, float dir_z);
	void set_up_dir(float dir_x, float dir_y, float dir_z);

	void set_controller(ModelToViewerBase &m2v) { controller = &m2v; }
};

#endif