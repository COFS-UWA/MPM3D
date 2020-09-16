#ifndef __Qt_Rigid_Rect_Object_h__
#define __Qt_Rigid_Rect_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

class QtRigidRectObject
{
public:
	struct NodeData
	{
		GLint type;
		GLfloat x, y;
	};

protected:
	QOpenGLFunctions_3_3_Core& gl;

	double hhx, hhy;
	NodeData node_datas[8];
	GLuint vao, vbo;
	QVector3D color;
	GLfloat line_width;

public:
	QtRigidRectObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtRigidRectObject();
	void clear();
	int init(double x, double y, double ang,
			 double hx, double hy,
			 QVector3D& c, GLfloat line_wd);
	int update(double x, double y, double ang);
	void draw(QOpenGLShaderProgram& shader);
};

#endif