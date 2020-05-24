#ifndef __Point_Buffer_h__
#define __Point_Buffer_h__

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

#include "Geometry.h"

class PointBuffer
{
protected:
	QOpenGLFunctions_3_3_Core &glwp;
	GLuint vao, vbo;
	GLsizei point_num;
	GLfloat point_size;
	QVector3D color;

public:
	PointBuffer(QOpenGLFunctions_3_3_Core &wp) : glwp(wp),
		vao(0), vbo(0), point_num(0), color(1.0f, 1.0f, 1.0f) {}
	~PointBuffer() { clear(); }
	void clear();

	inline void draw(QOpenGLShaderProgram& shader)
	{
		if (vao)
		{
			glwp.glBindVertexArray(vao);
			shader.setUniformValue("color", color);
			glwp.glPointSize(point_size);
			glwp.glDrawArrays(GL_POINTS, 0, point_num);
		}
	}

	int init(Point3D* pts, size_t pt_num, GLfloat pt_size, QVector3D& _color);
};

#endif