#ifndef __Qt_Rigid_Cube_Object_h__
#define __Qt_Rigid_Cube_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "Geometry3D.h"

class QtRigidCubeObject
{
public:
	struct NodeCoord { GLfloat x, y, z; };

	struct NodeData
	{
		GLuint type;
		GLfloat x, y, z;
		GLfloat nx, ny, nz;
	};

	NodeCoord node_coords[8];
	NodeData node_datas[12 * 3];

protected:
	QOpenGLFunctions_3_3_Core& gl;

	GLuint vao, vbo;
	GLuint vbo_index_num;
	QVector3D color;

	QMatrix4x4 model_mat;
	void form_model_mat(const Point3D& cen);
	void clear();

public:
	QtRigidCubeObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtRigidCubeObject();

	int init(double _x, double _y, double _z,
			 double _hx, double _hy, double _hz,
			 QVector3D& c);
	int update(double x, double y, double z);
	void draw(QOpenGLShaderProgram& shader);
};

#endif