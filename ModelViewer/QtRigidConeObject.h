#ifndef __Qt_Rigid_Cone_Object_h__
#define __Qt_Rigid_Cone_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "Geometry3D.h"

class QtRigidConeObject
{
public:
	struct NodeCoord { GLfloat x, y, z; };

	struct NodeData
	{
		GLuint type;
		GLfloat x, y, z;
		GLfloat nx, ny, nz;
	};

	NodeCoord node_coords[242];
	NodeData node_datas[480 * 3];

protected:
	QOpenGLFunctions_3_3_Core& gl;

	GLuint vao, vbo;
	GLuint vbo_index_num;
	QVector3D color;

	QMatrix4x4 model_mat;
	void form_model_mat(const Point3D &cen);
	void clear();

public:
	QtRigidConeObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtRigidConeObject();

	int init(double x, double y, double z, double h_tip,
			 double h_shaft, double r, QVector3D& c);
	int update(double x, double y, double z);
	void draw(QOpenGLShaderProgram& shader);
};

#endif