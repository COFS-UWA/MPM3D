#ifndef __Qt_Rigid_Circle_Object_h__
#define __Qt_Rigid_Circle_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

class QtRigidCircleObject
{
public:
	struct NodeData
	{
		GLint type;
		GLfloat x, y;
	};

protected:
	QOpenGLFunctions_3_3_Core& gl;

	GLuint vao, vbo, veo;
	size_t edge_node_num;
	GLfloat line_width;
	QVector3D color;

public:
	QtRigidCircleObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtRigidCircleObject();
	void clear();
	int init(double x, double y, double r,
			 QVector3D& c, GLfloat l_wd);
	int update(double x, double y, double r);
	void draw(QOpenGLShaderProgram& shader);
};

#endif