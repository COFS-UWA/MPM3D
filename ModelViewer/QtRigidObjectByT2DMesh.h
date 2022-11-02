#ifndef __Qt_Rigid_Object_By_T2D_Mesh_h__
#define __Qt_Rigid_Object_By_T2D_Mesh_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "Geometry2D.h"
#include "RigidObject/RigidObjectByT2DMesh.h"

class QtRigidObjectByT2DMesh
{	
protected:
	struct NodeData
	{
		GLint type;
		GLfloat x, y;
	};

	QOpenGLFunctions_3_3_Core& gl;
	
	GLuint vao, vbo;
	GLsizei vbo_index_num;
	QVector3D color;

	QMatrix4x4 model_mat;
	inline static void form_model_mat(
		const Point2D& cen,
		const Vector2D& ix,
		const Vector2D& iy,
		QMatrix4x4& model_mat
	) noexcept
	{
		float(*mm_data)[4] = reinterpret_cast<float(*)[4]>(model_mat.data());
		// opengl is column major
		// c++ is row major
		mm_data[0][0] = float(ix.x);
		mm_data[0][1] = float(ix.y);
		mm_data[0][2] = 0.0f;
		mm_data[0][3] = 0.0f;
		mm_data[1][0] = float(iy.x);
		mm_data[1][1] = float(iy.y);
		mm_data[1][2] = 0.0f;
		mm_data[1][3] = 0.0f;
		mm_data[2][0] = 0.0f;
		mm_data[2][1] = 0.0f;
		mm_data[2][2] = 1.0f;
		mm_data[2][3] = 0.0f;
		mm_data[3][0] = float(cen.x);
		mm_data[3][1] = float(cen.y);
		mm_data[3][2] = 0.0f;
		mm_data[3][3] = 1.0f;
	}

	void clear();
	
public:
	QtRigidObjectByT2DMesh(QOpenGLFunctions_3_3_Core& _gl);
	~QtRigidObjectByT2DMesh();

	int init(const RigidObjectByT2DMesh& rb, const QVector3D& c);
	int QtRigidObjectByT2DMesh::init_edges(
		const PointToLineDistance* pt_ln_dist, const size_t edge_num,
		double rb_x, double rb_y, double ang, const QVector3D& c);

	void draw(QOpenGLShaderProgram& shader);

	int update(double rb_x, double rb_y, double rb_ang);
};

#endif