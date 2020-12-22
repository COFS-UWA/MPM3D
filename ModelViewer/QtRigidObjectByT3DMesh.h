#ifndef __Qt_Rigid_Object_By_T3D_Mesh_h__
#define __Qt_Rigid_Object_By_T3D_Mesh_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "Geometry3D.h"
#include "RigidObject/RigidObjectByT3DMesh.h"

class QtRigidObjectByT3DMesh
{
public:
	enum DisplayMode : unsigned char
	{
		Invalid = 0,
		Surface = 1,
		LineFrame = 2
	};
	
protected:
	struct NodeData
	{
		GLuint type;
		GLfloat x, y, z;
		GLfloat nx, ny, nz;
	};

	QOpenGLFunctions_3_3_Core& gl;

	DisplayMode mode;
	QVector3D color;

	GLuint vao, vbo;
	GLsizei vbo_index_num;

	int init_face_data(const RigidObjectByT3DMesh& rb, const QVector3D& c);

	QMatrix4x4 model_mat;
	inline static void form_model_mat(
		const Point3D& cen,
		const Vector3D& ix,
		const Vector3D& iy,
		const Vector3D& iz,
		QMatrix4x4 &model_mat) noexcept
	{
		float(*mm_data)[4] = reinterpret_cast<float(*)[4]>(model_mat.data());
		// opengl is column major
		// c++ is row major
		mm_data[0][0] = float(ix.x);
		mm_data[0][1] = float(ix.y);
		mm_data[0][2] = float(ix.z);
		mm_data[0][3] = 0.0f;
		mm_data[1][0] = float(iy.x);
		mm_data[1][1] = float(iy.y);
		mm_data[1][2] = float(iy.z);
		mm_data[1][3] = 0.0f;
		mm_data[2][0] = float(iz.x);
		mm_data[2][1] = float(iz.y);
		mm_data[2][2] = float(iz.z);
		mm_data[2][3] = 0.0f;
		mm_data[3][0] = float(cen.x);
		mm_data[3][1] = float(cen.y);
		mm_data[3][2] = float(cen.z);
		mm_data[3][3] = 1.0f;		
	}

	void clear();
	
public:
	QtRigidObjectByT3DMesh(QOpenGLFunctions_3_3_Core& _gl);
	~QtRigidObjectByT3DMesh();

	inline DisplayMode get_display_mode() const noexcept { return mode; }

	// pre process
	inline int init_faces(const RigidObjectByT3DMesh &rb, const QVector3D &c)
	{
		mode = DisplayMode::Surface;
		return init_face_data(rb, c);
	}
	inline int init_line_frame(const RigidObjectByT3DMesh &rb, const QVector3D &c)
	{
		mode = DisplayMode::LineFrame;
		return init_line_frame(rb, c);
	}
	int update(const RigidObjectByT3DMesh &rb);

	//// post process
	//template <typename Node, typename Element>
	//int init_faces(const Node* nodes, const size_t node_num,
	//	const Element* elems, const size_t elem_num,
	//	const Point3D& cen, const Vector3D& ang, const QVector3D& c);
	//template <typename Node, typename Element>
	//int init_line_frame(const Node* nodes, const size_t node_num,
	//	const Element* elems, const size_t elem_num,
	//	const Point3D& cen, const Vector3D& ang, const QVector3D& c);
	//int update(const Point3D& cen, const Vector3D& ang);

	void draw(QOpenGLShaderProgram& shader);
};

#endif