#ifndef __Qt_Rigid_Tetrahedron_Mesh_GL_Object_h__
#define __Qt_Rigid_Tetrahedron_Mesh_GL_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "Geometry3D.h"
#include "RigidBody/RigidTetrahedronMesh.h"

class QtRigidTetrahedronMeshGLObject
{
public:
	enum DisplayMode : unsigned char
	{
		Invalid = 0,
		Surface = 1,
		LineFrame = 2
	};

protected:
	struct FaceNodeData
	{
		GLuint type;
		GLfloat x, y, z;
		GLfloat nx, ny, nz;
	};

	struct LineNodeData
	{
		GLuint type;
		GLfloat x, y, z;
	};

	typedef RigidTetrahedronMesh::Node RbNode;
	typedef RigidTetrahedronMesh::Face RbFace;

	QOpenGLFunctions_3_3_Core& gl;

	DisplayMode mode;
	QVector3D color;
	
	GLuint vao, vbo, veo;
	GLsizei vbo_index_num;

	QMatrix4x4 model_mat;
	void form_model_mat(const Point3D &cen, const Vector3D &ix,
						const Vector3D &iy, const Vector3D &iz);

	void clear();
	
public:
	QtRigidTetrahedronMeshGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtRigidTetrahedronMeshGLObject();

	inline DisplayMode get_display_mode() const noexcept { return mode; }

	// pre process
	int init_faces(const RigidTetrahedronMesh& rb, const QVector3D& c);
	int init_line_frame(const RigidTetrahedronMesh &rb, const QVector3D &c);
	int update(const RigidTetrahedronMesh& rb);

	// post process
	template <typename Node, typename Element>
	int init_faces(const Node *nodes, const size_t node_num,
				   const Element *elems, const size_t elem_num,
				   const Point3D& cen, const Vector3D& ang, const QVector3D& c);
	template <typename Node, typename Element>
	int init_line_frame(const Node* nodes, const size_t node_num,
						const Element* elems, const size_t elem_num,
						const Point3D& cen, const Vector3D& ang, const QVector3D& c);
	int update(const Point3D &cen, const Vector3D &ang);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename Node, typename Element>
int QtRigidTetrahedronMeshGLObject::init_faces(
	const Node* nodes,
	const size_t node_num,
	const Element* elems,
	const size_t elem_num,
	const Point3D& cen,
	const Vector3D& ang,
	const QVector3D& c
	)
{
	if (!nodes || !node_num || !elems || !elem_num)
		return -1;
	
	struct Face { size_t n1, n2, n3; };
	size_t bface_num;
	const Face *bfaces = extract_edges_from_triangles<Face>(elems, elem_num, bface_num);
	if (!bfaces || !bface_num)
		return -1;

	mode = Surface;
	color = c;
	vbo_index_num = bface_num * 3;

	gl.glGenVertexArrays(1, &vao);
	if (vao == 0)
		return -2;
	gl.glBindVertexArray(vao);

	gl.glGenBuffers(1, &vbo);
	if (vbo == 0)
		return -2;
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);

	FaceNodeData* node_datas = new FaceNodeData[vbo_index_num];
	FaceNodeData* cur_nd = node_datas;
	Vector3D face_edge21, face_edge31, face_normal;
	for (size_t f_id = 0; f_id < bface_num; ++f_id)
	{
		const Face& f = bfaces[f_id];
		const Node& n1 = nodes[f.n1];
		const Node& n2 = nodes[f.n2];
		const Node& n3 = nodes[f.n3];
		// cal normal
		face_edge21.substract(n2.x, n2.y, n2.z, n1.x, n1.y, n1.z);
		face_edge31.substract(n3.x, n3.y, n3.z, n1.x, n1.y, n1.z);
		face_normal.cross(face_edge31, face_edge21);
		face_normal.normalize();
		// swap n2 and n3 to be counter-clockwise
		// for external surface
		// n1
		cur_nd->type = 0;
		cur_nd->x = GLfloat(n1.x);
		cur_nd->y = GLfloat(n1.y);
		cur_nd->z = GLfloat(n1.z);
		cur_nd->nx = GLfloat(face_normal.x);
		cur_nd->ny = GLfloat(face_normal.y);
		cur_nd->nz = GLfloat(face_normal.z);
		++cur_nd;
		// n3
		cur_nd->type = 0;
		cur_nd->x = GLfloat(n3.x);
		cur_nd->y = GLfloat(n3.y);
		cur_nd->z = GLfloat(n3.z);
		cur_nd->nx = GLfloat(face_normal.x);
		cur_nd->ny = GLfloat(face_normal.y);
		cur_nd->nz = GLfloat(face_normal.z);
		++cur_nd;
		// n2
		cur_nd->type = 0;
		cur_nd->x = GLfloat(n2.x);
		cur_nd->y = GLfloat(n2.y);
		cur_nd->z = GLfloat(n2.z);
		cur_nd->nx = GLfloat(face_normal.x);
		cur_nd->ny = GLfloat(face_normal.y);
		cur_nd->nz = GLfloat(face_normal.z);
		++cur_nd;
	}
	gl.glBufferData(
		GL_ARRAY_BUFFER,
		node_num * sizeof(FaceNodeData),
		node_datas,
		GL_STATIC_DRAW
		);
	delete[] node_datas;
	delete[] bfaces;

	// v_type
	gl.glVertexAttribIPointer(0,
		1, GL_UNSIGNED_INT,
		sizeof(FaceNodeData),
		(GLvoid*)offsetof(FaceNodeData, type)
		);
	gl.glEnableVertexAttribArray(0);
	// v_pos
	gl.glVertexAttribPointer(1,
		3, GL_FLOAT, GL_FALSE,
		sizeof(FaceNodeData),
		(GLvoid*)offsetof(FaceNodeData, x)
		);
	gl.glEnableVertexAttribArray(1);
	// v_normal
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		sizeof(FaceNodeData),
		(GLvoid*)offsetof(FaceNodeData, nx)
		);
	gl.glEnableVertexAttribArray(2);

	Vector3D ix(1.0, 0.0, 0.0);
	Vector3D iy(0.0, 1.0, 0.0);
	Vector3D iz(0.0, 0.0, 1.0);
	rotate_axses_by_angle(ang, ix, iy, iz);
	form_model_mat(cen, ix, iy, iz);

	return 0;
}

template <typename Node, typename Element>
int QtRigidTetrahedronMeshGLObject::init_line_frame(
	const Node* nodes,
	const size_t node_num,
	const Element* elems,
	const size_t elem_num,
	const Point3D& cen,
	const Vector3D& ang,
	const QVector3D& c
	)
{
	int res = init_face<Node, Element>(nodes, node_num, elems, elem_num, cen, ang, c);
	if (res)
		return res;
	mode = LineFrame;
	return 0;
}

#endif