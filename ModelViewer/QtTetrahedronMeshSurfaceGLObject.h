#ifndef __Qt_Tetrahedron_Mesh_Surface_GL_Object_h__
#define __Qt_Tetrahedron_Mesh_Surface_GL_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "Geometry3D.h"

class QtTetrahedronMeshSurfaceGLObject
{
protected:
	QOpenGLFunctions_3_3_Core &gl;

	struct NodeData
	{
		GLuint type;
		GLfloat x, y, z;
		GLfloat nx, ny, nz;
	};

	GLuint vao, vbo;
	size_t face_node_num;
	QVector3D color;

	void clear();
	
public:
	QtTetrahedronMeshSurfaceGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtTetrahedronMeshSurfaceGLObject();

	// Node has members x, y, z
	// Edge has members n1, n2, n3
	template <typename Node, typename Face>
	int init_from_faces(Node* nodes, size_t node_num, Face* faces, size_t face_num, QVector3D& c);
	
	void draw(QOpenGLShaderProgram& shader);
};

template <typename Node, typename Face>
int QtTetrahedronMeshSurfaceGLObject::init_from_faces(
	Node* nodes, size_t node_num,
	Face* faces, size_t face_num,
	QVector3D& c
	)
{
	if (!nodes || !node_num ||
		!faces || !face_num)
		return -1;

	color = c;
	face_node_num = face_num * 3;

	gl.glGenVertexArrays(1, &vao);
	if (vao == 0)
		return -2;
	gl.glBindVertexArray(vao);

	gl.glGenBuffers(1, &vbo);
	if (vbo == 0)
		return -2;
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);
	NodeData* node_datas = new NodeData[face_node_num];
	NodeData* cur_nd = node_datas;
	QVector3D p21, p31, normal;
	for (size_t f_id = 0; f_id < face_num; ++f_id)
	{
		Face &f = faces[f_id];
		Node& n1 = nodes[f.n1];
		Node& n2 = nodes[f.n2];
		Node& n3 = nodes[f.n3];
		// cal normal
		p21.setX(GLfloat(n2.x - n1.x));
		p21.setY(GLfloat(n2.y - n1.y));
		p21.setZ(GLfloat(n2.z - n1.z));
		p31.setX(GLfloat(n3.x - n1.x));
		p31.setY(GLfloat(n3.y - n1.y));
		p31.setZ(GLfloat(n3.z - n1.z));
		normal = QVector3D::crossProduct(p31, p21);
		normal.normalize();
		// swap n2 and n3 to be counter-clockwise
		// for external surface
		// n1
		cur_nd->type = 0;
		cur_nd->x = GLfloat(n1.x);
		cur_nd->y = GLfloat(n1.y);
		cur_nd->z = GLfloat(n1.z);
		cur_nd->nx = normal.x();
		cur_nd->ny = normal.y();
		cur_nd->nz = normal.z();
		++cur_nd;
		// n3
		cur_nd->type = 0;
		cur_nd->x = GLfloat(n3.x);
		cur_nd->y = GLfloat(n3.y);
		cur_nd->z = GLfloat(n3.z);
		cur_nd->nx = normal.x();
		cur_nd->ny = normal.y();
		cur_nd->nz = normal.z();
		++cur_nd;
		// n2
		cur_nd->type = 0;
		cur_nd->x = GLfloat(n2.x);
		cur_nd->y = GLfloat(n2.y);
		cur_nd->z = GLfloat(n2.z);
		cur_nd->nx = normal.x();
		cur_nd->ny = normal.y();
		cur_nd->nz = normal.z();
		++cur_nd;
	}
	gl.glBufferData(
		GL_ARRAY_BUFFER,
		face_node_num * sizeof(NodeData),
		node_datas,
		GL_STATIC_DRAW
		);
	delete[] node_datas;

	// v_type
	gl.glVertexAttribIPointer(0,
		1, GL_UNSIGNED_INT,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, type)
		);
	gl.glEnableVertexAttribArray(0);
	// v_pos
	gl.glVertexAttribPointer(1,
		3, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid *)offsetof(NodeData, x)
		);
	gl.glEnableVertexAttribArray(1);
	// v_normal
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid *)offsetof(NodeData, nx)
		);
	gl.glEnableVertexAttribArray(2);

	return 0;
}

#endif