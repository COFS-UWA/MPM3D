#include "ModelViewer_pcp.h"

#include <iostream>

#include "QtRigidObjectByT3DMesh.h"

QtRigidObjectByT3DMesh::QtRigidObjectByT3DMesh(
	QOpenGLFunctions_3_3_Core& _gl) : gl(_gl),
	color(1.0f, 1.0f, 1.0f), vao(0), vbo(0),
	vbo_index_num(0), mode(DisplayMode::Invalid) {}

QtRigidObjectByT3DMesh::~QtRigidObjectByT3DMesh() { clear(); }

void QtRigidObjectByT3DMesh::clear()
{
	mode = DisplayMode::Invalid;
	if (vbo)
	{
		gl.glDeleteBuffers(1, &vbo);
		vbo = 0;
	}
	if (vao)
	{
		gl.glDeleteVertexArrays(1, &vao);
		vao = 0;
	}
	vbo_index_num = 0;
}

int QtRigidObjectByT3DMesh::init_face_data(
	const RigidObjectByT3DMesh& rb, const QVector3D& c)
{
	const size_t face_num = rb.get_face_num();
	const PointToTriangleDistance* pt_tri_dist = rb.get_pt_tri_dist();
	if (!pt_tri_dist || !face_num)
		return -1;

	color = c;
	vbo_index_num = face_num * 3;

	gl.glGenVertexArrays(1, &vao);
	if (vao == 0)
		return -2;
	gl.glBindVertexArray(vao);

	gl.glGenBuffers(1, &vbo);
	if (vbo == 0)
		return -2;
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);

	NodeData* node_datas = new NodeData[vbo_index_num];
	NodeData* cur_nd = node_datas;
	Vector3D face_edge21, face_edge31, face_normal;
	for (size_t f_id = 0; f_id < face_num; ++f_id)
	{
		const PointToTriangleDistance& ptd = pt_tri_dist[f_id];
		const Point3D& n1 = ptd.n1;
		const Point3D& n2 = ptd.n2;
		const Point3D& n3 = ptd.n3;
		// cal normal
		face_edge21.substract(n2.x, n2.y, n2.z, n1.x, n1.y, n1.z);
		face_edge31.substract(n3.x, n3.y, n3.z, n1.x, n1.y, n1.z);
		face_normal.cross(face_edge31, face_edge21);
		face_normal.normalize();
		// swap n2 and n3 to be counter-clockwise
		// for external surface
		// n1
		cur_nd->type = 1;
		cur_nd->x = GLfloat(n1.x);
		cur_nd->y = GLfloat(n1.y);
		cur_nd->z = GLfloat(n1.z);
		cur_nd->nx = GLfloat(face_normal.x);
		cur_nd->ny = GLfloat(face_normal.y);
		cur_nd->nz = GLfloat(face_normal.z);
		++cur_nd;
		// n3
		cur_nd->type = 1;
		cur_nd->x = GLfloat(n3.x);
		cur_nd->y = GLfloat(n3.y);
		cur_nd->z = GLfloat(n3.z);
		cur_nd->nx = GLfloat(face_normal.x);
		cur_nd->ny = GLfloat(face_normal.y);
		cur_nd->nz = GLfloat(face_normal.z);
		++cur_nd;
		// n2
		cur_nd->type = 1;
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
		vbo_index_num * sizeof(NodeData),
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
		(GLvoid*)offsetof(NodeData, x)
		);
	gl.glEnableVertexAttribArray(1);
	// v_normal
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, nx)
		);
	gl.glEnableVertexAttribArray(2);

	form_model_mat(
		rb.get_pos(),
		rb.get_ix(),
		rb.get_iy(),
		rb.get_iz(),
		model_mat);

	return 0;
}

int QtRigidObjectByT3DMesh::update(const RigidObjectByT3DMesh& rb)
{
	if (mode == DisplayMode::Invalid || !vbo)
		return 0;
	form_model_mat(
		rb.get_pos(),
		rb.get_ix(),
		rb.get_iy(),
		rb.get_iz(),
		model_mat);
	return 0;
}

void QtRigidObjectByT3DMesh::draw(QOpenGLShaderProgram& shader)
{
	if (mode == DisplayMode::Invalid || !vao)
		return;

	shader.bind();
	shader.setUniformValue("g_color", color);
	shader.setUniformValue("model_mat", model_mat);
	gl.glBindVertexArray(vao);

	// need to adjust color of tetrahedron mesh
	if (mode == DisplayMode::Surface)
	{
		//gl.glFrontFace(GL_CW);
		gl.glCullFace(GL_BACK);
		gl.glDrawArrays(GL_TRIANGLES, 0, vbo_index_num);
	}
	else if (mode == DisplayMode::LineFrame)
	{
		gl.glPolygonMode(GL_FRONT, GL_LINE);
		gl.glDrawArrays(GL_TRIANGLES, 0, vbo_index_num);
		gl.glPolygonMode(GL_FRONT, GL_FILL);
	}
}

int QtRigidObjectByT3DMesh::init_faces(
	const PointToTriangleDistance* pt_tri_dist,
	const size_t face_num,
	const Point3D& cen,
	const Vector3D& ang,
	const QVector3D& c)
{
	if (!pt_tri_dist || !face_num)
		return -1;

	mode = DisplayMode::Surface;
	color = c;
	vbo_index_num = face_num * 3;

	gl.glGenVertexArrays(1, &vao);
	if (vao == 0)
		return -2;
	gl.glBindVertexArray(vao);

	gl.glGenBuffers(1, &vbo);
	if (vbo == 0)
		return -2;
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);

	NodeData* node_datas = new NodeData[vbo_index_num];
	NodeData* cur_nd = node_datas;
	Vector3D face_edge21, face_edge31, face_normal;
	for (size_t f_id = 0; f_id < face_num; ++f_id)
	{
		const PointToTriangleDistance& ptd = pt_tri_dist[f_id];
		const Point3D& n1 = ptd.n1;
		const Point3D& n2 = ptd.n2;
		const Point3D& n3 = ptd.n3;
		// cal normal
		face_edge21.substract(n2.x, n2.y, n2.z, n1.x, n1.y, n1.z);
		face_edge31.substract(n3.x, n3.y, n3.z, n1.x, n1.y, n1.z);
		face_normal.cross(face_edge31, face_edge21);
		face_normal.normalize();
		// swap n2 and n3 to be counter-clockwise
		// for external surface
		// n1
		cur_nd->type = 1;
		cur_nd->x = GLfloat(n1.x);
		cur_nd->y = GLfloat(n1.y);
		cur_nd->z = GLfloat(n1.z);
		cur_nd->nx = GLfloat(face_normal.x);
		cur_nd->ny = GLfloat(face_normal.y);
		cur_nd->nz = GLfloat(face_normal.z);
		++cur_nd;
		// n3
		cur_nd->type = 1;
		cur_nd->x = GLfloat(n3.x);
		cur_nd->y = GLfloat(n3.y);
		cur_nd->z = GLfloat(n3.z);
		cur_nd->nx = GLfloat(face_normal.x);
		cur_nd->ny = GLfloat(face_normal.y);
		cur_nd->nz = GLfloat(face_normal.z);
		++cur_nd;
		// n2
		cur_nd->type = 1;
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
		vbo_index_num * sizeof(NodeData),
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
		(GLvoid*)offsetof(NodeData, x)
	);
	gl.glEnableVertexAttribArray(1);
	// v_normal
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, nx)
	);
	gl.glEnableVertexAttribArray(2);

	Vector3D ix, iy, iz;
	rotate_axses_by_angle(ang, ix, iy, iz);
	form_model_mat(cen, ix, iy, iz, model_mat);
	return 0;
}

int QtRigidObjectByT3DMesh::update(
	const Point3D& cen,
	const Vector3D& ang)
{
	if (mode == DisplayMode::Invalid || !vbo)
		return 0;
	Vector3D ix, iy, iz;
	rotate_axses_by_angle(ang, ix, iy, iz);
	form_model_mat(cen, ix, iy, iz, model_mat);
	return 0;
}
