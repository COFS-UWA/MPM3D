#include "ModelViewer_pcp.h"

#include <iostream>

#include "QtRigidObjectByT2DMesh.h"

QtRigidObjectByT2DMesh::QtRigidObjectByT2DMesh(
	QOpenGLFunctions_3_3_Core& _gl) : gl(_gl),
	color(1.0f, 1.0f, 1.0f), vbo_index_num(0), vao(0), vbo(0) {}

QtRigidObjectByT2DMesh::~QtRigidObjectByT2DMesh() { clear(); }

void QtRigidObjectByT2DMesh::clear()
{
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

int QtRigidObjectByT2DMesh::init(
	const RigidObjectByT2DMesh& rb,
	const QVector3D& c)
{
	const size_t edge_num = rb.get_edge_num();
	const PointToLineDistance* pt_ln_dist = rb.get_pt_ln_dist();
	if (!pt_ln_dist || !edge_num)
		return -1;

	color = c;
	vbo_index_num = edge_num * 2;

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
	for (size_t e_id = 0; e_id < edge_num; ++e_id)
	{
		const PointToLineDistance &pld = pt_ln_dist[e_id];
		const Point2D& n1 = pld.n1;
		const Point2D& n2 = pld.n2;
		// n1
		cur_nd->type = 0; // monocolor
		cur_nd->x = GLfloat(n1.x);
		cur_nd->y = GLfloat(n1.y);
		++cur_nd;
		// n2
		cur_nd->type = 0; // monocolor
		cur_nd->x = GLfloat(n2.x);
		cur_nd->y = GLfloat(n2.y);
		++cur_nd;
	}
	gl.glBufferData(
		GL_ARRAY_BUFFER,
		vbo_index_num * sizeof(NodeData),
		node_datas,
		GL_STATIC_DRAW);
	delete[] node_datas;

	// v_type
	gl.glVertexAttribIPointer(0,
		1, GL_INT,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, type));
	gl.glEnableVertexAttribArray(0);
	// v_pos
	gl.glVertexAttribPointer(1,
		2, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, x));
	gl.glEnableVertexAttribArray(1);
	// v_color (not used)
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		0, (GLvoid*)0);
	gl.glEnableVertexAttribArray(2);

	form_model_mat(
		rb.get_pos(),
		rb.get_ix(),
		rb.get_iy(),
		model_mat);

	return 0;
}

int QtRigidObjectByT2DMesh::init_edges(
	const PointToLineDistance* pt_ln_dist,
	const size_t edge_num,
	double rb_x,
	double rb_y,
	double ang,
	const QVector3D& c)
{
	if (!pt_ln_dist || !edge_num)
		return -1;

	color = c;
	vbo_index_num = edge_num * 2;

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
	for (size_t e_id = 0; e_id < edge_num; ++e_id)
	{
		const PointToLineDistance& pld = pt_ln_dist[e_id];
		const Point2D& n1 = pld.n1;
		const Point2D& n2 = pld.n2;
		// n1
		cur_nd->type = 0; // monocolor
		cur_nd->x = GLfloat(n1.x);
		cur_nd->y = GLfloat(n1.y);
		++cur_nd;
		// n2
		cur_nd->type = 0; // monocolor
		cur_nd->x = GLfloat(n2.x);
		cur_nd->y = GLfloat(n2.y);
		++cur_nd;
	}
	gl.glBufferData(
		GL_ARRAY_BUFFER,
		vbo_index_num * sizeof(NodeData),
		node_datas,
		GL_STATIC_DRAW);
	delete[] node_datas;

	// v_type
	gl.glVertexAttribIPointer(0,
		1, GL_INT,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, type));
	gl.glEnableVertexAttribArray(0);
	// v_pos
	gl.glVertexAttribPointer(1,
		2, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, x));
	gl.glEnableVertexAttribArray(1);
	// v_color (not used)
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		0, (GLvoid*)0);
	gl.glEnableVertexAttribArray(2);
	
	Point2D pos;
	pos.x = rb_x;
	pos.y = rb_y;
	Vector2D ix, iy;
	rotate_axses_by_angle(ang, ix, iy);
	form_model_mat(pos, ix, iy, model_mat);
	return 0;
}

int QtRigidObjectByT2DMesh::update(double rb_x, double rb_y, double rb_ang)
{
	if (!vbo)
		return 0;
	Point2D pos;
	pos.x = rb_x;
	pos.y = rb_y;
	Vector2D ix, iy;
	rotate_axses_by_angle(rb_ang, ix, iy);
	form_model_mat(pos, ix, iy, model_mat);
	return 0;
}

void QtRigidObjectByT2DMesh::draw(QOpenGLShaderProgram& shader)
{
	if (!vao)
		return;
	shader.bind();
	shader.setUniformValue("g_color", color);
	shader.setUniformValue("model_mat", model_mat);
	gl.glBindVertexArray(vao);
	gl.glDrawArrays(GL_LINES, 0, vbo_index_num);
}
