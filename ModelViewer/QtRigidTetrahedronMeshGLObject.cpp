#include "ModelViewer_pcp.h"

#include "QtRigidTetrahedronMeshGLObject.h"

QtRigidTetrahedronMeshGLObject::QtRigidTetrahedronMeshGLObject(QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), mode(Invalid), color(1.0f, 1.0f, 1.0f),
	vao(0), vbo(0), veo(0), vbo_index_num(0) {}

QtRigidTetrahedronMeshGLObject::~QtRigidTetrahedronMeshGLObject() { clear(); }

void QtRigidTetrahedronMeshGLObject::clear()
{
	if (veo)
	{
		gl.glDeleteBuffers(1, &veo);
		veo = 0;
	}
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
	mode = Invalid;
}

int QtRigidTetrahedronMeshGLObject::init_faces(
	const RigidTetrahedronMesh &rb,
	const QVector3D &c
	)
{
	size_t node_num = rb.get_node_num();
	const RbNode* nodes = rb.get_nodes();
	size_t bface_num = rb.get_bface_num();
	const RbFace* bfaces = rb.get_bfaces();
	if (!nodes  || !node_num ||
		!bfaces || !bface_num)
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

	FaceNodeData *node_datas = new FaceNodeData[vbo_index_num];
	FaceNodeData *cur_nd = node_datas;
	Vector3D face_edge21, face_edge31, face_normal;
	for (size_t f_id = 0; f_id < bface_num; ++f_id)
	{
		const RbFace& f = bfaces[f_id];
		const RbNode& n1 = nodes[f.n1];
		const RbNode& n2 = nodes[f.n2];
		const RbNode& n3 = nodes[f.n3];
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
		vbo_index_num * sizeof(FaceNodeData),
		node_datas,
		GL_STATIC_DRAW
		);
	delete[] node_datas;

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
	
	const Point3D& cen = rb.get_centre();
	const Vector3D& ix = rb.get_ix();
	const Vector3D& iy = rb.get_iy();
	const Vector3D& iz = rb.get_iz();
	form_model_mat(cen, ix, iy, iz);
	
	return 0;
}

int QtRigidTetrahedronMeshGLObject::init_line_frame(
	const RigidTetrahedronMesh& rb,
	const QVector3D& c
	)
{
	size_t node_num = rb.get_node_num();
	const RbNode* nodes = rb.get_nodes();
	size_t bface_num = rb.get_bface_num();
	const RbFace* bfaces = rb.get_bfaces();
	if (!nodes || !node_num ||
		!bfaces || !bface_num)
		return -1;

	struct Line { size_t n1, n2; };
	size_t line_num;
	Line* line_datas = extract_edges_from_triangles<Line, RbFace>(bfaces, bface_num, line_num);
	if (!line_datas || !line_num)
		return 0;

	mode = LineFrame;
	color = c;
	vbo_index_num = line_num * 2;

	gl.glGenVertexArrays(1, &vao);
	if (vao == 0)
		return -2;
	gl.glBindVertexArray(vao);

	gl.glGenBuffers(1, &vbo);
	if (vbo == 0)
		return -2;
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);

	LineNodeData* node_datas = new LineNodeData[vbo_index_num];
	LineNodeData* cur_nd = node_datas;
	for (size_t l_id = 0; l_id < line_num; ++l_id)
	{
		const Line& l = line_datas[l_id];
		const RbNode& n1 = nodes[l.n1];
		const RbNode& n2 = nodes[l.n2];
		// n1
		cur_nd->type = 2;
		cur_nd->x = GLfloat(n1.x);
		cur_nd->y = GLfloat(n1.y);
		cur_nd->z = GLfloat(n1.z);
		++cur_nd;
		// n2
		cur_nd->type = 2;
		cur_nd->x = GLfloat(n2.x);
		cur_nd->y = GLfloat(n2.y);
		cur_nd->z = GLfloat(n2.z);
		++cur_nd;
	}
	gl.glBufferData(
		GL_ARRAY_BUFFER,
		vbo_index_num * sizeof(LineNodeData),
		node_datas,
		GL_STATIC_DRAW
		);
	delete[] node_datas;
	delete[] line_datas;

	// v_type
	gl.glVertexAttribIPointer(0,
		1, GL_UNSIGNED_INT,
		sizeof(LineNodeData),
		(GLvoid*)offsetof(LineNodeData, type)
		);
	gl.glEnableVertexAttribArray(0);
	// v_pos
	gl.glVertexAttribPointer(1,
		3, GL_FLOAT, GL_FALSE,
		sizeof(LineNodeData),
		(GLvoid*)offsetof(LineNodeData, x)
		);
	gl.glEnableVertexAttribArray(1);

	const Point3D& cen = rb.get_centre();
	const Vector3D& ix = rb.get_ix();
	const Vector3D& iy = rb.get_iy();
	const Vector3D& iz = rb.get_iz();
	form_model_mat(cen, ix, iy, iz);

	return 0;
}

int QtRigidTetrahedronMeshGLObject::update(const RigidTetrahedronMesh& rb)
{
	if (mode == Invalid || !vbo)
		return 0;

	const Point3D& cen = rb.get_centre();
	const Vector3D &ix = rb.get_ix();
	const Vector3D &iy = rb.get_iy();
	const Vector3D &iz = rb.get_iz();
	form_model_mat(cen, ix, iy, iz);
	return 0;
}

int QtRigidTetrahedronMeshGLObject::update(
	const Point3D& cen,
	const Vector3D& ang
	)
{
	if (mode == Invalid || !vbo)
		return 0;

	Vector3D ix(1.0, 0.0, 0.0);
	Vector3D iy(0.0, 1.0, 0.0);
	Vector3D iz(0.0, 0.0, 1.0);
	rotate_axses_by_angle(ang, ix, iy, iz);
	form_model_mat(cen, ix, iy, iz);
	return 0;
}

void QtRigidTetrahedronMeshGLObject::draw(QOpenGLShaderProgram& shader)
{
	if (mode == Invalid || !vao)
		return;

	shader.setUniformValue("model_mat", model_mat);
	shader.setUniformValue("g_color", color);
	gl.glBindVertexArray(vao);

	// need to adjust color of tetrahedron mesh
	if (mode == Surface)
	{
		//gl.glFrontFace(GL_CW);
		gl.glCullFace(GL_BACK);
		gl.glDrawArrays(GL_TRIANGLES, 0, vbo_index_num);
	}
	else if (mode == LineFrame)
	{
		gl.glPolygonMode(GL_FRONT, GL_LINE);
		gl.glDrawArrays(GL_LINES, 0, vbo_index_num);
		gl.glPolygonMode(GL_FRONT, GL_FILL);
	}
}

void QtRigidTetrahedronMeshGLObject::form_model_mat(
	const Point3D& cen,
	const Vector3D& ix,
	const Vector3D& iy, 
	const Vector3D& iz
	)
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
