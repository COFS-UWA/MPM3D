#include "ModelViewer_pcp.h"

#include "PointBuffer.h"

void PointBuffer::clear()
{
	if (vao)
	{
		glwp.glDeleteVertexArrays(1, &vao);
		vao = 0;
	}
	if (vbo)
	{
		glwp.glDeleteBuffers(1, &vbo);
		vbo = 0;
	}
	point_num = 0;
}

int PointBuffer::init(
	Point3D* pts,
	size_t pt_num,
	GLfloat pt_size,
	QVector3D& _color
	)
{
	clear();
	if (!pts || pt_num == 0)
		return -1;

	point_size = pt_size;
	color = _color;

	glwp.glGenVertexArrays(1, &vao);
	glwp.glBindVertexArray(vao);

	glwp.glGenBuffers(1, &vbo);
	glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
	point_num = pt_num;
	GLfloat* coords = new GLfloat[point_num * 3];
	GLfloat* p = coords;
	for (size_t p_id = 0; p_id < point_num; ++p_id)
	{
		Point3D& pt = pts[p_id];
		p[0] = GLfloat(pt.x);
		p[1] = GLfloat(pt.y);
		p[2] = GLfloat(pt.z);
		p += 3;
	}
	glwp.glBufferData(
		GL_ARRAY_BUFFER,
		point_num * 3 * sizeof(GLfloat),
		coords,
		GL_STREAM_DRAW
		);
	delete[] coords;

	glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	glwp.glEnableVertexAttribArray(0);

	return 0;
}

