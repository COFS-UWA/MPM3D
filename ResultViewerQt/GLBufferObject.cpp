#include "ResultViewerQt_pcp.h"

#include "GLBufferObject.h"

GLBufferObject::GLBufferObject(QOpenGLFunctions_3_3_Core *pa) :
	parent(*pa),
	vao(0), vbo(0), veo(0),
	vert_data_num(0), elem_data_num(0) {}

GLBufferObject::~GLBufferObject()
{
	clear();
}

void GLBufferObject::init_vert(GLfloat *node_data, size_t data_num)
{
	parent.glGenVertexArrays(1, &vao);
	parent.glBindVertexArray(vao);
	parent.glGenBuffers(1, &vbo);
	parent.glBindBuffer(GL_ARRAY_BUFFER, vbo);
	parent.glBufferData(
		GL_ARRAY_BUFFER,
		data_num * sizeof(GLfloat),
		node_data,
		GL_STATIC_DRAW
		);
}

void GLBufferObject::set_attr_pointer()
{
	parent.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 6, (GLvoid *)0);
	parent.glEnableVertexAttribArray(0);
}

void GLBufferObject::init_elem(GLuint *indices, size_t indices_num)
{
	parent.glBindVertexArray(vao);
	parent.glGenBuffers(1, &veo);
	parent.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo);
	parent.glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		indices_num * sizeof(GLuint),
		indices,
		GL_STATIC_DRAW
		);
}

void GLBufferObject::clear(void)
{
	if (vao)
	{
		parent.glDeleteVertexArrays(1, &vao);
		vao = 0;
	}
	if (vbo)
	{
		parent.glDeleteBuffers(1, &vbo);
		vbo = 0;
		vert_data_num = 0;
	}
	if (veo)
	{
		parent.glDeleteBuffers(1, &veo);
		veo = 0;
		elem_data_num = 0;
	}
}
