#include "ModelViewer_pcp.h"

#include "QtTriangleMeshGLObject.h"

QtTriangleMeshGLObject::QtTriangleMeshGLObject(QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), vao(0), vbo(0), veo(0),
	edge_node_num(0) {}

QtTriangleMeshGLObject::~QtTriangleMeshGLObject() { clear(); }

void QtTriangleMeshGLObject::clear()
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
}

void QtTriangleMeshGLObject::draw(QOpenGLShaderProgram& shader)
{
	shader.setUniformValue("g_color", color);

	gl.glLineWidth(2.0f);

	gl.glBindVertexArray(vao);
	gl.glDrawElements(GL_LINES, edge_node_num, GL_UNSIGNED_INT, (GLvoid *)0);
}
