#include "ModelViewer_pcp.h"

#include "QtTetrahedronMeshGLObject.h"

QtTetrahedronMeshGLObject::QtTetrahedronMeshGLObject(QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), vao(0), vbo(0), veo(0), index_num(0), color(1.0f, 1.0f, 1.0f) {}

QtTetrahedronMeshGLObject::~QtTetrahedronMeshGLObject() { clear(); }

void QtTetrahedronMeshGLObject::clear()
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

	index_num = 0;
}

void QtTetrahedronMeshGLObject::draw(QOpenGLShaderProgram& shader)
{
	if (vao)
	{
		shader.setUniformValue("color", color);
		gl.glBindVertexArray(vao);
		gl.glLineWidth(1.0f);
		gl.glDrawElements(GL_LINES, index_num, GL_UNSIGNED_INT, 0);
	}
}
