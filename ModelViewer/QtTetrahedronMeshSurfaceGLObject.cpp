#include "ModelViewer_pcp.h"

#include "QtTetrahedronMeshSurfaceGLObject.h"

QtTetrahedronMeshSurfaceGLObject::QtTetrahedronMeshSurfaceGLObject(QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), vao(0), vbo(0), face_node_num(0), color(1.0f, 1.0f, 1.0f) {}

QtTetrahedronMeshSurfaceGLObject::~QtTetrahedronMeshSurfaceGLObject() { clear(); }

void QtTetrahedronMeshSurfaceGLObject::clear()
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
	face_node_num = 0;
}

void QtTetrahedronMeshSurfaceGLObject::draw(QOpenGLShaderProgram& shader)
{
	if (vao)
	{
		//gl.glFrontFace(GL_CW);
		//gl.glCullFace(GL_BACK);
		shader.setUniformValue("g_color", color);
		gl.glBindVertexArray(vao);
		gl.glDrawArrays(GL_TRIANGLES, 0, face_node_num);
	}
}
