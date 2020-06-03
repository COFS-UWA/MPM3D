#include "ModelViewer_pcp.h"

#include "ball_mesh_data.h"
#include "BallParticleBuffer.h"

BallParticleBuffer::BallParticleBuffer(QOpenGLFunctions_3_3_Core& wp) :
	glwp(wp),
	ball_vao(0), ball_vbo(0), ball_veo(0), index_num(0),
	inst_vbo(0), inst_num(0)
{

}

BallParticleBuffer::~BallParticleBuffer()
{
	clear_ball_buffer();
	clear_inst_buffer();
}

void BallParticleBuffer::gen_ball_buffer()
{
	if (ball_vao)
		return;

	glwp.glGenVertexArrays(1, &ball_vao);
	glwp.glBindVertexArray(ball_vao);

	// node data
	glwp.glGenBuffers(1, &ball_vbo);
	glwp.glBindBuffer(GL_ARRAY_BUFFER, ball_vbo);
	glwp.glBufferData(
		GL_ARRAY_BUFFER,
		ball_node_num * 3 * sizeof(GLuint),
		ball_nodes,
		GL_STREAM_DRAW
		);

	// coordinates / normal
	glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)0);
	glwp.glEnableVertexAttribArray(0);

	// element data
	index_num = ball_elem_num * 3;
	glwp.glGenBuffers(1, &ball_veo);
	glwp.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ball_veo);
	glwp.glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		ball_elem_num * 3 * sizeof(GLuint),
		ball_elems,
		GL_STREAM_DRAW
		);
}

void BallParticleBuffer::clear_ball_buffer()
{
	if (ball_veo)
	{
		glwp.glDeleteBuffers(1, &ball_veo);
		ball_veo = 0;
	}
	if (ball_vbo)
	{
		glwp.glDeleteBuffers(1, &ball_vbo);
		ball_vbo = 0;
	}
	if (ball_vao)
	{
		glwp.glDeleteVertexArrays(1, &ball_vao);
		ball_vao = 0;
	}
	index_num = 0;
}

void BallParticleBuffer::clear_inst_buffer()
{
	if (inst_vbo)
	{
		glwp.glDeleteBuffers(1, &inst_vbo);
		inst_vbo = 0;
	}
	inst_num = 0;
}

void BallParticleBuffer::draw()
{
	if (ball_vao && inst_vbo)
	{
		glwp.glPolygonMode(GL_FRONT, GL_FILL);
		glwp.glEnable(GL_CULL_FACE);
		glwp.glBindVertexArray(ball_vao);
		glwp.glDrawArraysInstanced(GL_TRIANGLES, 0, index_num, inst_num);
	}
}
