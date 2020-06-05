#include "ModelViewer_pcp.h"

#include <cstddef>

#include "ball_mesh_data.h"
#include "BallParticleBuffer.h"

BallParticleBuffer::BallParticleBuffer(QOpenGLFunctions_3_3_Core& wp) :
	glwp(wp),
	ball_vao(0), ball_vbo(0), ball_veo(0), index_num(0),
	inst_vbo(0), inst_num(0) {}

BallParticleBuffer::~BallParticleBuffer()
{
	clear_inst_buffer();
	clear_ball_buffer();
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
		sizeof(ball_nodes),
		ball_nodes,
		GL_STREAM_DRAW
	);

	// coordinates / normal
	glwp.glEnableVertexAttribArray(0);
	glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);

	// element data
	glwp.glGenBuffers(1, &ball_veo);
	glwp.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ball_veo);
	glwp.glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		sizeof(ball_elems),
		ball_elems,
		GL_STREAM_DRAW
	);

	index_num = ball_elem_num * 3;
}

void BallParticleBuffer::gen_inst_buffer(
	GLvoid* inst_buf,
	size_t inst_buf_size,
	size_t _inst_num
	)
{
	gen_ball_buffer();

	glwp.glBindVertexArray(ball_vao);

	glwp.glGenBuffers(1, &inst_vbo);
	glwp.glBindBuffer(GL_ARRAY_BUFFER, inst_vbo);
	glwp.glBufferData(GL_ARRAY_BUFFER,
		inst_buf_size,
		inst_buf,
		GL_STREAM_DRAW
	);

	// pcl coordinates
	glwp.glEnableVertexAttribArray(1);
	glwp.glVertexAttribPointer(1,
		3, GL_FLOAT, GL_FALSE,
		sizeof(BallInstance), (GLvoid*)0
		);
	glwp.glVertexAttribDivisor(1, 1);
	// radius
	glwp.glEnableVertexAttribArray(2);
	glwp.glVertexAttribPointer(2,
		1, GL_FLOAT, GL_FALSE,
		sizeof(BallInstance),
		(GLvoid *)offsetof(BallInstance, radius)
		);
	glwp.glVertexAttribDivisor(2, 1);
	// color
	glwp.glEnableVertexAttribArray(3);
	glwp.glVertexAttribPointer(3,
		3, GL_FLOAT, GL_FALSE,
		sizeof(BallInstance),
		(GLvoid *)offsetof(BallInstance, color)
		);
	glwp.glVertexAttribDivisor(3, 1);

	inst_num = _inst_num;
}

void BallParticleBuffer::draw()
{
	if (ball_vao)
	{
		glwp.glEnable(GL_CULL_FACE);
		glwp.glCullFace(GL_BACK);
		glwp.glPolygonMode(GL_FRONT, GL_FILL);
		glwp.glBindVertexArray(ball_vao);
		glwp.glDrawElementsInstanced(
			GL_TRIANGLES,
			index_num,
			GL_UNSIGNED_INT,
			0, inst_num
			);
	}
}
