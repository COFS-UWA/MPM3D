#include "ModelViewer_pcp.h"

#include "PhongParticleBuffer.h"

void PhongParticleBuffer::clear()
{
	if (veo)
	{
		glwp.glDeleteBuffers(1, &veo);
		veo = 0;
	}
	if (vbo)
	{
		glwp.glDeleteBuffers(1, &vbo);
		vbo = 0;
	}
	if (vao)
	{
		glwp.glDeleteVertexArrays(1, &vao);
		vao = 0;
	}
	index_num = 0;
}

void PhongParticleBuffer::draw()
{
	if (vao)
	{
		glwp.glPolygonMode(GL_FRONT, GL_FILL);
		glwp.glEnable(GL_CULL_FACE);
		glwp.glCullFace(GL_BACK);
		glwp.glBindVertexArray(vao);
		glwp.glDrawElements(GL_TRIANGLES, index_num, GL_UNSIGNED_INT, 0);
	}
}

int PhongParticleBuffer::create_buffer(
	GLvoid* vbo_data,
	GLsizei vbo_size,
	GLvoid* veo_data,
	GLsizei veo_size
	)
{
	glwp.glGenVertexArrays(1, &vao);
	glwp.glBindVertexArray(vao);

	// node data
	glwp.glGenBuffers(1, &vbo);
	glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glwp.glBufferData(
		GL_ARRAY_BUFFER,
		vbo_size,
		vbo_data,
		GL_STREAM_DRAW
		);

	size_t v_off1 = sizeof(GLfloat) * 3;
	size_t v_off2 = sizeof(GLfloat) * 6;
	GLsizei v_strd = sizeof(GLfloat) * 9;
	// coordinates
	glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, v_strd, (GLvoid*)0);
	glwp.glEnableVertexAttribArray(0);
	// color
	glwp.glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, v_strd, (GLvoid*)v_off1);
	glwp.glEnableVertexAttribArray(1);
	// normal
	glwp.glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, v_strd, (GLvoid*)v_off2);
	glwp.glEnableVertexAttribArray(2);

	// element data
	glwp.glGenBuffers(1, &veo);
	glwp.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo);
	index_num = veo_size / sizeof(GLuint);
	glwp.glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		veo_size,
		veo_data,
		GL_STREAM_DRAW
		);
	
	return 0;
}

int PhongParticleBuffer::update_vbo_buffer(
	GLvoid* vbo_data,
	GLsizei vbo_size
	)
{
	glwp.glBindVertexArray(vao);

	glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glwp.glBufferData(
		GL_ARRAY_BUFFER,
		vbo_size,
		vbo_data,
		GL_STREAM_DRAW
		);
	
	return 0;
}

// monocolor particle version
int PhongParticleBuffer::init_data(
	char* pcls_data, size_t pcl_size, size_t pcl_num,
	size_t x_off, size_t y_off, size_t z_off,
	size_t vol_off, float vol_scale, Color& color
	)
{
	clear();
	if (!pcls_data || pcl_num == 0)
		return -1;

	pcl_sys.init_data(
		pcls_data, pcl_size, pcl_num,
		x_off, y_off, z_off,
		vol_off, vol_scale, color
	);

	create_buffer(
		pcl_sys.get_vert_data(),
		pcl_sys.get_vert_data_size(),
		pcl_sys.get_elem_data(),
		pcl_sys.get_elem_data_size()
	);

	return 0;
}

int PhongParticleBuffer::update_data(char* pcls_data, Color& color)
{
	pcl_sys.update_data(pcls_data, color);

	update_vbo_buffer(
		pcl_sys.get_vert_data(),
		pcl_sys.get_vert_data_size()
	);

	return 0;
}