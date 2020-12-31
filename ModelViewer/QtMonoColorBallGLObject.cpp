#include "ModelViewer_pcp.h"

#include <cstddef>
#include <iostream>

#include "ball_mesh_data.h"
#include "QtMonoColorBallGLObject.h"

QtMonoColorBallGLObject::QtMonoColorBallGLObject(QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), vao(0), vbo_cs(0), veo_cs(0), vbo_pts(0),
	c_elem_node_num(0), pt_num(0),
	max_pcl_num_per_drawcall(1500000), pt_data(nullptr) {}

QtMonoColorBallGLObject::~QtMonoColorBallGLObject() { clear(); }

void QtMonoColorBallGLObject::clear()
{
	if (veo_cs)
	{
		gl.glDeleteBuffers(1, &veo_cs);
		veo_cs = 0;
	}
	if (vbo_cs)
	{
		gl.glDeleteBuffers(1, &vbo_cs);
		vbo_cs = 0;
	}
	if (vbo_pts)
	{
		gl.glDeleteBuffers(1, &vbo_pts);
		vbo_pts = 0;
	}
	if (vao)
	{
		gl.glDeleteVertexArrays(1, &vao);
		vao = 0;
	}
	if (pt_data)
	{
		delete[] pt_data;
		pt_data = nullptr;
	}
}

void QtMonoColorBallGLObject::draw(QOpenGLShaderProgram& shader)
{
	shader.bind();
	shader.setUniformValue("g_color", color);

	gl.glBindVertexArray(vao);

	if (pt_num > max_pcl_num_per_drawcall)
	{
		gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_pts);
		size_t pt_start_id = 0;
		size_t pt_end_id = pt_start_id + max_pcl_num_per_drawcall;
		while (pt_end_id < pt_num)
		{
			gl.glBufferSubData(GL_ARRAY_BUFFER, 0,
				max_pcl_num_per_drawcall * sizeof(PointData),
				(GLvoid*)(pt_data + pt_start_id));
			gl.glDrawElementsInstanced(
				GL_TRIANGLES,
				c_elem_node_num,
				GL_UNSIGNED_INT,
				(GLvoid*)0,
				max_pcl_num_per_drawcall);
			gl.glFinish();
			pt_start_id = pt_end_id;
			pt_end_id += max_pcl_num_per_drawcall;
		}
		gl.glBufferSubData(GL_ARRAY_BUFFER, 0,
			(pt_num - pt_start_id) * sizeof(PointData),
			(GLvoid*)(pt_data + pt_start_id));
		gl.glDrawElementsInstanced(
			GL_TRIANGLES,
			c_elem_node_num,
			GL_UNSIGNED_INT,
			(GLvoid*)0,
			pt_num - pt_start_id);
		gl.glFinish();
	}
	else
	{
		gl.glDrawElementsInstanced(
			GL_TRIANGLES,
			c_elem_node_num,
			GL_UNSIGNED_INT,
			(GLvoid*)0,
			pt_num);
	}
}

static const char* get_gl_error_info(GLenum errorId)
{
	switch (errorId)
	{
	case GL_INVALID_ENUM:
		return ("GL Invalid Enum\n");
	case GL_INVALID_VALUE:
		return ("GL Invalid Value\n");
	case GL_INVALID_OPERATION:
		return ("GL Invalid Operation\n");
	case GL_OUT_OF_MEMORY:
		return ("GL Out Of Memory\n");
		//case GL_INVALID_FRAMEBUFFER_OPERATION:
		//	return ("GL Invalid FrameBuffer Operation\n");
	case  GL_STACK_OVERFLOW:
		return ("GL Stack Overflow\n");
	case GL_STACK_UNDERFLOW:
		return ("GL Stack Underflow\n");
		//case GL_TABLE_TOO_LARGE:
		//	return ("GL Table Too Large\n");
	};

	return "GL Undefined Error";
}

int QtMonoColorBallGLObject::init_gl_buffer(
	PointData *pds,
	size_t pd_num
	)
{
	pt_num = pd_num;
	gl.glGenVertexArrays(1, &vao);
	if (vao == 0)
		return -1;
	gl.glBindVertexArray(vao);

	int res = init_ball_data();
	if (res)
		return res;

	gl.glGenBuffers(1, &vbo_pts);
	if (vbo_pts == 0)
		return -1;
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_pts);
	gl.glBufferData(GL_ARRAY_BUFFER,
		pd_num * sizeof(PointData),
		pds, GL_STREAM_DRAW);

	//for (GLenum err; (err = gl.glGetError()) != GL_NO_ERROR;)
	//	std::cout << get_gl_error_info(err) << "\n";

	// pt_type
	gl.glVertexAttribIPointer(1,
		1, GL_UNSIGNED_INT,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, type)
		);
	gl.glEnableVertexAttribArray(1);
	gl.glVertexAttribDivisor(1, 1);
	// pt_pos
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, x)
		);
	gl.glEnableVertexAttribArray(2);
	gl.glVertexAttribDivisor(2, 1);
	// pt_radius
	gl.glVertexAttribPointer(3,
		1, GL_FLOAT, GL_FALSE,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, radius)
		);
	gl.glEnableVertexAttribArray(3);
	gl.glVertexAttribDivisor(3, 1);

	return 0;
}

int QtMonoColorBallGLObject::init_ball_data()
{
	gl.glGenBuffers(1, &vbo_cs);
	if (vbo_cs == 0)
		return -1;
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_cs);
	gl.glBufferData(GL_ARRAY_BUFFER,
		sizeof(ball_nodes),
		ball_nodes,
		GL_STATIC_DRAW
		);

	gl.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
	gl.glEnableVertexAttribArray(0);

	gl.glGenBuffers(1, &veo_cs);
	if (veo_cs == 0)
		return -1;
	gl.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo_cs);
	gl.glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		sizeof(ball_elems),
		ball_elems,
		GL_STATIC_DRAW
		);

	c_elem_node_num = sizeof(ball_elems) / sizeof(ball_elems[0]);

	return 0;
}
