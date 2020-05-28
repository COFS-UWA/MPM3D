#ifndef __Multi_Color_Particle_Buffer_h__
#define __Multi_Color_Particle_Buffer_h__

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

#include "ValueToColor.h"
#include "MultiColorCubeParticleSystem.h"

class MultiColorParticleBuffer
{
protected:
	// opengl functions wrapper
	QOpenGLFunctions_3_3_Core& glwp;

	GLuint vao, vbo, veo;
	GLsizei index_num;
	MultiColorCubeParticleSystem pcl_sys;

public:
	MultiColorParticleBuffer(QOpenGLFunctions_3_3_Core& wp) : glwp(wp),
		vao(0), vbo(0), veo(0), index_num(0) {}
	~MultiColorParticleBuffer() { clear(); }
	void clear()
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

	void draw()
	{
		if (vao)
		{
			glwp.glBindVertexArray(vao);
			glwp.glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glwp.glDrawElements(GL_TRIANGLES, index_num, GL_UNSIGNED_INT, 0);
		}
	}

	template <typename FieldType>
	int init_data(
		char* pcls_data, size_t pcl_size, size_t pcl_num,
		size_t x_off, size_t y_off, size_t z_off,
		size_t vol_off, float vol_scale,
		size_t fld_off, ValueToColor& v2c
		)
	{
		clear();
		if (!pcls_data || pcl_num == 0)
			return -1;

		glwp.glGenVertexArrays(1, &vao);
		glwp.glBindVertexArray(vao);

		pcl_sys.init_data<FieldType>(
			pcls_data, pcl_size, pcl_num, 
			x_off, y_off, z_off,
			vol_off, vol_scale,
			fld_off, v2c
			);

		// node data
		glwp.glGenBuffers(1, &vbo);
		glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glwp.glBufferData(
			GL_ARRAY_BUFFER,
			pcl_sys.get_vert_data_size(),
			pcl_sys.get_vert_data(),
			GL_STREAM_DRAW
			);

		glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat), (GLvoid*)0);
		glwp.glEnableVertexAttribArray(0);
		glwp.glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GLfloat), (GLvoid*)(3*sizeof(GLfloat)));
		glwp.glEnableVertexAttribArray(1);

		// element data
		glwp.glGenBuffers(1, &veo);
		glwp.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo);
		index_num = pcl_sys.get_elem_data_size() / sizeof(GLuint);
		glwp.glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			pcl_sys.get_elem_data_size(),
			pcl_sys.get_elem_data(),
			GL_STREAM_DRAW
			);

		return 0;
	}

	template <typename FieldType>
	int update_data(char* pcls_data)
	{
		int res;
		res = pcl_sys.update_data<FieldType>(pcls_data);
		if (res < 0)
			return res;

		glwp.glBindVertexArray(vao);
		
		glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glwp.glBufferData(
			GL_ARRAY_BUFFER,
			pcl_sys.get_vert_data_size(),
			pcl_sys.get_vert_data(),
			GL_STREAM_DRAW
			);

		return 0;
	}
};

#endif