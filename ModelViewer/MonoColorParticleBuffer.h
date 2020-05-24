#ifndef __Mono_Color_Particle_Buffer_h__
#define __Mono_Color_Particle_Buffer_h__

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

#include "MonoColorCubeParticleSystem.h"

class MonoColorParticleBuffer
{
protected:
	// opengl functions wrapper
	QOpenGLFunctions_3_3_Core &glwp;

	GLuint vao, vbo, veo;
	GLsizei index_num;
	QVector3D color;
	MonoColorCubeParticleSystem pcl_sys;

public:
	MonoColorParticleBuffer(QOpenGLFunctions_3_3_Core &wp) : glwp(wp),
		vao(0), vbo(0), veo(0), index_num(0), color(1.0f, 1.0f, 1.0f) {}
	~MonoColorParticleBuffer() { clear(); }
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

	// use uniform color shader
	void draw(QOpenGLShaderProgram& shader)
	{
		if (vao)
		{
			shader.setUniformValue("color", color);

			glwp.glBindVertexArray(vao);
			glwp.glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			glwp.glDrawElements(GL_TRIANGLES, index_num, GL_UNSIGNED_INT, 0);
		}
	}

	template <typename Particle>
	int init_pcl_data(Particle *pcls, size_t pcl_num, QVector3D &_color)
	{
		clear();
		if (!pcls || pcl_num == 0)
			return -1;

		color = _color;

		glwp.glGenVertexArrays(1, &vao);
		glwp.glBindVertexArray(vao);

		pcl_sys.init_data<Particle>(pcls, pcl_num);

		size_t aa;
		aa = pcl_sys.get_vert_data_size();
		aa = pcl_sys.get_elem_data_size();

		// node data
		glwp.glGenBuffers(1, &vbo);
		glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glwp.glBufferData(
			GL_ARRAY_BUFFER,
			pcl_sys.get_vert_data_size(),
			pcl_sys.get_vert_data(),
			GL_STREAM_DRAW
			);

		glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid *)0);
		glwp.glEnableVertexAttribArray(0);

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

	template <typename Particle>
	int update_pcl_data(Particle* pcls)
	{
		if (!pcls)
			return -1;

		glwp.glBindVertexArray(vao);

		pcl_sys.update_data<Particle>(pcls);

		glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glwp.glBufferData(
			GL_ARRAY_BUFFER,
			pcl_sys.get_vert_data_size(),
			pcl_sys.get_vert_data(),
			GL_STREAM_DRAW
		);

		// need?
		//glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
		//glwp.glEnableVertexAttribArray(0);

		return 0;
	}
};

#endif