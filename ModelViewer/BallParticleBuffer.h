#ifndef __Ball_Particle_Buffer_h__
#define __Ball_Particle_Buffer_h__

#include <iostream>

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

#include "ItemArray.hpp"
#include "ValueToColor.h"

class BallParticleBuffer
{
protected:
	// opengl functions wrapper
	QOpenGLFunctions_3_3_Core &glwp;

	GLuint ball_vao, ball_vbo, ball_veo;
	GLsizei index_num;

	GLuint inst_vbo;
	GLsizei inst_num;

	typedef ValueToColor::Colorf Color;
	
	struct BallInstance
	{
		GLfloat x, y, z;
		GLfloat radius;
		union // color
		{
			struct { GLfloat r, g, b; };
			Color color;
		};
	};
	MemoryUtils::ItemArray<BallInstance> ball_inst_buf;

	size_t pcl_size, pcl_num;
	size_t x_off, y_off, z_off;
	size_t vol_off;
	float vol_scale;
	size_t fld_off;

public:
	BallParticleBuffer(QOpenGLFunctions_3_3_Core& wp);
	~BallParticleBuffer();

	void gen_ball_buffer();
	void clear_ball_buffer();

	void gen_inst_buffer(
		GLvoid *inst_buf, size_t inst_buf_size,
		size_t _inst_num);
	void clear_inst_buffer();

	void draw();

	// mono color particles
	template <typename Particle>
	int init_data(Particle* pcls, size_t _pcl_num,
				  float _vol_scale, Color& _color)
	{
		clear_inst_buffer();
		if (!pcls || _pcl_num == 0)
			return -1;

		pcl_num = _pcl_num;
		vol_scale = _vol_scale;

		GLfloat pcl_vol;
		ball_inst_buf.reserve(pcl_num);
		BallInstance* ball_insts = ball_inst_buf.get_mem();
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			Particle& pcl = pcls[pcl_id];
			BallInstance &bi = ball_insts[pcl_id];
			bi.x = (GLfloat)pcl.x;
			bi.y = (GLfloat)pcl.y;
			bi.z = (GLfloat)pcl.z;
			pcl_vol = (GLfloat)pcl.get_vol();
			bi.radius = pow(pcl_vol * vol_scale * 3.0f / (4.0f * 3.14159265359f), 0.3333333);
			bi.color = _color;
		}

		gen_inst_buffer(ball_insts, pcl_num * sizeof(BallInstance), pcl_num);

		return 0;
	}

	// colorful particles
	template <typename FieldType>
	int init_data(
		char *pcls_data, size_t _pcl_size, size_t _pcl_num,
		size_t _x_off, size_t _y_off, size_t _z_off,
		size_t _vol_off, float _vol_scale,
		size_t _fld_off, ValueToColor &v2c
		)
	{
		clear_inst_buffer();
		if (!pcls_data || _pcl_num == 0)
			return -1;

		pcl_size = _pcl_size;
		pcl_num = _pcl_num;
		x_off = _x_off;
		y_off = _y_off;
		z_off = _z_off;
		vol_off = _vol_off;
		vol_scale = _vol_scale;
		fld_off = _fld_off;
		
		char *pcl_data = pcls_data;
		double pcl_vol, pcl_fld;
		ball_inst_buf.reserve(pcl_num);
		BallInstance *ball_insts = ball_inst_buf.get_mem();
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			BallInstance& bi = ball_insts[pcl_id];
			bi.x = *(double*)(pcl_data + x_off);
			bi.y = *(double*)(pcl_data + y_off);
			bi.z = *(double*)(pcl_data + z_off);
			pcl_vol = *(double*)(pcl_data + vol_off);
			bi.radius = pow(pcl_vol*vol_scale * 3.0f/(4.0f*3.14159265359f), 0.3333333);
			pcl_fld = double(*(FieldType*)(pcl_data + fld_off));
			bi.color = v2c.get_color(pcl_fld);
			pcl_data += pcl_size;
		}

		gen_inst_buffer(ball_insts, pcl_num * sizeof(BallInstance), pcl_num);
		
		return 0;
	}

	template <typename FieldType>
	int update_data(char *pcls_data, ValueToColor &v2c)
	{
		if (!pcls_data || pcl_num == 0 || !inst_vbo)
			return -1;

		char* pcl_data = pcls_data;
		double pcl_vol, pcl_fld;
		BallInstance* ball_insts = ball_inst_buf.get_mem();
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			BallInstance& bi = ball_insts[pcl_id];
			bi.x = *(double*)(pcl_data + x_off);
			bi.y = *(double*)(pcl_data + y_off);
			bi.z = *(double*)(pcl_data + z_off);
			pcl_vol = *(double*)(pcl_data + vol_off);
			bi.radius = pow(pcl_vol * vol_scale * 3.0f / (4.0f * 3.14159265359f), 0.3333333);
			pcl_fld = double(*(FieldType*)(pcl_data + fld_off));
			bi.color = v2c.get_color(pcl_fld);
			pcl_data += pcl_size;
		}

		glwp.glBindVertexArray(ball_vao);
		glwp.glBindBuffer(GL_ARRAY_BUFFER, inst_vbo);
		glwp.glBufferData(GL_ARRAY_BUFFER,
			pcl_num * sizeof(BallInstance),
			ball_insts,
			GL_STREAM_DRAW
		);

		return 0;
	}

};

#endif