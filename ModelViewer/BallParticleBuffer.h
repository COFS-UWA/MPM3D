#ifndef __Ball_Particle_Buffer_h__
#define __Ball_Particle_Buffer_h__

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
	void clear_inst_buffer();

	void draw();

	// colorful particle version
	template <typename FieldType>
	int init_data(
		char *pcls_data, size_t _pcl_size, size_t _pcl_num,
		size_t _x_off, size_t _y_off, size_t _z_off,
		size_t _vol_off, float _vol_scale,
		size_t _fld_off, ValueToColor &v2c
		)
	{
		clear_inst_buffer();
		if (!pcls_data || pcl_num == 0)
			return -1;
		
		gen_ball_buffer();

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
			bi.Color = v2c.get_color(pcl_fld);
			pcl_data += pcl_size;
		}

		glwp.glGenBuffers(1, &inst_vbo);
		glwp.glBindBuffer(GL_ARRAY_BUFFER, inst_vbo);
		glwp.glBufferData(GL_ARRAY_BUFFER,
			pcl_num * sizeof(BallInstance),
			ball_insts,
			GL_STREAM_DRAW
			);

		// coordinates
		glwp.glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(BallInstance), (GLvoid *)0);
		glwp.glEnableVertexAttribArray(1);
		glwp.glVertexAttribDivisor(1, 1);
		// radius
		glwp.glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(BallInstance), offsetof(BallInstance, radius));
		glwp.glEnableVertexAttribArray(3);
		glwp.glVertexAttribDivisor(2, 1);
		// color
		glwp.glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(BallInstance), offsetof(BallInstance, color));
		glwp.glEnableVertexAttribArray(4);
		glwp.glVertexAttribDivisor(3, 1);

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
			bi.Color = v2c.get_color(pcl_fld);
			pcl_data += pcl_size;
		}

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