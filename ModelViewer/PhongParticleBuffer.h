#ifndef __Phong_Particle_Buffer_h__
#define __Phong_Particle_Buffer_h__

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

#include "PhongCubeParticleSystem.h"

class PhongParticleBuffer
{
protected:
	// opengl functions wrapper
	QOpenGLFunctions_3_3_Core &glwp;

	GLuint vao, vbo, veo;
	GLsizei index_num;
	PhongCubeParticleSystem pcl_sys;

	typedef ValueToColor::Colorf Color;

public:
	PhongParticleBuffer(QOpenGLFunctions_3_3_Core& wp) :
		glwp(wp), vao(0), vbo(0), veo(0), index_num(0) {}
	~PhongParticleBuffer() { clear(); }
	void clear();

	void draw();

protected:
	int create_buffer(GLvoid *vbo_data, GLsizei vbo_size,
					  GLvoid *veo_data, GLsizei veo_size);
	int update_vbo_buffer(GLvoid* vbo_data, GLsizei vbo_size);

public:
	// monocolor particle version
	int init_data(char *pcls_data, size_t pcl_size, size_t pcl_num,
				  size_t x_off, size_t y_off, size_t z_off,
				  size_t vol_off, float vol_scale, Color &color);
	int update_data(char* pcls_data, Color& color);

	int init_points(Point3D* pcls, size_t pcl_num,
					float pcl_vol, Color& color);
	
	template <typename Particle>
	int init_data(Particle *pcls, size_t pcl_num,
				  float vol_scale, Color &color)
	{
		clear();
		if (!pcls || pcl_num == 0)
			return -1;

		pcl_sys.init_data<Particle>(
			pcls, pcl_num,
			vol_scale, color 
			);

		create_buffer(
			pcl_sys.get_vert_data(),
			pcl_sys.get_vert_data_size(),
			pcl_sys.get_elem_data(),
			pcl_sys.get_elem_data_size()
		);
		
		return 0;
	}

	// colorful particle version
	template <typename FieldType>
	int init_data(
		char* pcls_data, size_t pcl_size, size_t pcl_num,
		size_t x_off, size_t y_off, size_t z_off,
		size_t vol_off, float vol_scale,
		size_t fld_off, ValueToColor &v2c
		)
	{
		clear();
		if (!pcls_data || pcl_num == 0)
			return -1;

		pcl_sys.init_data<FieldType>(
			pcls_data, pcl_size, pcl_num,
			x_off, y_off, z_off,
			vol_off, vol_scale,
			fld_off, v2c
			);

		create_buffer(
			pcl_sys.get_vert_data(),
			pcl_sys.get_vert_data_size(),
			pcl_sys.get_elem_data(),
			pcl_sys.get_elem_data_size()
		);

		return 0;
	}

	template <typename FieldType>
	int update_data(char *pcls_data, ValueToColor& v2c)
	{
		pcl_sys.update_data<FieldType>(pcls_data, v2c);

		update_vbo_buffer(
			pcl_sys.get_vert_data(),
			pcl_sys.get_vert_data_size()
		);

		return 0;
	}

};

#endif