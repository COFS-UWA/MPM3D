#ifndef __Qt_Multi_Color_Circle_GL_Object_h__
#define __Qt_Multi_Color_Circle_GL_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"

class QtMultiColorCircleGLObject
{
protected:
	QOpenGLFunctions_3_3_Core& gl;

	struct PointData
	{
		GLuint type;
		GLfloat x, y;
		GLfloat radius;
		GLfloat fld;
	};
	MemoryUtils::ItemArray<PointData> pt_data_mem;

	GLuint vao, vbo_cs, veo_cs, vbo_pts;
	size_t c_elem_node_num;
	size_t pt_num;

	void clear();

	int init_circle_data();
	int init_gl_buffer(PointData *pds, size_t pd_num);
	int update_gl_buffer(PointData *pds, size_t pd_num);

public:
	QtMultiColorCircleGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtMultiColorCircleGLObject();

	int init(size_t pcl_num, double* pcl_x_data, double* pcl_y_data,
		double* pcl_vol_data, double* pcl_fld_data, float radius_scale);

	int update(size_t pcl_num, double* pcl_x_data, double* pcl_y_data,
		double* pcl_vol_data, double* pcl_fld_data, float radius_scale);

	int init(size_t pcl_num, float* pcl_x_data, float* pcl_y_data,
		float *pcl_vol_data, float* pcl_fld_data, float radius_scale);

	int update(size_t pcl_num, float *pcl_x_data, float *pcl_y_data,
		float *pcl_vol_data, float *pcl_fld_data, float radius_scale);
	
	void draw(QOpenGLShaderProgram& shader);
};

#endif