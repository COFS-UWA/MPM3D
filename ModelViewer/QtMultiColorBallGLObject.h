#ifndef __Qt_Multi_Color_Ball_GL_Object_h__
#define __Qt_Multi_Color_Ball_GL_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"

class QtMultiColorBallGLObject
{
protected:
	QOpenGLFunctions_3_3_Core& gl;

	struct PointData
	{
		GLuint type;
		GLfloat x, y, z;
		GLfloat radius;
		GLfloat fld;
	};
	MemoryUtils::ItemArray<PointData> pt_data_mem;

	GLuint vao, vbo_cs, veo_cs, vbo_pts;
	size_t c_elem_node_num;
	size_t pt_num;

	void clear();

	int init_ball_data();
	int init_gl_buffer(PointData *pds, size_t pd_num);
	int update_gl_buffer(PointData *pds, size_t pd_num);

public:
	QtMultiColorBallGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtMultiColorBallGLObject();

	int init(size_t pcl_num, double *pcl_x_data, double *pcl_y_data, double *pcl_z_data,
		double *pcl_vol_data, double *pcl_fld_data, float radius_scale);

	int update(size_t pcl_num, double* pcl_x_data, double* pcl_y_data, double* pcl_z_data,
		double* pcl_vol_data, double *pcl_fld_data, float radius_scale);

	template <class DivisionSet>
	int init(size_t pcl_num, double* pcl_x_data, double* pcl_y_data, double* pcl_z_data,
		double* pcl_vol_data, double* pcl_fld_data,	float radius_scale, DivisionSet& div_set);

	template <class DivisionSet>
	int update(size_t pcl_num, double* pcl_x_data, double* pcl_y_data, double* pcl_z_data,
		double* pcl_vol_data, double* pcl_fld_data, float radius_scale, DivisionSet& div_set);

	void draw(QOpenGLShaderProgram& shader);
};

template <class DivisionSet>
int QtMultiColorBallGLObject::init(
	size_t pcl_num,
	double* pcl_x_data,
	double* pcl_y_data,
	double* pcl_z_data,
	double* pcl_vol_data,
	double* pcl_fld_data,
	float radius_scale,
	DivisionSet& div_set
	)
{
	pt_data_mem.reserve(pcl_num);
	PointData* pt_data = pt_data_mem.get_mem();
	size_t cur_p_id = 0;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		if (!div_set.inside(
			pcl_x_data[pcl_id],
			pcl_y_data[pcl_id],
			pcl_z_data[pcl_id]))
		{
			PointData& pd = pt_data[cur_p_id++];
			pd.type = 1; // multi color
			pd.x = GLfloat(pcl_x_data[pcl_id]);
			pd.y = GLfloat(pcl_y_data[pcl_id]);
			pd.z = GLfloat(pcl_z_data[pcl_id]);
			pd.radius = GLfloat(pow(3.0 * pcl_vol_data[pcl_id] / (4.0 * 3.14159265359), 1.0/3.0)) * radius_scale;
			pd.fld = GLfloat(pcl_fld_data[pcl_id]);
		}
	}
	int res = init_gl_buffer(pt_data, cur_p_id);
	return res;
}

template <class DivisionSet>
int QtMultiColorBallGLObject::update(
	size_t pcl_num,
	double* pcl_x_data,
	double* pcl_y_data,
	double* pcl_z_data,
	double* pcl_vol_data,
	double* pcl_fld_data,
	float radius_scale,
	DivisionSet& div_set
	)
{
	pt_data_mem.reserve(pcl_num);
	PointData* pt_data = pt_data_mem.get_mem();
	size_t cur_p_id = 0;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		if (!div_set.inside(
			pcl_x_data[pcl_id],
			pcl_y_data[pcl_id],
			pcl_z_data[pcl_id])
			)
		{
			PointData& pd = pt_data[cur_p_id++];
			pd.type = 1; // multi color
			pd.x = GLfloat(pcl_x_data[pcl_id]);
			pd.y = GLfloat(pcl_y_data[pcl_id]);
			pd.z = GLfloat(pcl_z_data[pcl_id]);
			pd.radius = GLfloat(pow(3.0 * pcl_vol_data[pcl_id] / (4.0 * 3.14159265359), 1.0/3.0)) * radius_scale;
			pd.fld = GLfloat(pcl_fld_data[pcl_id]);
		}
	}
	int res = update_gl_buffer(pt_data, cur_p_id);
	return res;
}

#endif