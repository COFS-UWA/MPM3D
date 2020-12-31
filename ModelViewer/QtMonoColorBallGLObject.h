#ifndef __Qt_Mono_Color_Ball_GL_Object_h__
#define __Qt_Mono_Color_Ball_GL_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "DivisionSet.h"

class QtMonoColorBallGLObject
{
protected:
	QOpenGLFunctions_3_3_Core &gl;

	struct PointData
	{
		GLint type;
		GLfloat x, y, z;
		GLfloat radius;
	};

	GLuint vao, vbo_cs, veo_cs, vbo_pts;
	size_t c_elem_node_num;
	size_t pt_num;
	QVector3D color;
	size_t max_pcl_num_per_drawcall;
	PointData* pt_data;

	void clear();

	int init_ball_data();
	int init_gl_buffer(PointData *pds, size_t pd_num);

public:
	QtMonoColorBallGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtMonoColorBallGLObject();

	inline void set_max_pcl_num_per_drawcall(size_t _num) noexcept { max_pcl_num_per_drawcall = _num; }
	inline void set_max_gl_buffer_size(size_t _size) noexcept { max_pcl_num_per_drawcall = _size / sizeof(PointData); }

	// Particle3D has member x, y, z and vol
	template <typename Particle3D>
	int init(const Particle3D* pcls, size_t pcl_num, QVector3D& c, float radius_scale = 0.5f);

	template <typename Particle3D>
	int init(const Particle3D* pcls, const double* pcl_m, const double* pcl_density,
		size_t pcl_num, QVector3D& c, float radius_scale);

	template <typename DivisionSet>
	int init(const double *pcl_x, const double *pcl_y,
			 const double *pcl_z, const double* pcl_vol,
			 const size_t pcl_num, 
			 const QVector3D& c, const float radius_scale,
			 const DivisionSet& div_set);
	
	// uniform size particles
	template <typename Point3D>
	int init(Point3D* pts, size_t pt_num, float pcl_radius, QVector3D& c);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename Particle3D>
int QtMonoColorBallGLObject::init(
	const Particle3D *pcls,
	size_t pcl_num,
	QVector3D& c,
	float radius_scale
	)
{
	clear();
	color = c;
	pt_data = new PointData[pcl_num];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		const Particle3D &pcl = pcls[p_id];
		PointData& pd = pt_data[p_id];
		pd.type = 0; // mono color
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
		pd.z = GLfloat(pcl.z);
		pd.radius = GLfloat(pow(3.0 * pcl.get_vol()/ (4.0 * 3.14159265359), 0.3333333))	* radius_scale;
	}
	if (pcl_num > max_pcl_num_per_drawcall)
	{
		init_gl_buffer(pt_data, max_pcl_num_per_drawcall);
		pt_num = pcl_num;
	}
	else
	{
		init_gl_buffer(pt_data, pcl_num);
		delete[] pt_data;
		pt_data = nullptr;
	}
	return 0;
}

template <typename Particle3D>
int QtMonoColorBallGLObject::init(
	const Particle3D* pcls,
	const double *pcl_m,
	const double *pcl_density,
	size_t pcl_num,
	QVector3D& c,
	float radius_scale
	)
{
	clear();
	color = c;
	pt_data = new PointData[pcl_num];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		const Particle3D& pcl = pcls[p_id];
		PointData& pd = pt_data[p_id];
		pd.type = 0; // mono color
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
		pd.z = GLfloat(pcl.z);
		pd.radius = GLfloat(pow(3.0 * pcl_m[p_id] / pcl_density[p_id] / (4.0 * 3.14159265359), 0.3333333)) * radius_scale;
	}
	if (pcl_num > max_pcl_num_per_drawcall)
	{
		init_gl_buffer(pt_data, max_pcl_num_per_drawcall);
		pt_num = pcl_num;
	}
	else
	{
		init_gl_buffer(pt_data, pcl_num);
		delete[] pt_data;
		pt_data = nullptr;
	}
	return 0;
}

template <typename Point3D>
int QtMonoColorBallGLObject::init(
	Point3D* pts,
	size_t _pt_num,
	float pt_radius,
	QVector3D& c
)
{
	clear();
	color = c;
	pt_data = new PointData[_pt_num];
	for (size_t p_id = 0; p_id < _pt_num; ++p_id)
	{
		Point3D& pt = pts[p_id];
		PointData& pd = pt_data[p_id];
		pd.type = 0; // mono color
		pd.x = GLfloat(pt.x);
		pd.y = GLfloat(pt.y);
		pd.z = GLfloat(pt.z);
		pd.radius = pt_radius;
	}
	if (_pt_num > max_pcl_num_per_drawcall)
	{
		init_gl_buffer(pt_data, max_pcl_num_per_drawcall);
		pt_num = _pt_num;
	}
	else
	{
		init_gl_buffer(pt_data, _pt_num);
		delete[] pt_data;
		pt_data = nullptr;
	}
	return 0;
}

template <typename DivisionSet>
int QtMonoColorBallGLObject::init(
	const double* pcl_x,
	const double* pcl_y,
	const double* pcl_z,
	const double* pcl_vol,
	const size_t pcl_num,
	const QVector3D& c,
	const float radius_scale,
	const DivisionSet &div_set
	)
{
	clear();
	color = c;
	PointData* pt_data = new PointData[pcl_num];
	size_t cur_p_id = 0;
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		if (!div_set.inside(pcl_x[p_id], pcl_y[p_id], pcl_z[p_id]))
		{
			PointData& pd = pt_data[cur_p_id++];
			pd.type = 0; // mono color
			pd.x = GLfloat(pcl_x[p_id]);
			pd.y = GLfloat(pcl_y[p_id]);
			pd.z = GLfloat(pcl_z[p_id]);
			pd.radius = GLfloat(pow(3.0 * pcl_vol[p_id] / (4.0 * 3.14159265359), 1.0 / 3.0)) * radius_scale;
		}
	}
	if (pcl_num > max_pcl_num_per_drawcall)
	{
		init_gl_buffer(pt_data, max_pcl_num_per_drawcall);
		pt_num = pcl_num;
	}
	else
	{
		init_gl_buffer(pt_data, cur_p_id);
		delete[] pt_data;
		pt_data = nullptr;
	}
	return 0;
}

#endif