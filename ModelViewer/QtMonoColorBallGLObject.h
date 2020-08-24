#ifndef __Qt_Mono_Color_Ball_GL_Object_h__
#define __Qt_Mono_Color_Ball_GL_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

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

	void clear();

	int init_ball_data();
	int init_gl_buffer(PointData* pds, size_t pd_num);

public:
	QtMonoColorBallGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtMonoColorBallGLObject();

	// Particle3D has member x, y, z and area
	template <typename Particle3D>
	int init(Particle3D* pcls, size_t pcl_num, QVector3D& c, float radius_scale = 0.5f);

	// Point3D has member x, y and z
	template <typename Point3D>
	int init(Point3D* pts, size_t pt_num, float pcl_radius, QVector3D& c);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename Particle3D>
int QtMonoColorBallGLObject::init(
	Particle3D *pcls,
	size_t pcl_num,
	QVector3D& c,
	float radius_scale
	)
{
	clear();
	color = c;
	PointData* pt_data = new PointData[pcl_num];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle3D &pcl = pcls[p_id];
		PointData& pd = pt_data[p_id];
		pd.type = 0; // mono color
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
		pd.z = GLfloat(pcl.z);
		pd.radius = GLfloat(pow(3.0 * pcl.get_vol()/ (4.0 * 3.14159265359), 0.3333333))	* radius_scale;
	}
	int res = init_gl_buffer(pt_data, pcl_num);
	delete[] pt_data;
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
	PointData* pt_data = new PointData[_pt_num];
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
	int res = init_gl_buffer(pt_data, _pt_num);
	delete[] pt_data;
	return res;
}

#endif