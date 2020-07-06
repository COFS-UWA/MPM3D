#ifndef __Qt_Mono_Colour_Circle_GL_Object_h__
#define __Qt_Mono_Colour_Circle_GL_Object_h__

#include <QColor>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

class QtMonoColourCircleGLObject
{
protected:
	QOpenGLFunctions_3_3_Core& gl;

	struct PointData
	{
		GLint type;
		GLfloat x, y;
		GLfloat radius;
	};

	GLuint vao, vbo_cs, veo_cs, vbo_pts;
	size_t c_elem_node_num;
	size_t pt_num;
	QVector3D color;

	void clear();

	int init_circle_data();
	int init_gl_buffer(PointData* pds, size_t pd_num);

public:
	QtMonoColourCircleGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtMonoColourCircleGLObject();

	// Particle2D has member x, y and area
	template <typename Particle2D>
	int init(Particle2D* pcls, size_t pcl_num, QVector3D &c, float radius_scale = 0.5f);

	// Point2D has member x and y
	template <typename Point2D>
	int init(Point2D *pts, size_t pt_num, float pcl_radius, QVector3D &c);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename Particle2D>
int QtMonoColourCircleGLObject::init(
	Particle2D* pcls,
	size_t pcl_num,
	QVector3D &c,
	float radius_scale
	)
{
	clear();

	color = c;

	PointData* pt_datas = new PointData[pcl_num];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle2D& pcl = pcls[p_id];
		PointData& pd = pt_datas[p_id];
		pd.type = 0; // mono color
		pd.x = GLfloat(pcl.x);
		pd.y = GLfloat(pcl.y);
		pd.radius = sqrt(GLfloat(pcl.get_vol()) / 3.14159265359f)
						* radius_scale;
	}
	int res = init_gl_buffer(pt_datas, pcl_num);
	delete[] pt_datas;
	return 0;
}

template <typename Point2D>
int QtMonoColourCircleGLObject::init(
	Point2D* pts,
	size_t _pt_num,
	float pt_radius,
	QVector3D& c
)
{
	clear();

	color = c;

	PointData* pt_datas = new PointData[_pt_num];
	for (size_t p_id = 0; p_id < _pt_num; ++p_id)
	{
		Point2D& pt = pts[p_id];
		PointData& pd = pt_datas[p_id];
		pd.type = 0; // mono color
		pd.x = GLfloat(pt.x);
		pd.y = GLfloat(pt.y);
		pd.radius = pt_radius;
	}
	int res = init_gl_buffer(pt_datas, _pt_num);
	delete[] pt_datas;
	return res;
}

#endif