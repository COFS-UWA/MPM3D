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

	gl.glGenVertexArrays(1, &vao);
	gl.glBindVertexArray(vao);
	
	int res = init_circle_data();
	if (res) return res;

	color = c;
	pt_num = pcl_num;
	
	gl.glGenBuffers(1, &vbo_cs);
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_cs);
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
	gl.glBufferData(GL_ARRAY_BUFFER,
		pcl_num * sizeof(PointData),
		pt_datas,
		GL_STREAM_DRAW
		);
	delete[] pt_datas;

	// pt_type
	gl.glVertexAttribPointer(1,
		1, GL_INT, GL_FALSE,
		sizeof(PointData),
		(GLvoid *)offsetof(PointData, type)
		);
	gl.glEnableVertexAttribArray(1);
	gl.glVertexAttribDivisor(1, 1);
	// pt_pos
	gl.glVertexAttribPointer(2,
		2, GL_FLOAT, GL_FALSE,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, x)
		);
	gl.glEnableVertexAttribArray(2);
	gl.glVertexAttribDivisor(2, 1);
	// pt_radius
	gl.glVertexAttribPointer(3,
		1, GL_FLOAT, GL_FALSE,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, radius)
		);
	gl.glEnableVertexAttribArray(3);
	gl.glVertexAttribDivisor(3, 1);
	// point value (not used)
	gl.glVertexAttribPointer(4,
		1, GL_FLOAT, GL_FALSE,
		0, (GLvoid*)0
		);
	gl.glEnableVertexAttribArray(4);
	gl.glVertexAttribDivisor(4, 1);

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

	pt_num = _pt_num;
	color = c;

	gl.glGenVertexArrays(1, &vao);
	gl.glBindVertexArray(vao);

	int res = init_circle_data();
	if (res)
		return res;

	gl.glGenBuffers(1, &vbo_cs);
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_cs);
	PointData* pt_datas = new PointData[pt_num];
	for (size_t p_id = 0; p_id < pt_num; ++p_id)
	{
		Point2D& pt = pts[p_id];
		PointData& pd = pt_datas[p_id];
		pd.type = 0; // mono color
		pd.x = GLfloat(pt.x);
		pd.y = GLfloat(pt.y);
		pd.radius = pt_radius;
	}
	gl.glBufferData(GL_ARRAY_BUFFER,
		pt_num * sizeof(PointData),
		pt_datas,
		GL_STREAM_DRAW
	);
	delete[] pt_datas;

	// pt_type
	gl.glVertexAttribPointer(1,
		1, GL_INT, GL_FALSE,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, type)
	);
	gl.glEnableVertexAttribArray(1);
	gl.glVertexAttribDivisor(1, 1);
	// pt_pos
	gl.glVertexAttribPointer(2,
		2, GL_FLOAT, GL_FALSE,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, x)
	);
	gl.glEnableVertexAttribArray(2);
	gl.glVertexAttribDivisor(2, 1);
	// pt_radius
	gl.glVertexAttribPointer(3,
		1, GL_FLOAT, GL_FALSE,
		sizeof(PointData),
		(GLvoid*)offsetof(PointData, radius)
	);
	gl.glEnableVertexAttribArray(3);
	gl.glVertexAttribDivisor(3, 1);
	// point value (not used)
	gl.glVertexAttribPointer(4,
		1, GL_FLOAT, GL_FALSE,
		0, (GLvoid*)0
	);
	gl.glEnableVertexAttribArray(4);
	gl.glVertexAttribDivisor(4, 1);

	return 0;
}

#endif