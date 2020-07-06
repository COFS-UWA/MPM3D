#ifndef __Qt_Multi_Colour_Circle_GL_Object_h__
#define __Qt_Multi_Colour_Circle_GL_Object_h__

#include <QColor>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

class QtMultiColourCircleGLObject
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

	GLuint vao, vbo_cs, veo_cs, vbo_pts;
	size_t c_elem_node_num;
	size_t pt_num;

	void clear();

	int init_circle_data();
	int init_gl_buffer(PointData *pds, size_t pd_num);

public:
	QtMultiColourCircleGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtMultiColourCircleGLObject();

	template <typename FieldDataType>
	int init(
		char* pcls_data, size_t pcl_size, size_t pcl_num,
		size_t x_off, size_t y_off,
		size_t vol_off, float vol_scale,
		size_t fld_off);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename FieldDataType>
int QtMultiColourCircleGLObject::init(
	char* pcls_data, size_t pcl_size, size_t pcl_num,
	size_t x_off, size_t y_off,
	size_t vol_off, float radius_scale,
	size_t fld_off)
{
	GLfloat pcl_vol;
	char* pcl_data = pcls_data;
	PointData* pt_datas = new PointData[pcl_num];
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		PointData& pd = pt_datas[pcl_id];
		pd.type = 1; // multi color
		pd.x = GLfloat(*(double *)(pcl_data + x_off));
		pd.y = GLfloat(*(double *)(pcl_data + y_off));
		pcl_vol = GLfloat(*(double *)(pcl_data + vol_off));
		pd.radius = sqrt(pcl_vol / 3.14159265359f) * radius_scale;
		pd.fld = GLfloat(*(FieldDataType*)(pcl_data + fld_off));
		pcl_data += pcl_size;
	}
	int res = init_gl_buffer(pt_datas, pcl_num);
	delete[] pt_datas;
	return res;
}

#endif