#ifndef __MPM_3D_Model_View_h__
#define __MPM_3D_Model_View_h__

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

#include "Geometry.h"
#include "TetrahedronMeshBuffer.h"
#include "MonoColorParticleBuffer.h"
#include "PointBuffer.h"

class TimeHistory_ModelView;

class MPM3DModelView : public QOpenGLWidget,
	public QOpenGLFunctions_3_3_Core
{
	Q_OBJECT
protected:
	// shader
	QOpenGLShaderProgram shader_unicolor;

	QMatrix4x4 view_mat;
	QMatrix4x4 proj_mat;

	// camera info
	GLfloat fov_angle;
	QVector3D view_dir;
	QVector3D up_dir;

	// background mesh
	TetrahedronMeshBuffer bg_mesh_buf;
	bool need_to_paint_bg_mesh;
	// bounding circle
	float mh_radius;
	QVector3D mh_centre;
	float view_dist_scale;

	// particle data
	MonoColorParticleBuffer pcl_buf;
	bool need_to_paint_pcl_buf;

	// point data
	PointBuffer point_buf;
	bool need_to_paint_point_buf;

	// model view timehistory
	friend TimeHistory_ModelView;
	TimeHistory_ModelView *th_mv;

	void update_view_mat();
	void update_proj_mat();

public:
	explicit MPM3DModelView(QWidget* parent = Q_NULLPTR);
	~MPM3DModelView();
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);

	inline void set_display_bg_mesh(bool op = true) { need_to_paint_bg_mesh = op; }
	inline void set_display_pcls(bool op = true) { need_to_paint_pcl_buf = op; }
	inline void set_display_points(bool op = true) { need_to_paint_point_buf = op; }

	inline void set_view_dir(QVector3D& _view_dir) { view_dir = _view_dir; }
	inline void set_view_dir(float x, float y, float z)
	{
		view_dir.setX(x);
		view_dir.setY(y);
		view_dir.setZ(z);
	}
	// theta and fai are in degree, 0 < theta < 360, -90 < fai < 90
	inline void set_view_dir(float theta, float fai)
	{
		theta *= 3.14159265359f / 180.0f;
		fai *= 3.14159265359f / 180.0f;
		double cos_fai = cos(fai);
		double sin_fai = sin(fai);
		double vd_x = cos_fai * cos(theta);
		double vd_y = cos_fai * sin(theta);
		double vd_z = sin_fai;
		set_view_dir(QVector3D(-vd_x, -vd_y, -vd_z));
	}
	inline void set_view_dist_scale(float scale) { view_dist_scale = scale; }
	
	template <typename TMesh>
	inline int init_bg_mesh(TMesh &mesh, QVector3D&_color)
	{
		// init mesh_range
		Cube mh_bbox = mesh.get_bounding_box();
		mh_centre.setX(float(mh_bbox.xl + mh_bbox.xu) * 0.5f);
		mh_centre.setY(float(mh_bbox.yl + mh_bbox.yu) * 0.5f);
		mh_centre.setZ(float(mh_bbox.zl + mh_bbox.zu) * 0.5f);
		float dx = mh_bbox.xu - mh_bbox.xl;
		float dy = mh_bbox.yu - mh_bbox.yl;
		float dz = mh_bbox.zu - mh_bbox.zl;
		mh_radius = sqrt(dx * dx + dy * dy + dz * dz) * 0.5f;
		
		return bg_mesh_buf.init(mesh, _color);
	}
	
	template <typename Particle>
	inline int init_monocolor_pcl_data(Particle* pcls, size_t pcl_num, QVector3D& _color)
	{
		return pcl_buf.init_pcl_data(pcls, pcl_num, _color);
	}

	template <typename Particle>
	inline int update_monocolor_pcl_data(Particle *pcls)
	{
		return pcl_buf.update_pcl_data<Particle>(pcls);
	}

	inline int init_point_data(Point3D* points, size_t point_num,
		GLfloat point_size,	QVector3D& points_color)
	{
		return point_buf.init(points, point_num, point_size, points_color);
	}
};

#endif