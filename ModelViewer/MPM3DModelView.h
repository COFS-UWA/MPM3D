#ifndef __MPM_3D_Model_View_h__
#define __MPM_3D_Model_View_h__

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>
#include <QTimer>

#include "Geometry.h"
#include "TetrahedronMeshBuffer.h"
#include "MonoColorParticleBuffer.h"
#include "ValueToColor.h"
#include "MultiColorParticleBuffer.h"
#include "PointBuffer.h"

class MPM3DModelView : public QOpenGLWidget,
	public QOpenGLFunctions_3_3_Core
{
	Q_OBJECT
public:
	typedef ValueToColor::Colori Colori;
	typedef ValueToColor::Colorf Colorf;

protected:
	// camera info
	GLfloat fov_angle;
	QVector3D view_dir;
	QVector3D up_dir;
	// bounding circle
	float md_radius;
	QVector3D md_centre;
	float view_dist_scale;

	QMatrix4x4 view_mat;
	QMatrix4x4 proj_mat;
	void update_view_mat();
	void update_proj_mat();

	// shader
	QOpenGLShaderProgram shader_unicolor;
	QOpenGLShaderProgram shader_multicolor;

	// background mesh
	TetrahedronMeshBuffer bg_mesh_buf;
	bool need_to_paint_bg_mesh;

	// particle data
	MonoColorParticleBuffer pcl_buf;
	ValueToColor color_scale;
	MultiColorParticleBuffer cpcl_buf;
	bool is_monocolor_pcl;
	bool need_to_paint_pcl_buf;

	// point data
	PointBuffer point_buf;
	bool need_to_paint_point_buf;

public:
	// Controller of this window behavior
	class Controller
	{
	protected:
		MPM3DModelView *view;

	public:
		Controller(MPM3DModelView &v) : view(&v)
		{
			view->controller = this;
		}
		virtual ~Controller() {}

		inline void set_view(MPM3DModelView &v)
		{
			view = &v;
			view->controller = this;
		}

		// called in initializeGL()
		virtual int initialize_model_view_data() = 0;
		// callled in paintGL()
		virtual int before_render() { return 0; }
		virtual int after_render() { return 0; }
	};

protected:
	Controller *controller;

public:
	explicit MPM3DModelView(QWidget* parent = Q_NULLPTR);
	~MPM3DModelView();
	void initializeGL();
	void paintGL();
	void resizeGL(int width, int height);

	inline ValueToColor &get_color_scale() { return color_scale; }
	
	inline void set_display_bg_mesh(bool op = true) { need_to_paint_bg_mesh = op; }
	inline void set_display_pcls(bool op = true) { need_to_paint_pcl_buf = op; }
	inline void set_display_points(bool op = true) { need_to_paint_point_buf = op; }

	inline void set_view_dir(QVector3D &_view_dir) { view_dir = _view_dir; }
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
	
	// ==================== init model data buffer =================
	template <typename TMesh>
	inline int init_bg_mesh(TMesh &mesh, QVector3D &_color)
	{
		// init mesh_range
		Cube mh_bbox = mesh.get_bounding_box();
		md_centre.setX(float(mh_bbox.xl + mh_bbox.xu) * 0.5f);
		md_centre.setY(float(mh_bbox.yl + mh_bbox.yu) * 0.5f);
		md_centre.setZ(float(mh_bbox.zl + mh_bbox.zu) * 0.5f);
		float dx = mh_bbox.xu - mh_bbox.xl;
		float dy = mh_bbox.yu - mh_bbox.yl;
		float dz = mh_bbox.zu - mh_bbox.zl;
		md_radius = sqrt(dx * dx + dy * dy + dz * dz) * 0.5f;
		
		return bg_mesh_buf.init(mesh, _color);
	}

	template <typename Node, typename Element>
	inline int init_bg_mesh(
		Node *nodes, size_t node_num,
		Element *elems, size_t elem_num,
		QVector3D& _color
		)
	{
		if (!nodes || node_num == 0 ||
			!elems || elem_num == 0)
			return 0;

		// init mesh_range
		Cube mh_bbox;
		mh_bbox.xl = nodes[0].x;
		mh_bbox.xu = mh_bbox.xl;
		mh_bbox.yl = nodes[0].y;
		mh_bbox.yu = mh_bbox.yl;
		mh_bbox.zl = nodes[0].z;
		mh_bbox.zu = mh_bbox.zl;
		for (size_t n_id = 1; n_id < node_num; ++n_id)
		{
			Node& n = nodes[n_id];
			if (mh_bbox.xl > n.x)
				mh_bbox.xl = n.x;
			if (mh_bbox.xu < n.x)
				mh_bbox.xu = n.x;
			if (mh_bbox.yl > n.y)
				mh_bbox.yl = n.y;
			if (mh_bbox.yu < n.y)
				mh_bbox.yu = n.y;
			if (mh_bbox.zl > n.z)
				mh_bbox.zl = n.z;
			if (mh_bbox.zu < n.z)
				mh_bbox.zu = n.z;
		}
		md_centre.setX(float(mh_bbox.xl + mh_bbox.xu) * 0.5f);
		md_centre.setY(float(mh_bbox.yl + mh_bbox.yu) * 0.5f);
		md_centre.setZ(float(mh_bbox.zl + mh_bbox.zu) * 0.5f);
		float dx = mh_bbox.xu - mh_bbox.xl;
		float dy = mh_bbox.yu - mh_bbox.yl;
		float dz = mh_bbox.zu - mh_bbox.zl;
		md_radius = sqrt(dx * dx + dy * dy + dz * dz) * 0.5f;

		return bg_mesh_buf.init<Node, Element>(nodes, node_num, elems, elem_num, _color);
	}

	template <typename Particle>
	inline int init_monocolor_pcl_data(Particle *pcls, size_t pcl_num, QVector3D& _color)
	{
		is_monocolor_pcl = true;
		return pcl_buf.init_pcl_data<Particle>(pcls, pcl_num, _color);
	}

	template <typename Particle>
	inline int update_monocolor_pcl_data(Particle *pcls)
	{
		is_monocolor_pcl = true;
		return pcl_buf.update_pcl_data<Particle>(pcls);
	}

	inline int init_color_scale(
		double lower, double upper,
		Colori *colors, size_t color_num,
		bool out_of_bound_color = true
		)
	{
		return color_scale.init(lower, upper, 
			colors, color_num, out_of_bound_color);
	}
	
	template <typename FieldType = double>
	inline int init_multicolor_pcl_data(
		char* pcls_data, size_t pcl_size, size_t pcl_num,
		size_t x_off, size_t y_off, size_t z_off,
		size_t vol_off, float vol_scale,
		size_t fld_off, ValueToColor& color_scale
		)
	{
		is_monocolor_pcl = false;
		return cpcl_buf.init_data<FieldType>(
							pcls_data, pcl_size, pcl_num,
							x_off, y_off, z_off,
							vol_off, vol_scale,
							fld_off, color_scale);
	}

	template <typename FieldType = double>
	inline int update_multicolor_pcl_data(char* pcls)
	{
		is_monocolor_pcl = false;
		return cpcl_buf.update_data<FieldType>(pcls);
	}

	inline int init_point_data(Point3D* points, size_t point_num,
		GLfloat point_size,	QVector3D& points_color)
	{
		return point_buf.init(points, point_num, point_size, points_color);
	}

public: // whether window is fully loaded
	// display model and take screenshot after window is fully loaded
	inline bool is_fully_loaded() { return win_is_fully_loaded; }
private:
	bool win_is_fully_loaded;
	QTimer init_timer;
private slots:
	void init_timer_func();
};

#endif