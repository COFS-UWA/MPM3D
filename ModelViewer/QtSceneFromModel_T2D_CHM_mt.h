#ifndef __Qt_Scene_From_Model_T2D_CHM_mt_h__
#define __Qt_Scene_From_Model_T2D_CHM_mt_h__

#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"
#include "Model_T2D_CHM_mt.h"
#include "QtTriangleMeshGLObject.h"
#include "QtRigidCircleObject.h"
#include "QtMonoColorCircleGLObject.h"
#include "QtSceneFromModel.h"

class QtSceneFromModel_T2D_CHM_mt :
	public QtSceneFromModel
{
protected:
	struct PointData { GLfloat x, y; };

	Model_T2D_CHM_mt* model;
	size_t pt_num;
	PointData* pts;
	float pt_radius;
	MemoryUtils::ItemArray<PointData> pts_mem;

	QOpenGLShaderProgram shader_plain2D;
	QOpenGLShaderProgram shader_circles;

	bool display_bg_mesh;
	bool display_pcls;
	bool display_rigid_circle;
	bool display_pts;

	QtTriangleMeshGLObject bg_mesh_obj;
	QtMonoColorCircleGLObject pcls_obj;
	QtRigidCircleObject rc_obj;
	// points to be high lighted
	QtMonoColorCircleGLObject pts_obj;

	bool display_whole_model;
	GLfloat padding_ratio;
	Rect display_bbox;

	// viewport info
	GLint vp_x_pos, vp_y_pos;
	GLsizei vp_x_size, vp_y_size;

	QMatrix4x4 view_mat;

	void set_viewport(int wd, int ht, GLfloat xlen, GLfloat ylen);

	QVector3D bg_color;
	inline void set_bg_color(GLfloat r, GLfloat g, GLfloat b)
	{ bg_color[0] = r; bg_color[1] = g; bg_color[2] = b; }

public:
	typedef Model_T2D_CHM_mt::ElemNodeIndex ElemNodeIndex;

	QtSceneFromModel_T2D_CHM_mt(QOpenGLFunctions_3_3_Core &_gl);
	~QtSceneFromModel_T2D_CHM_mt();

	inline void set_display_bg_mesh(bool op = true) { display_bg_mesh = op; }
	inline void set_display_pcls(bool op = true) { display_pcls = op; }
	inline void set_display_rc(bool op = true) { display_rigid_circle = op; }
	inline void set_display_pts(bool op = true) { display_pts = op; }

	inline void set_display_whole_model() { display_whole_model = true; }
	inline void set_display_range(double _xl, double _xu, double _yl, double _yu)
	{
		display_whole_model = false;
		display_bbox.xl = _xl; display_bbox.xu = _xu;
		display_bbox.yl = _yl; display_bbox.yu = _yu;
	}

	inline void set_model(Model_T2D_CHM_mt& _model) { model = &_model; }
	int set_pts_from_pcl_id(size_t* ids, size_t id_num, float radius);
	int set_pts_from_node_id(size_t* ids, size_t id_num, float radius);
	template<typename Point2D>
	int set_pts(Point2D* _pts, size_t _pt_num, float radius)
	{
		pt_num = _pt_num;
		pt_radius = radius;
		pts_mem.reserve(pt_num);
		pts = pts_mem.get_mem();
		for (size_t p_id = 0; p_id < pt_num; ++p_id)
		{
			Point2D& p = _pts[p_id];
			PointData& pd = pts[p_id];
			pd.x = GLfloat(p.x);
			pd.y = GLfloat(p.y);
		}
		return 0;
	}
	int set_pts_from_vx_s_bc(float radius);
	int set_pts_from_vy_s_bc(float radius);
	int set_pts_from_vx_f_bc(float radius);
	int set_pts_from_vy_f_bc(float radius);

	int initialize(int wd, int ht);
	void draw();
	void resize(int wd, int ht);
};

#endif