#ifndef __Qt_Scene_From_Model_T3D_ME_mt_h__
#define __Qt_Scene_From_Model_T3D_ME_mt_h__

#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"
#include "Model_T3D_ME_mt.h"
#include "QtTetrahedronMeshGLObject.h"
#include "QtMonoColorBallGLObject.h"
#include "QtRigidCylinderObject.h"
#include "QtSceneFromModel.h"

class QtSceneFromModel_T3D_ME_mt : public QtSceneFromModel
{
protected:
	struct PointData { GLfloat x, y, z; };

	Model_T3D_ME_mt *model;

	size_t pt_num;
	PointData* pts;
	float pt_radius;
	MemoryUtils::ItemArray<PointData> pts_mem;

	bool display_bg_mesh;
	bool display_pcls;
	bool display_pts;
	bool display_rcy, has_rcy;

	QtTetrahedronMeshGLObject bg_mesh_obj;
	QtMonoColorBallGLObject pcls_obj;
	QtMonoColorBallGLObject pts_obj;
	QtRigidCylinderObject rcy_obj;

	// camera info
	GLfloat fov_angle;
	QVector3D view_dir;
	QVector3D up_dir;
	QVector3D view_pos;
	// bounding circle
	float md_radius;
	QVector3D md_centre;
	float view_dist_scale;

	int width, height;

	QMatrix4x4 view_mat;
	QMatrix4x4 proj_mat;
	void update_view_mat();
	void update_proj_mat();

	QVector3D bg_color;

	QVector3D light_dir;
	float light_dist_scale;
	QVector3D light_pos;
	void update_light_pos();

	float fog_coef;
	QVector3D fog_color;

	QVector3D light_color;
	float amb_coef;
	float diff_coef;
	float spec_coef;
	float spec_shininess;

	QOpenGLShaderProgram shader_plain3D;
	QOpenGLShaderProgram shader_balls;
	QOpenGLShaderProgram shader_rigid_mesh;
	
public:
	typedef Model_T3D_ME_mt::Position Position;
	typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;

	explicit QtSceneFromModel_T3D_ME_mt(QOpenGLFunctions_3_3_Core &_gl);
	~QtSceneFromModel_T3D_ME_mt();

	inline void set_view_dir(QVector3D& _dir) { view_dir = _dir; }
	inline void set_view_dir(float x, float y, float z)
	{
		view_dir[0] = x; view_dir[1] = y; view_dir[2] = z;
	}
	// theta and fai are in degree, 0 < theta < 360, -90 < fai < 90
	inline void set_view_dir(float theta, float fai)
	{
		theta *= 3.14159265359f / 180.0f;
		fai *= 3.14159265359f / 180.0f;
		double cos_fai = cos(fai);
		double vd_x = cos_fai * cos(theta);
		double vd_y = cos_fai * sin(theta);
		double vd_z = sin(fai);
		set_view_dir(-vd_x, -vd_y, -vd_z);
	}
	inline void set_view_dist_scale(float scale) { view_dist_scale = scale; }

	inline void set_light_dir(QVector3D& _dir) { light_dir = _dir; }
	inline void set_light_dir(float x, float y, float z)
	{
		light_dir[0] = x; light_dir[1] = y; light_dir[2] = z;
	}
	// theta and fai are in degree, 0 < theta < 360, -90 < fai < 90
	inline void set_light_dir(float theta, float fai)
	{
		theta *= 3.14159265359f / 180.0f;
		fai *= 3.14159265359f / 180.0f;
		double cos_fai = cos(fai);
		double vd_x = cos_fai * cos(theta);
		double vd_y = cos_fai * sin(theta);
		double vd_z = sin(fai);
		set_light_dir(-vd_x, -vd_y, -vd_z);
	}
	inline void set_light_dist_scale(float scale) { light_dist_scale = scale; }

	inline void set_bg_color(QVector3D& color) { bg_color = color; }
	inline void set_bg_color(GLfloat r, GLfloat g, GLfloat b)
	{
		bg_color[0] = r; bg_color[1] = g; bg_color[2] = b;
	}

	inline void set_fog_coef(float coef) { fog_coef = coef; }
	inline void set_fog_color(QVector3D& color) { fog_color = color; }

	inline void set_light_color(QVector3D& color) { light_color = color; }
	inline void set_amb_coef(float coef) { amb_coef = coef; }
	inline void set_diff_coef(float coef) { diff_coef = coef; }
	inline void set_spec_coef(float coef) { spec_coef = coef; }
	inline void set_spec_shininess(float shininess) { spec_shininess = shininess; }

	// ============================= model data ============================
	inline void set_display_bg_mesh(bool op = true) { display_bg_mesh = op; }
	inline void set_display_pcls(bool op = true) { display_pcls = op; }
	inline void set_display_pts(bool op = true) { display_pts = op; }

	inline void set_model(Model_T3D_ME_mt &_model) { model = &_model; }
	int set_pts_from_pcl_id(size_t* ids, size_t id_num, float radius);
	int set_pts_from_node_id(size_t* ids, size_t id_num, float radius);
	template<typename Point3D>
	int set_pts(Point3D* _pts, size_t _pt_num, float radius)
	{
		pt_radius = radius;
		pt_num = _pt_num;
		pts_mem.reserve(pt_num);
		pts = pts_mem.get_mem();
		for (size_t p_id = 0; p_id < pt_num; ++p_id)
		{
			Point3D& p = _pts[p_id];
			PointData& pd = pts[p_id];
			pd.x = GLfloat(p.x);
			pd.y = GLfloat(p.y);
			pd.z = GLfloat(p.z);
		}
		return 0;
	}

	int initialize(int wd, int ht);
	void draw();
	void resize(int wd, int ht);
};

#endif