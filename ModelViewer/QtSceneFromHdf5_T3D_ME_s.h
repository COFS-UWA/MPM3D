#ifndef __Qt_Scene_From_Hdf5_T3D_ME_s_h__
#define __Qt_Scene_From_Hdf5_T3D_ME_s_h__

#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"
#include "ResultFile_hdf5.h"
#include "Hdf5Field.h"
#include "QtTetrahedronMeshGLObject.h"
#include "QtMultiColorBallGLObject.h"
#include "UniformColorMap_Abaqus.h"
#include "QtUniformColorMapObject.h"
#include "QtSceneFromHdf5.h"

class QtSceneFromHdf5_T3D_ME_s : public QtSceneFromHdf5
{
protected:
	struct PointData
	{
		GLfloat x, y, z;
		GLfloat va;
	};

	// hdf5 result file data infos
	ResultFile_hdf5* res_file;
	Hdf5DataLoader data_loader;
	Hdf5FieldExtraction_x x_fld;
	Hdf5FieldExtraction_y y_fld;
	Hdf5FieldExtraction_z z_fld;
	Hdf5FieldExtraction_vol vol_fld;
	Hdf5FieldExtraction *pfld;
	MemoryUtils::ItemArray<double> pcl_x_fld_mem;
	MemoryUtils::ItemArray<double> pcl_y_fld_mem;
	MemoryUtils::ItemArray<double> pcl_z_fld_mem;
	MemoryUtils::ItemArray<double> pcl_fld_mem;

	std::string field_name;
	hid_t th_id;

	bool display_bg_mesh;
	bool display_pcls;

	QtTetrahedronMeshGLObject bg_mesh_obj;

	UniformColorMap_Abaqus color_map;
	GLuint color_map_texture;
	QtMultiColorBallGLObject pcls_obj;
	
	bool has_color_map;
	float cm_xpos, cm_ypos, cm_ht;
	QtUniformColorMapObject color_map_obj;

	// camera info
	GLfloat fov_angle;
	QVector3D view_dir;
	QVector3D up_dir;
	QVector3D view_pos;
	// bounding circle
	float md_radius;
	QVector3D md_centre;
	float view_dist_scale;

	int win_wd, win_ht;
	
	QMatrix4x4 view_mat;
	QMatrix4x4 hud_view_mat;
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
	QOpenGLShaderProgram shader_plain2D;
	QOpenGLShaderProgram shader_char;
	
	void clear();

public:
	explicit QtSceneFromHdf5_T3D_ME_s(QOpenGLFunctions_3_3_Core &_gl);
	~QtSceneFromHdf5_T3D_ME_s();
	void close_file();

	inline void set_view_dir(QVector3D& _dir) { view_dir = _dir; }
	inline void set_view_dir(float x, float y, float z)
	{ view_dir[0] = x; view_dir[1] = y; view_dir[2] = z; }
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
	{ light_dir[0] = x; light_dir[1] = y; light_dir[2] = z; }
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
	
	inline void set_bg_color(QVector3D &color)
	{ bg_color = color; }
	inline void set_bg_color(GLfloat r, GLfloat g, GLfloat b)
	{ bg_color[0] = r; bg_color[1] = g; bg_color[2] = b; }

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
	
	inline void set_color_map_fld_range(double min, double max)
	{ color_map.set_range(min, max); }
	inline void set_color_map_geometry(float xpos, float ypos, float ht)
	{
		has_color_map = true;
		cm_xpos = xpos;	cm_ypos = ypos; cm_ht = ht;
	}

	size_t get_frame_num();
	double get_frame_time(size_t frame_id);
	
	int set_res_file(ResultFile_hdf5& rf,
		const char* th_name, Hdf5Field::FieldType fld_type);

	// create the scene, including bg mesh and pcls
	int init_scene(int wd, int ht, size_t frame_id) override;
	// only update pcls, for animation
	void update_scene(size_t frame_id) override;
	void draw() override;
	void resize(int wd, int ht) override;
};

#endif