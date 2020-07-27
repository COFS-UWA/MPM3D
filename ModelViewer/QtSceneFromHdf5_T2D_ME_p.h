#ifndef __Qt_Scene_From_Hdf5_T2D_ME_p_h__
#define __Qt_Scene_From_Hdf5_T2D_ME_p_h__

#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"
#include "ResultFile_hdf5.h"
#include "QtTriangleMeshGLObject.h"
#include "UniformColorMap_Abaqus.h"
#include "QtMultiColorCircleGLObject.h"
#include "QtRigidCircleObject.h"
#include "QtUniformColorMapObject.h"
#include "QtSceneFromHdf5.h"

class QtSceneFromHdf5_T2D_ME_p : public QtSceneFromHdf5
{
protected:
	struct PointData
	{
		GLfloat x, y;
		GLfloat va;
	};

	// hdf5 result file data infos
	ResultFile_hdf5* res_file;
	hid_t th_id;
	hid_t frame_grp_id;
	size_t pcl_num;
	hid_t pcl_dt_id;
	size_t pcl_size;
	size_t pcl_x_off;
	size_t pcl_y_off;
	size_t pcl_vol_off;
	std::string field_name;
	size_t pcl_fld_off;
	hid_t pcl_fld_type;
	MemoryUtils::ItemArray<char> pcls_data_mem;

	QOpenGLShaderProgram shader_plain2D;
	QOpenGLShaderProgram shader_circles;
	QOpenGLShaderProgram shader_char;

	bool display_bg_mesh;
	bool display_pcls;
	bool display_rc;

	QtTriangleMeshGLObject bg_mesh_obj;

	UniformColorMap_Abaqus color_map;
	GLuint color_map_texture;
	QtMultiColorCircleGLObject pcls_obj;
	
	bool has_rc_obj;
	QtRigidCircleObject rc_obj;
	
	bool has_color_map;
	float cm_xpos, cm_ypos, cm_ht;
	QtUniformColorMapObject color_map_obj;

	bool display_whole_model;
	GLfloat xl, xu, yl, yu, padding_ratio;

	// viewport
	GLint vp_x_pos, vp_y_pos;
	GLsizei vp_x_size, vp_y_size;
	// hud viewport
	GLint win_wd, win_ht;

	void set_viewport(int wd, int ht, GLfloat xlen, GLfloat ylen);
	
	QMatrix4x4 view_mat;
	QMatrix4x4 hud_view_mat;

	QVector3D bg_color;
	inline void set_bg_color(GLfloat r, GLfloat g, GLfloat b)
	{
		bg_color[0] = r;
		bg_color[1] = g;
		bg_color[2] = b;
	}

	void clear();

public:
	QtSceneFromHdf5_T2D_ME_p(QOpenGLFunctions_3_3_Core &_gl);
	~QtSceneFromHdf5_T2D_ME_p();
	void close_file();

	inline void set_display_bg_mesh(bool op = true) { display_bg_mesh = op; }
	inline void set_display_pcls(bool op = true) { display_pcls = op; }
	inline void set_display_rc(bool op = true) { display_rc = op; }
	
	inline void set_display_whole_model() { display_whole_model = true; }
	inline void set_display_range(double _xl, double _xu, double _yl, double _yu)
	{
		display_whole_model = false;
		xl = _xl; xu = _xu; yl = _yl; yu = _yu;
	}

	inline void set_fld_range(double min, double max)
	{ color_map.set_range(min, max); }
	inline void set_color_map_pos(float xpos, float ypos, float ht)
	{
		has_color_map = true;
		cm_xpos = xpos;	cm_ypos = ypos; cm_ht = ht;
	}

	size_t get_frame_num();
	double get_frame_time(size_t frame_id);
	
	int set_res_file(
		ResultFile_hdf5& rf,
		const char* th_na,
		const char* field_na
		);

public:
	// create the scene, including bg mesh and pcls
	int init_scene(int wd, int ht, size_t frame_id) override;
	// only update pcls, for animation
	void update_scene(size_t frame_id) override;
	void draw() override;
	void resize(int wd, int ht) override;
};

#endif