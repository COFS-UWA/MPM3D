#ifndef __Qt_Scene_From_Hdf5_T2D_CHM_mt_h__
#define __Qt_Scene_From_Hdf5_T2D_CHM_mt_h__

#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"
#include "Geometry2D.h"
#include "ResultFile_hdf5.h"
#include "Hdf5Field.h"
#include "QtTriangleMeshGLObject.h"
#include "UniformColorMap_Abaqus.h"
#include "QtMultiColorCircleGLObject.h"
#include "QtRigidCircleObject.h"
#include "QtUniformColorMapObject.h"
#include "QtSceneFromHdf5.h"

class QtSceneFromHdf5_T2D_CHM_mt :
	public QtSceneFromHdf5
{
protected:
	struct PointData { GLfloat x, y, va; };

	// hdf5 result file data infos
	ResultFile_hdf5* res_file;
	Hdf5DataLoader data_loader;
	Hdf5FieldExtraction_x x_fld;
	Hdf5FieldExtraction_y y_fld;
	Hdf5FieldExtraction_vol_m vol_fld;
	Hdf5FieldExtraction* pfld;
	MemoryUtils::ItemArray<double> pcl_fld_mem;

	std::string field_name;
	hid_t th_id;
	
	QOpenGLShaderProgram shader_plain2D;
	QOpenGLShaderProgram shader_circles;
	QOpenGLShaderProgram shader_char;

	bool display_bg_mesh;
	bool display_pcls;
	bool display_rc;
	bool has_rc_obj;

	QtTriangleMeshGLObject bg_mesh_obj;

	UniformColorMap_Abaqus color_map;
	GLuint color_map_texture;
	QtMultiColorCircleGLObject pcls_obj;
	
	QtRigidCircleObject rc_obj;
	
	bool has_color_map;
	float cm_xpos, cm_ypos, cm_ht;
	QtUniformColorMapObject color_map_obj;

	bool display_whole_model;
	Rect bbox;
	GLfloat padding_ratio;

	// viewport
	GLint vp_x_pos, vp_y_pos;
	GLsizei vp_x_size, vp_y_size;
	// hud viewport
	GLint win_wd, win_ht;

	void set_viewport(int wd, int ht, GLfloat xlen, GLfloat ylen);
	
	QMatrix4x4 view_mat;
	QMatrix4x4 hud_view_mat;

	QVector3D bg_color;

	void clear();

public:
	explicit QtSceneFromHdf5_T2D_CHM_mt(QOpenGLFunctions_3_3_Core &_gl);
	~QtSceneFromHdf5_T2D_CHM_mt();
	void close_file();

	inline void set_bg_color(GLfloat r, GLfloat g, GLfloat b)
	{ bg_color[0] = r; bg_color[1] = g; bg_color[2] = b; }
	inline void set_char_color(GLfloat r, GLfloat g, GLfloat b)
	{ color_map_obj.set_char_color(r, g, b); }
	inline void set_display_bg_mesh(bool op = true) { display_bg_mesh = op; }
	inline void set_display_pcls(bool op = true) { display_pcls = op; }
	inline void set_display_rc(bool op = true) { display_rc = op; }
	
	inline void set_display_whole_model() { display_whole_model = true; }
	inline void set_display_range(double _xl, double _xu, double _yl, double _yu)
	{
		display_whole_model = false;
		bbox.xl = _xl; bbox.xu = _xu;
		bbox.yl = _yl; bbox.yu = _yu;
	}

	inline void set_color_map_fld_range(double min, double max)
	{ color_map.set_range(min, max); }
	inline void set_color_map_geometry(float xpos, float ypos, float ht)
	{
		has_color_map = true;
		cm_xpos = xpos;	cm_ypos = ypos; cm_ht = ht;
	}

	size_t get_frame_num();
	double get_frame_time(size_t frame_id);
	
	int set_res_file(
		ResultFile_hdf5& rf,
		const char* th_name,
		Hdf5Field::FieldType fld_type
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