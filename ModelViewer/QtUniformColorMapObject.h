#ifndef __Qt_Uniform_Color_Map_Object_h__
#define __Qt_Uniform_Color_Map_Object_h__

#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "ItemArray.hpp"
#include "UniformColorMap.h"
#include "QtCharBitmapLoader.h"

class QtUniformColorMapObject
{
protected:
	QOpenGLFunctions_3_3_Core& gl;
	
	// geometry parameters
	GLfloat ht_margin_ratio, wd_margin_ratio;
	
	GLfloat gap_ratio; // gap between title and content
	GLfloat title_ratio; // ratio between title height and number

	GLfloat cbox_as;
	GLfloat cbox_tick_ratio;
	GLfloat cbox_gap_ratio;

	// geometry
	GLfloat xpos, ypos;
	GLfloat cmap_wd, cmap_ht;
	GLfloat content_wd, content_ht;
	GLfloat title_wd, title_ht;
	GLfloat title_gap_ht;
	GLfloat cbox_wd, cbox_ht;
	GLfloat cbox_tick_wd, cbox_gap_wd;
	GLfloat num_wd, num_ht;

	// strings in color map
	QtCharBitmapLoader char_loader;

	// position of each component
	Rect cmap_rect;
	Rect content_rect;
	Rect title_rect;
	Rect cbar_rect;
	std::vector<Rect> cboxes_rect;
	std::vector<Rect> nums_rect;

	struct StrInfo
	{
		std::string str;
		float number;
		GLuint bm_wd, bm_ht;
		GLfloat bm_as;
		GLfloat wd, ht;
		GLfloat norm_ratio; // normalized ratio
	};
	
	StrInfo title_str_info; // color map title
	std::vector<StrInfo> num_str_infos; // color map number
	
	void get_str_size(StrInfo& strif, QtCharBitmapLoader& loader);

	// mono-color node data
	struct MCNodeData
	{
		GLint type;
		GLfloat x, y;
	};

	struct NodeData
	{
		GLint type;
		GLfloat x, y;
		union
		{
			struct { GLfloat r, g, b; };
			struct { GLint tex_id; GLfloat tex_x, tex_y; };
		};
	};

	MemoryUtils::ItemArray<MCNodeData> line_node_mem;
	MemoryUtils::ItemArray<GLuint> line_elem_mem;
	void add_line(GLfloat x0, GLfloat y0, GLfloat x1, GLfloat y1);
	
	MemoryUtils::ItemArray<NodeData> tri_node_mem;
	MemoryUtils::ItemArray<GLuint> tri_elem_mem;
	void add_rect_color(Rect &rect, QVector3D &color);
	void add_char_texture(Rect &rect, GLint tex_id);
	void add_string_texture(GLfloat xpos, GLfloat ypos, StrInfo& strif);

	GLfloat line_width;
	GLuint line_vao, line_vbo, line_veo;
	GLuint tri_vao, tri_vbo, tri_veo;

public:
	QtUniformColorMapObject(QOpenGLFunctions_3_3_Core &_gl);
	~QtUniformColorMapObject();
	void clear_gl_buffer();

	int init(GLfloat _xpos, GLfloat _ypos, GLfloat ht,
		UniformColorMap& color_map,
		const char* title,
		const char* num_format,
		const char* ttf_filename,
		int win_wd, int win_ht);

	void draw(QOpenGLShaderProgram &shader);
};

#endif