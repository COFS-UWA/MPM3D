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

	QtCharBitmapLoader char_loader;
	// cal font_wd (width), font_ht (height)
	void cal_str_geometry(const char *str, QtCharBitmapLoader& loader,
						  GLint &font_wd, GLint &font_ht);

	// get position of each component
	struct Rect
	{
		GLfloat xl, xu, yl, yu;
		inline GLfloat width() { return xu - xl; }
		inline GLfloat height() { return yu - yl; }
	};
	
	Rect cmap_rect;
	Rect content_rect;
	Rect title_rect;
	Rect cbar_rect;

	struct StrInfo
	{
		std::string str;
		GLfloat as; // aspect ratio
		GLuint font_pixel_wd;
		GLfloat nr; // normalized ratio
	};

	StrInfo title_str_info;
	Rect title_str_rect;

	struct EntryInfo
	{
		Rect cbox_rect;
		Rect num_rect;

		float number;
		StrInfo num_str_info;
		Rect num_str_rect;
	};
	std::vector<EntryInfo> entry_infos;

	void scale_rect_height(Rect& ori_rect, GLfloat sc, Rect& new_rect);

	void get_mid_align_str_pos(Rect &pos_rect, GLfloat str_as, Rect &str_rect);
	void get_left_align_str_pos(Rect& pos_rect, GLfloat str_as, Rect& str_rect);

	// form opengl buffer, single color
	struct Type0NodeData
	{
		GLint type;
		GLfloat x, y;
	};

	MemoryUtils::ItemArray<Type0NodeData> line_node_mem;
	MemoryUtils::ItemArray<GLuint> line_elem_mem;
	void add_line(GLfloat x0, GLfloat y0, GLfloat x1, GLfloat y1);

	// node data of type 1, multi-color
	struct Type1NodeData
	{
		GLint type;
		GLfloat x, y;
		GLfloat r, g, b;
	};
	
	MemoryUtils::ItemArray<Type1NodeData> tri_node_mem;
	MemoryUtils::ItemArray<GLuint> tri_elem_mem;
	void add_rect(Rect &rect, QVector3D &color);

	struct CharNodeData
	{
		GLfloat x, y;
		GLfloat tex_x, tex_y; // texture coordinates
	};
	size_t char_num;
	MemoryUtils::ItemArray<CharNodeData> char_node_mem;
	MemoryUtils::ItemArray<GLuint> char_textures;
	void add_char_texture(Rect &rect, GLuint tex_id);
	void add_string_texture(GLfloat xpos, GLfloat ypos, GLfloat nr, const char *str);

	GLfloat line_width;

	union
	{
		struct
		{
			GLuint line_vao;
			GLuint tri_vao;
			GLuint char_vao;
		};
		GLuint vaos[3];
	};
	union
	{
		struct
		{
			GLuint line_vbo, line_veo;
			GLuint tri_vbo, tri_veo;
			GLuint char_vbo;
		};
		GLuint vbos[5];
	};

public:
	QtUniformColorMapObject(QOpenGLFunctions_3_3_Core &_gl);
	~QtUniformColorMapObject();
	void clear_gl_buffer();

	int init(GLfloat _xpos, GLfloat _ypos, GLfloat ht,
		UniformColorMap& color_map,
		const char* title,
		const char* num_format,
		const char* ttf_filename);

	void draw(QOpenGLShaderProgram& shader_plain2D,
			  QOpenGLShaderProgram& shader_char);
};

#endif