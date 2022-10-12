#include "ModelViewer_pcp.h"

#include <iostream>

#include <ft2build.h>
#include FT_FREETYPE_H

#include "QtUniformcolorMapObject.h"

// "empirical factor"
#define Title_Font_Size_Scale 1.0f // 0.7f
#define Num_Font_Size_Scale 0.55f

QtUniformColorMapObject::QtUniformColorMapObject(
	QOpenGLFunctions_3_3_Core& _gl) : gl(_gl), char_loader(_gl),
	char_color(1.0f, 1.0f, 1.0f), 
	// geometry parameters
	ht_margin_ratio(0.05f), wd_margin_ratio(0.05f),
	gap_ratio(0.5f), title_ratio(0.7f),
	cbox_as(1.3f), cbox_tick_ratio(0.2f), cbox_gap_ratio(0.3f),
	// gl draw infos
	line_width(1.0f), char_num(0),
	line_vao(0), line_vbo(0), line_veo(0),
	tri_vao(0), tri_vbo(0), tri_veo(0),
	char_vao(0), char_vbo(0) {}

QtUniformColorMapObject::~QtUniformColorMapObject() { clear_gl_buffer(); }

void QtUniformColorMapObject::clear_gl_buffer()
{
	if (vbos[0] || vbos[1] || vbos[2] || vbos[3] || vbos[4])
	{
		gl.glDeleteBuffers(5, vbos);
		vbos[0] = 0;
		vbos[1] = 0;
		vbos[2] = 0;
		vbos[3] = 0;
		vbos[4] = 0;
	}
	if (vaos[0] || vaos[1] || vaos[2])
	{
		gl.glDeleteVertexArrays(3, vaos);
		vaos[0] = 0;
		vaos[1] = 0;
		vaos[2] = 0;
	}
}

void QtUniformColorMapObject::cal_str_geometry(
	const char *str,
	QtCharBitmapLoader& loader,
	GLint& font_wd, GLint& font_ht
	)
{
	const QtCharBitmapLoader::CharData* pref_char 
				= loader.get_char_data('1');
	font_ht = pref_char->get_advance_y();
	size_t str_len = strlen(str);
	font_wd = 0;
	for (size_t c_id = 0; c_id < str_len; ++c_id)
	{
		const QtCharBitmapLoader::CharData* pcd
				= loader.get_char_data(str[c_id]);
		if (!pcd)
			continue;
		font_wd += pcd->get_advance_x();
	}
}

void QtUniformColorMapObject::scale_rect_height(
	Rect& ori_rect,
	GLfloat sc,
	Rect& new_rect
	)
{
	new_rect.xl = ori_rect.xl;
	new_rect.xu = ori_rect.xu;

	GLfloat ori_rect_ht = ori_rect.yu - ori_rect.yl;
	GLfloat new_rect_ht = ori_rect_ht * sc;
	GLfloat y_padding = (ori_rect_ht - new_rect_ht) * 0.5;
	new_rect.yl = ori_rect.yl + y_padding;
	new_rect.yu = ori_rect.yu - y_padding;
}

void QtUniformColorMapObject::get_mid_align_str_pos(
	Rect& pos_rect,
	GLfloat str_as,
	Rect& str_rect
	)
{
	GLfloat pos_wd = pos_rect.xu - pos_rect.xl;
	GLfloat pos_ht = pos_rect.yu - pos_rect.yl;
	GLfloat pos_as = pos_wd / pos_ht;
	GLfloat pad;
	if (pos_as > str_as)
	{
		pad = (pos_wd - pos_ht * str_as) * 0.5;
		str_rect.xl = pos_rect.xl + pad;
		str_rect.xu = pos_rect.xu - pad;
		str_rect.yl = pos_rect.yl;
		str_rect.yu = pos_rect.yu;
	}
	else
	{
		pad = (pos_ht - pos_wd / str_as) * 0.5;
		str_rect.xl = pos_rect.xl;
		str_rect.xu = pos_rect.xu;
		str_rect.yl = pos_rect.yl + pad;
		str_rect.yu = pos_rect.yu - pad;
	}
}

void QtUniformColorMapObject::get_left_align_str_pos(
	Rect& pos_rect,
	GLfloat str_as,
	Rect& str_rect
	)
{
	GLfloat pos_wd = pos_rect.xu - pos_rect.xl;
	GLfloat pos_ht = pos_rect.yu - pos_rect.yl;
	GLfloat pos_as = pos_wd / pos_ht;
	if (pos_as > str_as)
	{
		str_rect.xl = pos_rect.xl;
		str_rect.xu = pos_rect.xl + pos_ht * str_as;
		str_rect.yl = pos_rect.yl;
		str_rect.yu = pos_rect.yu;
	}
	else
	{
		GLfloat pad = (pos_ht - pos_wd / str_as) * 0.5;
		str_rect.xl = pos_rect.xl;
		str_rect.xu = pos_rect.xu;
		str_rect.yl = pos_rect.yl + pad;
		str_rect.yu = pos_rect.yu - pad;
	}
}

void QtUniformColorMapObject::add_line(
	GLfloat x0, GLfloat y0,
	GLfloat x1, GLfloat y1
)
{
	Type0NodeData nd;
	nd.type = 0;
	nd.x = x0;
	nd.y = y0;
	line_node_mem.add(nd);
	nd.x = x1;
	nd.y = y1;
	line_node_mem.add(nd);
	GLuint n_id = line_elem_mem.get_num();
	line_elem_mem.add(n_id);
	++n_id;
	line_elem_mem.add(n_id);
}

void QtUniformColorMapObject::add_rect(
	Rect& rect,
	QVector3D& color
)
{
	GLuint n_id0 = tri_node_mem.get_num();
	GLuint n_id1 = n_id0 + 1;
	GLuint n_id2 = n_id0 + 2;
	GLuint n_id3 = n_id0 + 3;
	// 0, 1, 2
	tri_elem_mem.add(n_id0);
	tri_elem_mem.add(n_id1);
	tri_elem_mem.add(n_id2);
	// 0, 2, 3
	tri_elem_mem.add(n_id0);
	tri_elem_mem.add(n_id2);
	tri_elem_mem.add(n_id3);

	Type1NodeData nd;
	nd.type = 1;
	nd.r = color.x();
	nd.g = color.y();
	nd.b = color.z();
	// n0
	nd.x = rect.xl;
	nd.y = rect.yl;
	tri_node_mem.add(nd);
	// n1
	nd.x = rect.xu;
	nd.y = rect.yl;
	tri_node_mem.add(nd);
	// n2
	nd.x = rect.xu;
	nd.y = rect.yu;
	tri_node_mem.add(nd);
	// n3
	nd.x = rect.xl;
	nd.y = rect.yu;
	tri_node_mem.add(nd);
}

void QtUniformColorMapObject::add_char_texture(
	Rect& rect,
	GLuint tex_id
	)
{
	CharNodeData nd;
	// n0
	nd.x = rect.xl;
	nd.y = rect.yu;
	nd.tex_x = 0.0f;
	nd.tex_y = 0.0f;
	char_node_mem.add(nd);
	// n1
	nd.x = rect.xl;
	nd.y = rect.yl;
	nd.tex_x = 0.0f;
	nd.tex_y = 1.0f;
	char_node_mem.add(nd);
	// n2
	nd.x = rect.xu;
	nd.y = rect.yl;
	nd.tex_x = 1.0f;
	nd.tex_y = 1.0f;
	char_node_mem.add(nd);
	// n3
	nd.x = rect.xl;
	nd.y = rect.yu;
	nd.tex_x = 0.0f;
	nd.tex_y = 0.0f;
	char_node_mem.add(nd);
	// n4
	nd.x = rect.xu;
	nd.y = rect.yl;
	nd.tex_x = 1.0f;
	nd.tex_y = 1.0f;
	char_node_mem.add(nd);
	// n5
	nd.x = rect.xu;
	nd.y = rect.yu;
	nd.tex_x = 1.0f;
	nd.tex_y = 0.0f;
	char_node_mem.add(nd);

	char_textures.add(tex_id);

	++char_num;
}

void QtUniformColorMapObject::add_string_texture(
	GLfloat xpos,
	GLfloat ypos,
	GLfloat nr,
	const char* str
	)
{
	Rect char_pos;
	size_t str_len = strlen(str);
	for (size_t c_id = 0; c_id < str_len; ++c_id)
	{
		const QtCharBitmapLoader::CharData* cd
			= char_loader.get_char_data(str[c_id]);
		if (!cd)
			continue;

		char_pos.xl = xpos + GLfloat(cd->get_bearing_x()) * nr;
		char_pos.yl = ypos + (GLfloat(cd->get_bearing_y()) - GLfloat(cd->get_height())) * nr;
		char_pos.xu = char_pos.xl + GLfloat(cd->get_width()) * nr;
		char_pos.yu = char_pos.yl + GLfloat(cd->get_height()) * nr;
		add_char_texture(char_pos, cd->get_texture_id());

		xpos += cd->get_advance_x() * nr;
	}
}

int QtUniformColorMapObject::init(
	GLfloat _xpos, GLfloat _ypos, GLfloat ht,
	UniformColorMap& color_map,
	const char* title,
	const char* num_format,
	const char* ttf_filename
	)
{
	// ==== calculate size of each components ====
	cmap_ht = ht;

	// height margin
	GLfloat margin_ht = cmap_ht * ht_margin_ratio;
	content_ht = cmap_ht - 2.0f * margin_ht;

	// height of title, color boxes and gap between them
	size_t color_num = color_map.get_color_num();
	GLfloat divider_tmp = title_ratio + gap_ratio + GLfloat(color_num);
	num_ht = content_ht / divider_tmp;
	cbox_ht = num_ht;
	title_ht = num_ht * title_ratio;
	title_gap_ht = num_ht * gap_ratio;

	// width color box
	cbox_wd = cbox_ht * cbox_as;
	cbox_tick_wd = cbox_wd * cbox_tick_ratio;
	cbox_gap_wd = cbox_wd * cbox_gap_ratio;

	// load string render to get string info
	char_loader.load_ttf_file(ttf_filename, 60);

	GLint str_wd, str_ht;
	// title string size
	title_str_info.str = title;
	cal_str_geometry(title, char_loader, str_wd, str_ht);
	title_str_info.as = GLfloat(str_wd) / GLfloat(str_ht);
	title_str_info.font_pixel_wd = str_wd;
	title_wd = title_ht * title_str_info.as;

	// number string size
	entry_infos.resize(color_num);
	char num_str_buf[200];
	float num_start = color_map.get_lower_bound();
	float num_inv = (color_map.get_upper_bound() - num_start)
					/ float(color_num - 1);
	GLfloat num_wd_tmp;
	num_wd = 0.0; // get longest number string
	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		EntryInfo &ei = entry_infos[c_id];
		StrInfo& si = ei.num_str_info;
		// form num str
		ei.number = num_start + float(c_id) * num_inv;
		snprintf(num_str_buf, 200, num_format, ei.number);
		si.str = num_str_buf;
		// get num str geometry
		cal_str_geometry(num_str_buf, char_loader, str_wd, str_ht);
		si.as = GLfloat(str_wd) / GLfloat(str_ht);
		si.font_pixel_wd = str_wd;
		num_wd_tmp = num_ht * si.as * Num_Font_Size_Scale;
		if (num_wd < num_wd_tmp)
			num_wd = num_wd_tmp;
	}

	content_wd = cbox_wd + cbox_tick_wd + cbox_gap_wd + num_wd;
	if (content_wd < title_wd)
		content_wd = title_wd;
	title_wd = content_wd;

	cmap_wd = title_wd / (1.0f - 2.0f * wd_margin_ratio);
	GLfloat margin_wd = cmap_wd * wd_margin_ratio;

	// ==== calculate position of each components ====
	xpos = _xpos;
	ypos = _ypos;
	if (xpos < 0.0f)
		xpos = 0.0f;
	if (ypos < 0.0f)
		ypos = 0.0f;

	cmap_rect.xl = xpos;
	cmap_rect.xu = xpos + cmap_wd;
	cmap_rect.yl = ypos;
	cmap_rect.yu = ypos + cmap_ht;
	content_rect.xl = cmap_rect.xl + margin_wd;
	content_rect.xu = cmap_rect.xu - margin_wd;
	content_rect.yl = cmap_rect.yl + margin_ht;
	content_rect.yu = cmap_rect.yu - margin_ht;

	title_rect.xl = content_rect.xl;
	title_rect.xu = content_rect.xu;
	title_rect.yu = content_rect.yu;
	title_rect.yl = title_rect.yu - title_ht;
	// title string
	Rect str_rect_tmp;
	scale_rect_height(title_rect, Title_Font_Size_Scale, str_rect_tmp);
	get_mid_align_str_pos(str_rect_tmp, title_str_info.as, title_str_rect);
	title_str_info.nr = title_str_rect.width() / GLfloat(title_str_info.font_pixel_wd);

	cbar_rect.xl = content_rect.xl;
	cbar_rect.xu = cbar_rect.xl + cbox_wd;
	cbar_rect.yl = content_rect.yl;
	cbar_rect.yu = cbar_rect.yl + GLfloat(color_num) * cbox_ht;

	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		EntryInfo& ei = entry_infos[c_id];

		Rect& cbox_rect = ei.cbox_rect;
		cbox_rect.xl = cbar_rect.xl;
		cbox_rect.xu = cbar_rect.xu;
		cbox_rect.yl = content_rect.yl + GLfloat(c_id) * cbox_ht;
		cbox_rect.yu = cbox_rect.yl + cbox_ht;

		Rect& num_rect = ei.num_rect;
		num_rect.xl = cbox_rect.xu + cbox_tick_wd + cbox_gap_wd;
		num_rect.xu = num_rect.xl + num_wd;
		num_rect.yl = cbox_rect.yl;
		num_rect.yu = cbox_rect.yu;

		Rect& num_str_rect = ei.num_str_rect;
		StrInfo& si = ei.num_str_info;
		scale_rect_height(num_rect, Num_Font_Size_Scale, str_rect_tmp);
		get_mid_align_str_pos(str_rect_tmp, si.as, num_str_rect);
		si.nr = num_str_rect.width() / GLfloat(si.font_pixel_wd);
	}

	// form gl buffer
	line_node_mem.reserve((7 + color_num) * 2);
	line_elem_mem.reserve((7 + color_num) * 2);
	tri_node_mem.reserve(color_num * 2 * 4);
	tri_elem_mem.reserve(color_num * 2 * 6);
	char_node_mem.reserve((color_num + 1) * 20 * 6);
	char_textures.reserve((color_num + 1) * 20);

	// outer frame
	add_line(cmap_rect.xl, cmap_rect.yl, cmap_rect.xu, cmap_rect.yl);
	add_line(cmap_rect.xu, cmap_rect.yl, cmap_rect.xu, cmap_rect.yu);
	add_line(cmap_rect.xu, cmap_rect.yu, cmap_rect.xl, cmap_rect.yu);
	add_line(cmap_rect.xl, cmap_rect.yu, cmap_rect.xl, cmap_rect.yl);

	// title string
	add_string_texture(
		title_str_rect.xl,
		title_str_rect.yl,
		title_str_info.nr,
		title_str_info.str.c_str()
		);
	
	// color boxes
	UniformColorMap::Color *colors = color_map.get_colors();
	GLfloat cbox_tick_xu = cbar_rect.xu + cbox_tick_wd;
	add_line(cbar_rect.xl, cbar_rect.yl, cbar_rect.xl, cbar_rect.yu);
	add_line(cbar_rect.xu, cbar_rect.yl, cbar_rect.xu, cbar_rect.yu);
	add_line(cbar_rect.xl, cbar_rect.yl, cbox_tick_xu, cbar_rect.yl);
	
	QVector3D color_tmp;
	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		EntryInfo& ei = entry_infos[c_id];

		// color boxes
		Rect& cbox_rect = ei.cbox_rect;
		add_line(cbar_rect.xl, cbox_rect.yu, cbox_tick_xu, cbox_rect.yu);
		
		UniformColorMap::Color &color = colors[c_id];
		color_tmp.setX(color.r);
		color_tmp.setY(color.g);
		color_tmp.setZ(color.b);
		add_rect(cbox_rect, color_tmp);

		// number strings
		add_string_texture(
			ei.num_str_rect.xl,
			ei.num_str_rect.yl,
			ei.num_str_info.nr,
			ei.num_str_info.str.c_str());
	}

	// set up opengl buffer
	clear_gl_buffer();

	gl.glGenVertexArrays(3, vaos);
	gl.glGenBuffers(5, vbos);

	// line gl buffer
	gl.glBindVertexArray(line_vao);
	gl.glBindBuffer(GL_ARRAY_BUFFER, line_vbo);
	gl.glBufferData(GL_ARRAY_BUFFER,
		line_node_mem.get_num() * sizeof(Type0NodeData),
		(GLvoid *)line_node_mem.get_mem(),
		GL_STREAM_DRAW);

	// type
	gl.glVertexAttribIPointer(0, 1, GL_INT,
		sizeof(Type0NodeData), (GLvoid*)0);
	gl.glEnableVertexAttribArray(0);
	// coordinates
	gl.glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE,
		sizeof(Type0NodeData), (GLvoid*)offsetof(Type0NodeData, x));
	gl.glEnableVertexAttribArray(1);
	
	gl.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, line_veo);
	gl.glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		line_elem_mem.get_num() * sizeof(GLuint),
		(GLvoid *)line_elem_mem.get_mem(),
		GL_STREAM_DRAW);

	// color squares gl buffer
	gl.glBindVertexArray(tri_vao);
	gl.glBindBuffer(GL_ARRAY_BUFFER, tri_vbo);
	gl.glBufferData(GL_ARRAY_BUFFER,
		tri_node_mem.get_num() * sizeof(Type1NodeData),
		(GLvoid*)tri_node_mem.get_mem(),
		GL_STREAM_DRAW);

	// type
	gl.glVertexAttribIPointer(0, 1, GL_INT,
		sizeof(Type1NodeData), (GLvoid*)0);
	gl.glEnableVertexAttribArray(0);
	// coordinates
	gl.glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE,
		sizeof(Type1NodeData), (GLvoid*)offsetof(Type1NodeData, x));
	gl.glEnableVertexAttribArray(1);
	// color
	gl.glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE,
		sizeof(Type1NodeData), (GLvoid*)offsetof(Type1NodeData, r));
	gl.glEnableVertexAttribArray(2);

	gl.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, tri_veo);
	gl.glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		tri_elem_mem.get_num() * sizeof(GLuint),
		(GLvoid*)tri_elem_mem.get_mem(),
		GL_STREAM_DRAW);

	// for character gl buffer
	CharNodeData cnd[6];
	gl.glBindVertexArray(char_vao);
	gl.glBindBuffer(GL_ARRAY_BUFFER, char_vbo);
	gl.glBufferData(GL_ARRAY_BUFFER,
		6 * sizeof(CharNodeData),
		(GLvoid *)&cnd,
		GL_DYNAMIC_DRAW);

	// coordinates
	gl.glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE,
		sizeof(CharNodeData), (GLvoid*)0);
	gl.glEnableVertexAttribArray(0);
	// texture coordinates
	gl.glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE,
		sizeof(CharNodeData), (GLvoid*)offsetof(CharNodeData, tex_x));
	gl.glEnableVertexAttribArray(1);

	return 0;
}

void QtUniformColorMapObject::draw(
	QOpenGLShaderProgram& shader_plain2D,
	QOpenGLShaderProgram& shader_char)
{
	// draw chars
	shader_char.bind();
	shader_char.setUniformValue("g_color", char_color);
	shader_char.setUniformValue("char_texture", 0);

	gl.glEnable(GL_BLEND);
	gl.glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	gl.glBindVertexArray(char_vao);
	gl.glBindBuffer(GL_ARRAY_BUFFER, char_vbo);
	gl.glActiveTexture(GL_TEXTURE_2D);
	CharNodeData* cnds = char_node_mem.get_mem();
	GLuint* ctexs = char_textures.get_mem();
	for (size_t c_id = 0; c_id < char_num; ++c_id)
	{
		gl.glBufferSubData(GL_ARRAY_BUFFER, 0,
			6 * sizeof(CharNodeData),
			(GLvoid *)cnds);
		gl.glBindTexture(GL_TEXTURE_2D, ctexs[c_id]);
		gl.glDrawArrays(GL_TRIANGLES, 0, 6);
		cnds += 6;
	}

	gl.glDisable(GL_BLEND);

	shader_plain2D.bind();
	shader_plain2D.setUniformValue("g_color", char_color);

	// draw rects
	gl.glBindVertexArray(tri_vao);
	gl.glDrawElements(GL_TRIANGLES,
		tri_elem_mem.get_num(),
		GL_UNSIGNED_INT,
		(GLvoid*)0
		);
	
	// draw lines
	gl.glLineWidth(line_width);
	gl.glBindVertexArray(line_vao);
	gl.glDrawElements(GL_LINES,
		line_elem_mem.get_num(),
		GL_UNSIGNED_INT,
		(GLvoid*)0
		);
}
