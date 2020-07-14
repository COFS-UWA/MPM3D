#include "ModelViewer_pcp.h"

#include <ft2build.h>
#include FT_FREETYPE_H

#include "Geometry.h"
#include "QtUniformcolorMapObject.h"

QtUniformColorMapObject::QtUniformColorMapObject(
	QOpenGLFunctions_3_3_Core& _gl) : gl(_gl),
	char_loader(_gl),
	ht_margin_ratio(0.05f), wd_margin_ratio(0.05f),
	gap_ratio(0.1f), title_ratio(1.5f),
	cbox_as(1.28f), cbox_tick_ratio(0.15f), cbox_gap_ratio(0.1f),
	line_width(2.0f),
	line_vao(0), line_vbo(0), line_veo(0),
	tri_vao(0), tri_vbo(0), tri_veo(0)
{

}

QtUniformColorMapObject::~QtUniformColorMapObject() { clear_gl_buffer(); }

void QtUniformColorMapObject::clear_gl_buffer()
{
	if (line_veo)
	{
		gl.glDeleteBuffers(1, &line_veo);
		line_veo = 0;
	}
	if (line_vbo)
	{
		gl.glDeleteBuffers(1, &line_vbo);
		line_vbo = 0;
	}
	if (line_vao)
	{
		gl.glDeleteVertexArrays(1, &line_vao);
		line_vao = 0;
	}
	if (tri_veo)
	{
		gl.glDeleteBuffers(1, &tri_veo);
		tri_veo = 0;
	}
	if (tri_vbo)
	{
		gl.glDeleteBuffers(1, &tri_vbo);
		tri_vbo = 0;
	}
	if (tri_vao)
	{
		gl.glDeleteVertexArrays(1, &tri_vao);
		tri_vao = 0;
	}
}

void QtUniformColorMapObject::add_line(
	GLfloat x0, GLfloat y0,
	GLfloat x1, GLfloat y1
	)
{
	MCNodeData nd;
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

void QtUniformColorMapObject::add_rect_color(
	Rect &rect,
	QVector3D& color
	)
{
	NodeData nd;
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
	tri_node_mem.add(nd);
	// n2
	nd.y = rect.yu;
	tri_node_mem.add(nd);
	// n3
	nd.x = rect.xl;
	tri_node_mem.add(nd);

	GLuint n_id0 = tri_elem_mem.get_num();
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
}

void QtUniformColorMapObject::add_char_texture(
	Rect &rect,
	GLint tex_id
	)
{
	NodeData nd;
	nd.type = 2;
	nd.tex_id = tex_id;
	// n0
	nd.x = rect.xl;
	nd.y = rect.yu;
	nd.tex_x = 0.0f;
	nd.tex_y = 0.0f;
	tri_node_mem.add(nd);
	// n1
	nd.y = rect.yl;
	nd.tex_x = 0.0f;
	nd.tex_y = 1.0f;
	tri_node_mem.add(nd);
	// n2
	nd.x = rect.xu;
	nd.tex_x = 1.0f;
	nd.tex_y = 1.0f;
	tri_node_mem.add(nd);
	// n3
	nd.y = rect.yu;
	nd.tex_x = 1.0f;
	nd.tex_y = 0.0f;
	tri_node_mem.add(nd);

	GLuint n_id0 = tri_elem_mem.get_num();
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
}

void QtUniformColorMapObject::add_string_texture(
	GLfloat xpos,
	GLfloat ypos,
	StrInfo& strif
	)
{
	Rect char_pos;
	GLfloat nr = strif.norm_ratio;
	size_t str_len = strif.str.size();
	const char* str = strif.str.c_str();
	for (size_t c_id = 0; c_id < str_len; ++c_id)
	{
		QtCharBitmapLoader::CharData *cd = char_loader.get_char_data(str[c_id]);
		if (!cd)
			continue;
		
		char_pos.xl = xpos + cd->get_bearing_x() * nr;
		char_pos.yl = ypos + (cd->get_bearing_y() - cd->get_height()) * nr;
		char_pos.xu = char_pos.xl + cd->get_width() * nr;
		char_pos.yu = char_pos.yl + cd->get_height() * nr;

		add_char_texture(char_pos, cd->get_texture_id());

		xpos += cd->get_advance_x() * strif.norm_ratio;
	}
}

void QtUniformColorMapObject::get_str_size(
	StrInfo& strif,
	QtCharBitmapLoader& loader
	)
{
	strif.bm_wd = 0;
	strif.bm_ht = 0;
	std::string& str = strif.str;
	size_t str_len = str.size();
	const char* str_chars = str.c_str();
	for (size_t c_id = 0; c_id < str_len; ++c_id)
	{
		QtCharBitmapLoader::CharData *pcd = loader.get_char_data(str_chars[c_id]);
		if (!pcd)
			continue;
		
		if (strif.bm_wd < pcd->get_advance_y())
			strif.bm_wd = pcd->get_advance_y();
		
		strif.bm_ht += pcd->get_advance_x();
	}
	strif.bm_as = GLfloat(strif.bm_wd) / GLfloat(strif.bm_ht);
}

int QtUniformColorMapObject::init(
	GLfloat _xpos, GLfloat _ypos, GLfloat ht,
	UniformColorMap& color_map,
	const char* title,
	const char* num_format,
	const char* ttf_filename,
	int win_wd, int win_ht
	)
{
	// ==== calculate size of each components ====
	cmap_ht = ht;

	// margin
	GLfloat margin_ht = cmap_ht * ht_margin_ratio;
	content_ht = cmap_ht - 2.0f * margin_ht;

	// height of title, color boxes and gap between them
	size_t color_num = color_map.get_color_num();
	GLfloat divider_tmp = title_ratio + gap_ratio + GLfloat(color_num);
	num_ht = content_ht / divider_tmp;
	title_ht = num_ht * title_ratio;
	title_gap_ht = num_ht * gap_ratio;
	cbox_ht = num_ht;

	// color box geometry
	cbox_wd = cbox_ht * cbox_as;
	cbox_tick_wd = cbox_wd * cbox_tick_ratio;
	cbox_gap_wd = cbox_wd * cbox_gap_ratio;

	// load string render
	char_loader.load_ttf_file(ttf_filename, cbox_ht * win_ht);
	(void)win_wd;

	// title string size
	title_str_info.str = title;
	get_str_size(title_str_info, char_loader);
	title_str_info.ht = title_ht;
	title_str_info.wd = title_ht * title_str_info.bm_as;

	// number string size
	num_str_infos.reserve(color_num);
	char num_str_buf[200];
	float num_start, num_inv;
	num_start = color_map.get_lower_bound();
	num_inv = (color_map.get_upper_bound() - num_start)
			/ float(color_num);
	GLfloat num_wd = 0.0; // get longest number string
	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		StrInfo& num_si = num_str_infos[c_id];
		// form str
		num_si.number = num_start + float(c_id) * num_inv;
		snprintf(num_str_buf, 200, num_format, num_si.number);
		num_si.str = num_str_buf;
		// get size
		get_str_size(num_si, char_loader);
		num_si.ht = cbox_ht;
		num_si.wd = num_si.ht * num_si.bm_as;
		num_si.norm_ratio = num_si.ht / GLfloat(num_si.bm_ht);
		if (num_wd < num_si.wd)
			num_wd = num_si.wd;
	}

	content_wd = cbox_wd + cbox_tick_wd + cbox_gap_wd + num_wd;
	content_wd = content_wd > title_wd ? content_wd : title_wd;
	title_wd = content_wd;

	GLfloat margin_wd = cmap_wd * wd_margin_ratio;
	cmap_wd = margin_wd * 2.0f + content_wd;

	// ==== calculate position of each components ====
	xpos = _xpos;
	ypos = _ypos;
	// adjust xpos, ypos
	if (xpos < 0.0f)
		xpos = 0.0f;
	if (xpos + cmap_wd > 1.0f)
		xpos = 1.0f - cmap_wd;
	if (ypos < 0.0f)
		ypos = 0.0f;
	if (ypos + cmap_ht > 1.0f)
		ypos = 1.0f - cmap_ht;

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
	cbar_rect.xl = content_rect.xl;
	cbar_rect.xu = cbar_rect.xl + cbox_wd;
	cbar_rect.yl = content_rect.yl;
	cbar_rect.yu = cbar_rect.yl + GLfloat(color_num) * cbox_ht;
	cboxes_rect.reserve(color_num);
	nums_rect.reserve(color_num);
	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		Rect& cbox_rect = cboxes_rect[c_id];
		cbox_rect.xl = cbar_rect.xl;
		cbox_rect.xu = cbar_rect.xu;
		cbox_rect.yl = content_rect.yl + GLfloat(c_id) * cbox_ht;
		cbox_rect.yu = cbox_rect.yl + cbox_ht;
		Rect& num_rect = nums_rect[c_id];
		num_rect.xl = cbox_rect.xu + cbox_tick_wd + cbox_gap_wd;
		num_rect.xu = num_rect.xl + num_wd;
		num_rect.yl = cbox_rect.yl;
		num_rect.yu = cbox_rect.yu;
	}

	// form gl buffer
	line_node_mem.reserve(color_num*10);
	line_elem_mem.reserve(color_num*10);
	tri_node_mem.reserve(color_num * 10);
	tri_elem_mem.reserve(color_num * 10);

	// outer frame
	add_line(cmap_rect.xl, cmap_rect.yl, cmap_rect.xu, cmap_rect.yl);
	add_line(cmap_rect.xu, cmap_rect.yl, cmap_rect.xu, cmap_rect.yu);
	add_line(cmap_rect.xu, cmap_rect.yu, cmap_rect.xl, cmap_rect.yu);
	add_line(cmap_rect.xl, cmap_rect.yu, cmap_rect.xl, cmap_rect.yl);

	// color boxes
	UniformColorMap::Color *colors = color_map.get_colors();
	GLfloat cbox_tick_xu = cbar_rect.xu + cbox_tick_wd;
	add_line(cbar_rect.xl, cbar_rect.yl, cbar_rect.xl, cbar_rect.yu);
	add_line(cbar_rect.xu, cbar_rect.yl, cbar_rect.xu, cbar_rect.yu);
	add_line(cbar_rect.xl, cbar_rect.yl, cbox_tick_xu, cbar_rect.yl);
	QVector3D color_tmp;
	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		Rect& cbox_rect = cboxes_rect[c_id];
		add_line(cbar_rect.xl, cbox_rect.yu, cbox_tick_xu, cbox_rect.yu);
		UniformColorMap::Color &color = colors[c_id];
		color_tmp.setX(color.r);
		color_tmp.setY(color.g);
		color_tmp.setZ(color.b);
		add_rect_color(cbox_rect, color_tmp);
	}

	// title
	Rect title_str_rect;
	// get_mid_align_str_pos(confined_rect, title_str_rect, 0.9f)
	add_string_texture(
		title_str_rect.xl,
		title_str_rect.yl,
		title_str_info
		);

	// numbers
	Rect num_str_rect;
	for (size_t c_id = 0; c_id < color_num; ++c_id)
	{
		// get_left_align_str_pos(confined_rect, num_str_rect, 0.9f);
		add_string_texture(
			num_str_rect.xl,
			num_str_rect.yl,
			num_str_infos[c_id]
		);
	}

	// set up opengl buffer
	clear_gl_buffer();

	gl.glGenVertexArrays(1, &line_vao);
	gl.glBindVertexArray(line_vao);
	// ...
	//GLuint line_vao, line_vbo, line_veo;
	//GLuint tri_vao, tri_vbo, tri_veo;


	return 0;
}

void QtUniformColorMapObject::draw(QOpenGLShaderProgram& shader)
{
	//shader.setUniformValue();
	QVector3D white(1.0f, 1.0f, 1.0f);

	
}
