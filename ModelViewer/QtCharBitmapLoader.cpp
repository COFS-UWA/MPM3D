#include "ModelViewer_pcp.h"

#define KEEP_ID_TO_CHAR_NUM
#include "QtCharBitmapLoader.h"

QtCharBitmapLoader::QtCharBitmapLoader(QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), char_num(0), ft_is_init(false)
{
	memset(id_to_char, 0, sizeof(CharData*) * ID_TO_CHAR_NUM);
}

QtCharBitmapLoader::~QtCharBitmapLoader() { clear(); }

void QtCharBitmapLoader::clear()
{
	// delete opengl textures
	for (CharDataMap::iterator iter = char_map.begin();
		iter != char_map.end(); ++iter)
		iter->second.clear(gl);

	if (ft_is_init)
	{
		FT_Done_Face(ft_face);
		FT_Done_FreeType(ft_lib);
		ft_is_init = false;
	}
}

int QtCharBitmapLoader::load_ttf_file(
	const char* filename,
	unsigned int font_ht
)
{
	if (FT_Init_FreeType(&ft_lib))
		return -1;

	if (FT_New_Face(ft_lib, filename, 0, &ft_face))
	{
		FT_Done_FreeType(ft_lib);
		return -2;
	}

	ft_font_ht = font_ht;
	if (FT_Set_Pixel_Sizes(ft_face, 0, font_ht)) // need better value
	{
		FT_Done_Face(ft_face);
		FT_Done_FreeType(ft_lib);
		return -3;
	}

	// completed
	ft_is_init = true;
	return 0;
}

QtCharBitmapLoader::CharData* QtCharBitmapLoader::load_char_data(char c)
{
	if (c <= 31 || c == 127)
		return nullptr;

	if (!ft_is_init)
		return nullptr;

	if (char_num >= ID_TO_CHAR_NUM) // exceed the uniform length limit
		return nullptr;

	FT_UInt char_index = FT_Get_Char_Index(ft_face, c);
	if (FT_Load_Char(ft_face, char_index, FT_LOAD_RENDER))
		return nullptr;
	CharData cd;
	cd.id = char_num;
	cd.c = c;
	FT_GlyphSlot& char_glyph = ft_face->glyph;
	cd.bearing_x = char_glyph->bitmap_left;
	cd.bearing_y = char_glyph->bitmap_top;
	FT_Bitmap& char_bitmap = char_glyph->bitmap;
	cd.height = char_bitmap.rows;
	cd.width = char_bitmap.width;
	FT_Vector& adv = char_glyph->advance;
	cd.advance_x = adv.x >> 6; // 1/64 pixels
	cd.advance_y = adv.y >> 6;
	// generate bitmap
	gl.glGenTextures(1, &cd.texture_id);
	gl.glBindTexture(GL_TEXTURE_2D, cd.texture_id);
	gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	gl.glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	gl.glTexImage2D(
		GL_TEXTURE_2D,
		0,
		GL_RED,
		char_bitmap.width,
		char_bitmap.rows,
		0,
		GL_RED,
		GL_UNSIGNED_BYTE,
		char_bitmap.buffer
		);
	gl.glGenerateMipmap(GL_TEXTURE_2D);

	glBindTexture(GL_TEXTURE_2D, 0);

	auto res = char_map.emplace(c, cd);
	CharData* res_cd = &res.first->second;
	id_to_char[char_num] = res_cd;
	++char_num;
	return res_cd;
}

QtCharBitmapLoader::CharData* QtCharBitmapLoader::get_char_data(char c)
{
	CharDataMap::iterator iter = char_map.find(c);
	if (iter != char_map.end())
		return &(iter->second);

	return load_char_data(c);
}
