#ifndef __Qt_Char_Bitmap_Loader_h__
#define __Qt_Char_Bitmap_Loader_h__

#include <unordered_map>

#include <QOpenGLFunctions_3_3_Core>

#include <ft2build.h>
#include FT_FREETYPE_H

class QtCharBitmapLoader
{
public:
	struct CharData
	{
		friend QtCharBitmapLoader;
	protected:
		unsigned char id;
		char c;
		GLuint texture_id;
		GLuint width, height;
		GLuint bearing_x, bearing_y;
		GLint advance_x, advance_y;
	public:
		CharData() : texture_id(0) {}
		CharData(const CharData& other)
		{
			id = other.id;
			c = other.c;
			texture_id = other.texture_id;
			width = other.width;
			height = other.height;
			bearing_x = other.bearing_x;
			bearing_y = other.bearing_y;
			advance_x = other.advance_x;
			advance_y = other.advance_y;
		}
		~CharData() {}
		inline void clear(QOpenGLFunctions_3_3_Core& gl)
		{
			if (texture_id)
			{
				gl.glDeleteTextures(1, &texture_id);
				texture_id = 0;
			}
		}
		inline unsigned char get_id() { return id; }
		inline char get_char() { return c; }
		inline GLuint get_texture_id() { return texture_id; }
		inline GLuint get_width() { return width; }
		inline GLuint get_height() { return height; }
		inline GLuint get_bearing_x() { return bearing_x; }
		inline GLuint get_bearing_y() { return bearing_y; }
		inline GLint get_advance_x() { return advance_x; }
		inline GLint get_advance_y() { return advance_y; }
	};

protected:
	QOpenGLFunctions_3_3_Core& gl;

	size_t char_num;
	typedef std::unordered_map<char, CharData> CharDataMap;
	CharDataMap char_map;
	CharData *id_to_char[256];

	bool ft_is_init;
	FT_Library ft_lib;
	FT_Face ft_face;
	FT_UInt ft_font_ht; // in pixel

	CharData* load_char_data(char c);

public:
	QtCharBitmapLoader(QOpenGLFunctions_3_3_Core& _gl);
	~QtCharBitmapLoader();
	void clear();

	int load_ttf_file(const char *filename, unsigned int font_ht);

	CharData *get_char_data(char c);
	
	inline CharData* get_char_data_by_id(size_t id)
	{
		if (id < char_num)
			return id_to_char[id];
		return nullptr;
	}
};

#endif