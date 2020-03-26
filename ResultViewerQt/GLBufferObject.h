#ifndef __GL_Buffer_Object_H__
#define __GL_Buffer_Object_H__

#include <QOpenGLFunctions_3_3_Core>

struct GLBufferObject
{
protected:
    QOpenGLFunctions_3_3_Core &parent;
    GLuint vao, vbo, veo;
	GLsizeiptr vert_data_num;
	GLsizeiptr elem_data_num;
	
public:
	GLBufferObject(QOpenGLFunctions_3_3_Core *pa);
	~GLBufferObject();

	void init_vert(GLfloat *node_data, size_t data_num);
	void set_attr_pointer();

	void init_elem(GLuint *indices, size_t indices_num);
	void clear(void);

    inline void use(void) const { parent.glBindVertexArray(vao); }
	inline void unuse(void) const
	{
		parent.glBindVertexArray(0);
		parent.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	}
};

#endif