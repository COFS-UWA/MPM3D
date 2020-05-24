#ifndef __Tetrahedron_Mesh_Buffer_h__
#define __Tetrahedron_Mesh_Buffer_h__

#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

class TetrahedronMeshBuffer
{
protected:
	QOpenGLFunctions_3_3_Core &glwp;
	GLuint vao, vbo, veo;
	GLsizei index_num;
	QVector3D color;

public:
	TetrahedronMeshBuffer(QOpenGLFunctions_3_3_Core &wp) : glwp(wp),
		vao(0), vbo(0), veo(0), index_num(0), color(1.0f, 1.0f, 1.0f) {}
	~TetrahedronMeshBuffer() { clear(); }

	void clear()
	{
		if (veo)
		{
			glwp.glDeleteBuffers(1, &veo);
			veo = 0;
		}

		if (vbo)
		{
			glwp.glDeleteBuffers(1, &vbo);
			vbo = 0;
		}

		if (vao)
		{
			glwp.glDeleteVertexArrays(1, &vao);
			vao = 0;
		}

		index_num = 0;
	}

	// use uniform color shader
	void draw(QOpenGLShaderProgram& shader)
	{
		if (vao)
		{
			shader.setUniformValue("color", color);

			glwp.glBindVertexArray(vao);
			glwp.glLineWidth(1.0f);
			glwp.glDrawElements(GL_LINES, index_num, GL_UNSIGNED_INT, 0);
		}
	}

	template <typename TMesh>
	int init(TMesh &mesh, QVector3D &_color)
	{
		typedef typename TMesh::Node Node;
		typedef typename TMesh::Element Element;
		clear();

		size_t node_num = mesh.get_node_num();
		Node* nodes = mesh.get_nodes();
		size_t elem_num = mesh.get_elem_num();
		Element* elems = mesh.get_elems();
		if (!nodes || node_num == 0 || !elems || elem_num == 0)
			return -1;

		color = _color;

		glwp.glGenVertexArrays(1, &vao);
		glwp.glBindVertexArray(vao);

		// node data
		glwp.glGenBuffers(1, &vbo);
		glwp.glBindBuffer(GL_ARRAY_BUFFER, vbo);
		size_t node_data_num = node_num * 3;
		GLfloat *node_data = new GLfloat[node_data_num];
		GLfloat *pn = node_data;
		for (size_t n_id = 0; n_id < node_num; ++n_id)
		{
			Node& n = nodes[n_id];
			pn[0] = GLfloat(n.x);
			pn[1] = GLfloat(n.y);
			pn[2] = GLfloat(n.z);
			pn += 3;
		}
		glwp.glBufferData(
			GL_ARRAY_BUFFER,
			node_data_num * sizeof(GLfloat),
			node_data,
			GL_STREAM_DRAW
		);
		delete[] node_data;

		glwp.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
		glwp.glEnableVertexAttribArray(0);

		// element data
		glwp.glGenBuffers(1, &veo);
		glwp.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo);
		index_num = elem_num * 12;
		GLuint *line_indices = new GLuint[index_num];
		GLuint *pl = line_indices;
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element& e = elems[e_id];
			// edge 1
			pl[0] = e.n1;
			pl[1] = e.n2;
			// edge 2
			pl[2] = e.n1;
			pl[3] = e.n3;
			// edge 3
			pl[4] = e.n1;
			pl[5] = e.n4;
			// edge 4
			pl[6] = e.n2;
			pl[7] = e.n3;
			// edge 5
			pl[8] = e.n2;
			pl[9] = e.n4;
			// edge 6
			pl[10] = e.n3;
			pl[11] = e.n4;
			pl += 12;
		}
		glwp.glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			index_num * sizeof(GLuint),
			line_indices,
			GL_STREAM_DRAW
		);
		delete[] line_indices;

		return 0;
	}
};

#endif