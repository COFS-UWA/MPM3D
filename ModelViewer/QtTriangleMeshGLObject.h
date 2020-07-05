#ifndef __Qt_Triangle_Mesh_GL_Object_h__
#define __Qt_Triangle_Mesh_GL_Object_h__

#include <QColor>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

class QtTriangleMeshGLObject
{
protected:
	QOpenGLFunctions_3_3_Core& gl;

	struct NodeData
	{
		GLint type;
		GLfloat x, y;
	};

	struct EdgeData { GLuint n1, n2; };

	GLuint vao, vbo, veo;
	size_t edge_node_num;
	QVector3D color;

	void clear();

public:
	QtTriangleMeshGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtTriangleMeshGLObject();

	// Node has members x, y
	// Edge has members n1, n2
	template <typename Node, typename Edge>
	int init(Node *nodes, size_t node_num, Edge *edges, size_t edge_num, QVector3D &c);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename Node, typename Edge>
int QtTriangleMeshGLObject::init(
	Node* nodes, size_t node_num,
	Edge* edges, size_t edge_num,
	QVector3D &c
	)
{
	color = c;
	edge_node_num = edge_num * 2;

	gl.glGenVertexArrays(1, &vao);
	gl.glBindVertexArray(vao);

	gl.glGenBuffers(1, &vbo);
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);
	NodeData* node_datas = new NodeData[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		NodeData& nd = node_datas[n_id];
		nd.type = 0; // monocolor
		nd.x = GLfloat(n.x);
		nd.y = GLfloat(n.y);
	}
	gl.glBufferData(
		GL_ARRAY_BUFFER,
		node_num * sizeof(NodeData),
		node_datas,
		GL_STATIC_DRAW
		);
	delete[] node_datas;
	
	// v_type
	gl.glVertexAttribPointer(0,
		1, GL_INT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid *)offsetof(NodeData, type)
		);
	gl.glEnableVertexAttribArray(0);
	// v_pos
	gl.glVertexAttribPointer(1,
		2, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid *)offsetof(NodeData, x)
		);
	gl.glEnableVertexAttribArray(1);
	// v_color (not used)
	gl.glVertexAttribPointer(2,
		3, GL_FLOAT, GL_FALSE,
		0, (GLvoid *)0
		);
	gl.glEnableVertexAttribArray(2);

	gl.glGenBuffers(1, &veo);
	gl.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo);
	EdgeData* edge_datas = new EdgeData[edge_num];
	for (size_t e_id = 0; e_id < edge_num; ++e_id)
	{
		Edge& e = edges[e_id];
		EdgeData& ed = edge_datas[e_id];
		ed.n1 = GLuint(e.n1);
		ed.n2 = GLuint(e.n2);
	}
	gl.glBufferData(
		GL_ELEMENT_ARRAY_BUFFER,
		edge_num * sizeof(EdgeData),
		edge_datas,
		GL_STATIC_DRAW
		);
	delete[] edge_datas;

	return 0;
}

#endif