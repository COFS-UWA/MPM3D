#ifndef __Qt_Triangle_Mesh_GL_Object_h__
#define __Qt_Triangle_Mesh_GL_Object_h__

#include <QColor>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "NumPairHashTable.hpp"

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

	inline static void sort_acc(size_t ids[3])
	{
		size_t tmp, min_id;
		min_id = 0;
		if (ids[0] > ids[1])
			min_id = 1;
		if (ids[0] > ids[2])
			min_id = 2;
		if (min_id != 0)
		{
			tmp = ids[0];
			ids[0] = ids[min_id];
			ids[min_id] = tmp;
		}
		if (ids[1] > ids[2])
		{
			tmp = ids[1];
			ids[1] = ids[min_id];
			ids[min_id] = tmp;
		}
	}

public:
	QtTriangleMeshGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtTriangleMeshGLObject();

	// Node has members x, y
	// Edge has members n1, n2
	template <typename Node, typename Edge>
	int init_from_edges(Node *nodes, size_t node_num, Edge *edges, size_t edge_num, QVector3D &c);
	
	// Node has members x, y
	// Edge has members n1, n2, n3
	template <typename Node, typename Element>
	int init_from_elements(Node* nodes, size_t node_num, Element* elems, size_t elem_num, QVector3D& c);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename Node, typename Edge>
int QtTriangleMeshGLObject::init_from_edges(
	Node* nodes, size_t node_num,
	Edge* edges, size_t edge_num,
	QVector3D &c
	)
{
	if (!nodes || !node_num ||
		!edges || !edge_num)
		return -1;

	color = c;
	edge_node_num = edge_num * 2;

	gl.glGenVertexArrays(1, &vao);
	if (vao == 0)
		return -2;
	gl.glBindVertexArray(vao);

	gl.glGenBuffers(1, &vbo);
	if (vbo == 0)
		return -2;
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
	gl.glVertexAttribIPointer(0,
		1, GL_INT,
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
	if (veo == 0)
		return -2;
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

template <typename Node, typename Element>
int QtTriangleMeshGLObject::init_from_elements(
	Node* nodes,
	size_t node_num,
	Element *elems,
	size_t elem_num,
	QVector3D& c
	)
{
	if (!nodes || !node_num ||
		!elems || !elem_num)
		return -1;

	NumPairHashTable<size_t> table(node_num * 2);
	union
	{
		struct { size_t n1_id, n2_id, n3_id; };
		size_t n_ids[3];
	};

	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		n1_id = e.n1;
		n2_id = e.n2;
		n3_id = e.n3;
		sort_acc(n_ids); // n1_id < n2_id < n3_id
		table.add_pair(n1_id, n2_id);
		table.add_pair(n1_id, n3_id);
		table.add_pair(n2_id, n3_id);
	}

	size_t edge_num = table.get_pair_num();
	if (!edge_num)
		return -1;
	struct EdgeData { size_t n1, n2; } *edges;
	edges = new EdgeData[edge_num];
	table.output_pairs((size_t*)edges);
	int res = init_from_edges(nodes, node_num, edges, edge_num, c);
	delete[] edges;
	return res;
}

#endif