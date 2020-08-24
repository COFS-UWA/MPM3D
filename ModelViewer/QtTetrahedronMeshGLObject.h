#ifndef __Qt_Tetrahedron_Mesh_GL_Object_h__
#define __Qt_Tetrahedron_Mesh_GL_Object_h__

#include <iostream>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

#include "NumPairHashTable.hpp"

class QtTetrahedronMeshGLObject
{
protected:
	QOpenGLFunctions_3_3_Core &gl;

	struct NodeData
	{
		GLint type;
		GLfloat x, y, z;
	};

	struct EdgeData { GLuint n1, n2; };

	GLuint vao, vbo, veo;
	size_t edge_node_num;
	QVector3D color;

	void clear();

	inline static void sort_acc(size_t ids[4])
	{
		size_t tmp, min_id;
		for (size_t i = 0; i < 3; ++i)
		{
			min_id = i;
			for (size_t j = i + 1; j < 4; ++j)
			{
				if (ids[j] < ids[min_id])
					min_id = j;
			}
			if (min_id != i)
			{
				tmp = ids[min_id];
				ids[min_id] = ids[i];
				ids[i] = tmp;
			}
		}
	}
	
public:
	QtTetrahedronMeshGLObject(QOpenGLFunctions_3_3_Core& _gl);
	~QtTetrahedronMeshGLObject();

	// Node has members x, y
	// Edge has members n1, n2
	template <typename Node, typename Edge>
	int init_from_edges(Node* nodes, size_t node_num, Edge* edges, size_t edge_num, QVector3D& c);

	// Node has members x, y
	// Element has members n1, n2, n3, n4
	template <typename Node, typename Element>
	int init_from_elements(Node* nodes, size_t node_num, Element* elems, size_t elem_num, QVector3D& c);

	void draw(QOpenGLShaderProgram& shader);
};

template <typename Node, typename Edge>
int QtTetrahedronMeshGLObject::init_from_edges(
	Node* nodes, size_t node_num,
	Edge* edges, size_t edge_num,
	QVector3D& c
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
		nd.z = GLfloat(n.z);
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
		(GLvoid*)offsetof(NodeData, type)
		);
	gl.glEnableVertexAttribArray(0);
	// v_pos
	gl.glVertexAttribPointer(1,
		3, GL_FLOAT, GL_FALSE,
		sizeof(NodeData),
		(GLvoid*)offsetof(NodeData, x)
		);
	gl.glEnableVertexAttribArray(1);

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
int QtTetrahedronMeshGLObject::init_from_elements(
	Node* nodes,
	size_t node_num,
	Element* elems,
	size_t elem_num,
	QVector3D& c
)
{
	if (!nodes || !node_num || !elems || !elem_num)
		return -1;

	NumPairHashTable<size_t> table(node_num * 3);
	union
	{
		struct { size_t n1_id, n2_id, n3_id, n4_id; };
		size_t n_ids[4];
	};

	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		n1_id = e.n1;
		n2_id = e.n2;
		n3_id = e.n3;
		n4_id = e.n4;
		sort_acc(n_ids); // n1_id < n2_id < n3_id < n4_id
		table.add_pair(n1_id, n2_id);
		table.add_pair(n1_id, n3_id);
		table.add_pair(n1_id, n4_id);
		table.add_pair(n2_id, n3_id);
		table.add_pair(n2_id, n4_id);
		table.add_pair(n3_id, n4_id);
	}

	size_t edge_num = table.get_pair_num();
	if (!edge_num)
		return -1;
	EdgeData *edges = new EdgeData[edge_num];
	table.output_pairs((GLuint *)edges);
	int res = init_from_edges(nodes, node_num, edges, edge_num, c);
	delete[] edges;

	return res;
}

#endif