#ifndef __Triangle_Mesh_H__
#define __Triangle_Mesh_H__

#include "TriangleMeshTemplate.hpp"

namespace TriangleMesh_Internal
{
	struct Node
	{
		size_t id;
		double x, y;
	};

	struct Element
	{
		size_t id;
		size_t n1, n2, n3;
		double vol;
	};

	struct Edge { size_t n1, n2; };
};

struct TriangleMesh : 
	public TriangleMeshTemplate<TriangleMesh_Internal::Node,
			TriangleMesh_Internal::Element, TriangleMesh_Internal::Edge>
{
public:
	typedef TriangleMesh_Internal::Node Node;
	typedef TriangleMesh_Internal::Element Element;
	typedef TriangleMesh_Internal::Edge Edge;
};

#endif