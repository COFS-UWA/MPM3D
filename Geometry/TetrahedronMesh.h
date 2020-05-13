#ifndef __Tetrahedron_Mesh_H__
#define __Tetrahedron_Mesh_H__

#include "TetrahedronMeshTemplate.hpp"

namespace TetrahedronMesh_Internal
{
	struct Node
	{
		size_t id;
		double x, y, z;
	};

	struct Element
	{
		size_t id;
		size_t n1, n2, n3, n4;
		double vol;
	};

	struct Edge { size_t n1, n2; };
};

struct TetrahedronMesh : 
	public TetrahedronMeshTemplate<TetrahedronMesh_Internal::Node,
			TetrahedronMesh_Internal::Element, TetrahedronMesh_Internal::Edge>
{
public:
	typedef TetrahedronMesh_Internal::Node Node;
	typedef TetrahedronMesh_Internal::Element Element;
	typedef TetrahedronMesh_Internal::Edge Edge;
};

#endif