#ifndef __Extract_Edge_From_T2D_Mesh_h__
#define __Extract_Edge_From_T2D_Mesh_h__

#include <unordered_map>

class ExtractEdgeFromT2DMesh
{
public:
	struct Edge
	{
		size_t n1, n2;
		inline bool operator==(const Edge& other) const noexcept
		{ return n1 == other.n1 && n2 == other.n2; }
	};
	
	explicit ExtractEdgeFromT2DMesh();
	~ExtractEdgeFromT2DMesh();

	inline void clear() noexcept
	{
		if (boundary_edges)
		{
			delete[] boundary_edges;
			boundary_edges = nullptr;
		}
		boundary_edge_num = 0;
	}

	inline Edge *get_boundary_edges() noexcept { return boundary_edges; }
	inline const Edge *get_boundary_edges() const noexcept { return boundary_edges; }
	inline size_t get_boundary_edge_num() const noexcept { return boundary_edge_num; }
	inline size_t get_edge_id(const Edge& e) const noexcept { return &e - boundary_edges; }

	// Tetrahedron has member n1, n2, n3
	template <typename Triangle>
	void init_from_triangle(Triangle *tris, size_t tri_num)
	{
		clear();
		EdgeMap edge_map;
		Edge edge_key;
		size_t inside_edge_num = 0;
		for (size_t t_id = 0; t_id < tri_num; ++t_id)
		{
			Triangle& tri = tris[t_id];
			// face 1
			edge_key.n1 = tri.n1;
			edge_key.n2 = tri.n2;
			inside_edge_num += try_adding_tri_to_map(edge_map, edge_key);
			// face 2
			edge_key.n1 = tri.n2;
			edge_key.n2 = tri.n3;
			inside_edge_num += try_adding_tri_to_map(edge_map, edge_key);
			// face 3
			edge_key.n1 = tri.n3;
			edge_key.n2 = tri.n1;
			inside_edge_num += try_adding_tri_to_map(edge_map, edge_key);
		}
		
		// the
		boundary_edge_num = tri_num * 3 - inside_edge_num * 2;
		boundary_edges = new Edge[boundary_edge_num];
		size_t e_id = 0;
		for (auto fiter = edge_map.begin(); fiter != edge_map.end(); ++fiter)
		{
			if ((fiter->second & 0x03) == 0)
			{
				Edge &e = boundary_edges[e_id];
				const size_t *e_ns = reinterpret_cast<const size_t *>(&fiter->first);
				e.n1 = e_ns[(fiter->second >> 2) & 0x3];
				e.n2 = e_ns[(fiter->second >> 4) & 0x3];
				++e_id;
			}
		}
	}

protected:
	Edge *boundary_edges;
	size_t boundary_edge_num;

	// hash edge node index as string
	struct EdgeHasher
	{
		size_t operator()(const Edge& e) const
		{
			size_t h = 0;
			h ^= std::hash<size_t>{}(e.n1) + 0x9e3779b9 + (h << 6) + (h >> 2);
			h ^= std::hash<size_t>{}(e.n2) + 0x9e3779b9 + (h << 6) + (h >> 2);
			return h;
		}
	};
	using EdgeMap = std::unordered_map<Edge, unsigned char, EdgeHasher>;
	
	// if edge already in map, return 1
	// otherwise return 0
	static size_t try_adding_tri_to_map(EdgeMap& edge_map, Edge& edge);
};

#endif