#ifndef __Extract_Surface_From_T3D_Mesh_h__
#define __Extract_Surface_From_T3D_Mesh_h__

#include <unordered_map>

class ExtractSurfaceFromT3DMesh
{
public:
	struct Face { size_t n1, n2, n3; };

	explicit ExtractSurfaceFromT3DMesh();
	~ExtractSurfaceFromT3DMesh();

	inline void clear() noexcept
	{
		if (faces)
		{
			delete[] faces;
			faces = nullptr;
		}
		face_num = 0;
	}

	inline Face* get_faces() noexcept { return faces; }
	inline const Face* get_faces() const noexcept { return faces; }
	inline size_t get_face_num() const noexcept { return face_num; }
	inline size_t get_face_id(const Face& f) const noexcept { return &f - faces; }

	// Tetrahedron has member n1, n2, n3, n4
	template <typename Tetrahedron>
	void init_from_tetrahedron(Tetrahedron *tehs, size_t teh_num)
	{
		clear();
		FaceMap face_map;
		union
		{
			Face face_key;
			char fkey_str[25];
		};
		fkey_str[24] = '\0';
		size_t f_id, inside_face_num = 0;
		for (f_id = 0; f_id < teh_num; ++f_id)
		{
			Tetrahedron &teh = tehs[f_id];
			// face 1
			face_key.n1 = teh.n1;
			face_key.n2 = teh.n2;
			face_key.n3 = teh.n3;
			inside_face_num += try_adding_teh_to_map(face_map, face_key);
			// face 2
			face_key.n1 = teh.n1;
			face_key.n2 = teh.n4;
			face_key.n3 = teh.n2;
			inside_face_num += try_adding_teh_to_map(face_map, face_key);
			// face 3
			face_key.n1 = teh.n2;
			face_key.n2 = teh.n4;
			face_key.n3 = teh.n3;
			inside_face_num += try_adding_teh_to_map(face_map, face_key);
			// face 4
			face_key.n1 = teh.n1;
			face_key.n2 = teh.n3;
			face_key.n3 = teh.n4;
			inside_face_num += try_adding_teh_to_map(face_map, face_key);
		}
		
		face_num = teh_num * 4 - inside_face_num * 2;
		faces = new Face[face_num];
		f_id = 0;
		for (auto fiter = face_map.begin(); fiter != face_map.end(); ++fiter)
		{
			if ((fiter->second & 0x03) == 0)
			{
				Face& f = faces[f_id];
				const size_t* ns = reinterpret_cast<const size_t*>(fiter->first.c_str());
				f.n1 = ns[(fiter->second >> 2) & 0x3];
				f.n2 = ns[(fiter->second >> 4) & 0x3];
				f.n3 = ns[(fiter->second >> 6) & 0x3];
				++f_id;
			}
		}
	}

protected:
	Face* faces;
	size_t face_num;

	using FaceMap = std::unordered_map<std::string, unsigned char>;
	
	// if face already in map, return 1
	// otherwise return 0
	static size_t try_adding_teh_to_map(FaceMap& face_map, Face& face);
};

#endif