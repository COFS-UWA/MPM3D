#ifndef __Mono_Color_Cube_Particle_System_h__
#define __Mono_Color_Cube_Particle_System_h__

#include <QOpenGLWidget> // for opengl type

#include "ItemArray.hpp"

// Generate opengl buffer for particles
class MonoColorCubeParticleSystem
{
protected:
	size_t pcl_num;
	GLfloat vol_scale;

	struct VertData
	{
		GLfloat x1, y1, z1;
		GLfloat x2, y2, z2;
		GLfloat x3, y3, z3;
		GLfloat x4, y4, z4;
		GLfloat x5, y5, z5;
		GLfloat x6, y6, z6;
		GLfloat x7, y7, z7;
		GLfloat x8, y8, z8;
		// Each particle has 8 nodes
		inline void init(
			GLfloat xc,
			GLfloat yc,
			GLfloat zc,
			GLfloat vol)
		{
			GLfloat half_len = pow(vol, 0.33333333f) * 0.5;
			// v1
			x1 = xc + half_len;
			y1 = yc + half_len;
			z1 = zc + half_len;
			// v2
			x2 = xc - half_len;
			y2 = yc + half_len;
			z2 = zc + half_len;
			// v3
			x3 = xc - half_len;
			y3 = yc - half_len;
			z3 = zc + half_len;
			// v4
			x4 = xc + half_len;
			y4 = yc - half_len;
			z4 = zc + half_len;
			// v5
			x5 = xc + half_len;
			y5 = yc + half_len;
			z5 = zc - half_len;
			// v6
			x6 = xc - half_len;
			y6 = yc + half_len;
			z6 = zc - half_len;
			// v7
			x7 = xc - half_len;
			y7 = yc - half_len;
			z7 = zc - half_len;
			// v8
			x8 = xc + half_len;
			y8 = yc - half_len;
			z8 = zc - half_len;
		}
	};
	MemoryUtils::ItemArray<VertData> vert_data;
	GLsizei vert_data_size;

	struct ElemData
	{
		GLuint n11, n12, n13;
		GLuint n21, n22, n23;
		GLuint n31, n32, n33;
		GLuint n41, n42, n43;
		GLuint n51, n52, n53;
		GLuint n61, n62, n63;
		GLuint n71, n72, n73;
		GLuint n81, n82, n83;
		GLuint n91, n92, n93;
		GLuint n101, n102, n103;
		GLuint n111, n112, n113;
		GLuint n121, n122, n123;
		// Each cube has 12 triangle faces
		inline void offset(GLuint off)
		{
			// tri 1
			n11 = 0 + off;
			n12 = 1 + off;
			n13 = 2 + off;
			// tri 2
			n21 = 0 + off;
			n22 = 2 + off;
			n23 = 3 + off;
			// tri 3
			n31 = 0 + off;
			n32 = 5 + off;
			n33 = 1 + off;
			// tri 4
			n41 = 0 + off;
			n42 = 4 + off;
			n43 = 5 + off;
			// tri 5
			n51 = 0 + off;
			n52 = 3 + off;
			n53 = 7 + off;
			// tri 6
			n61 = 0 + off;
			n62 = 7 + off;
			n63 = 4 + off;
			// tri 7
			n71 = 2 + off;
			n72 = 1 + off;
			n73 = 5 + off;
			// tri 8
			n81 = 2 + off;
			n82 = 5 + off;
			n83 = 6 + off;
			// tri 9
			n91 = 2 + off;
			n92 = 6 + off;
			n93 = 7 + off;
			// tri 10
			n101 = 2 + off;
			n102 = 7 + off;
			n103 = 3 + off;
			// tri 11
			n111 = 5 + off;
			n112 = 7 + off;
			n113 = 6 + off;
			// tri 12
			n121 = 5 + off;
			n122 = 4 + off;
			n123 = 7 + off;
		}
	};
	MemoryUtils::ItemArray<ElemData> elem_data; // face id
	GLsizei elem_data_size;

public:
	MonoColorCubeParticleSystem() : vert_data_size(0), elem_data_size(0) {}
	~MonoColorCubeParticleSystem() { clear(); }
	inline void clear()
	{
		vert_data.clear();
		elem_data.clear();
	}

	inline GLvoid *get_vert_data() { return (GLvoid *)vert_data.get_mem(); }
	inline GLsizei get_vert_data_size() { return vert_data_size; }
	inline GLvoid *get_elem_data() { return (GLvoid *)elem_data.get_mem(); }
	inline GLsizei get_elem_data_size() { return elem_data_size; }

	// Assumptions: Particle has memebers x, y, z, vol
	template <typename Particle>
	int init_data(
		Particle *pcls,
		size_t _pcl_num,
		GLfloat _vol_scale = 0.125f
		)
	{
		if (!pcls || _pcl_num == 0)
			return -1;

		pcl_num = _pcl_num;
		vol_scale = _vol_scale;

		vert_data_size = pcl_num * sizeof(VertData);
		vert_data.reserve(pcl_num);
		VertData *verts = vert_data.get_mem();

		elem_data_size = pcl_num * sizeof(ElemData);
		elem_data.reserve(pcl_num);
		ElemData *elems = elem_data.get_mem();
		
		size_t cur_id_off = 0;
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			Particle &pcl = pcls[pcl_id];
			GLfloat pcl_vol = GLfloat(pcl.get_vol());

			VertData &v = verts[pcl_id];
			v.init((GLfloat)pcl.x, (GLfloat)pcl.y, (GLfloat)pcl.z,
				   pcl_vol * vol_scale);
			
			ElemData &e = elems[pcl_id];
			e.offset(cur_id_off);

			cur_id_off += 8;
		}
	
		return 0;
	}

	// only update vertext buffer
	template <typename Particle>
	int update_data(Particle* pcls)
	{
		if (!pcls || pcl_num == 0)
			return -1;

		VertData* verts = vert_data.get_mem();
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			Particle& pcl = pcls[pcl_id];
			VertData& v = verts[pcl_id];
			v.init((GLfloat)pcl.x, (GLfloat)pcl.y, (GLfloat)pcl.z,
				   (GLfloat)pcl.vol * vol_scale);
		}
		return 0;
	}
};

#endif