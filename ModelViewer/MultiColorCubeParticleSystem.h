#ifndef __Multi_Color_Cube_Particle_System_h__
#define __Multi_Color_Cube_Particle_System_h__

#include <QOpenGLWidget> // for opengl type

#include "ItemArray.hpp"
#include "ValueToColor.h"

class MultiColorCubeParticleSystem
{
protected:
	typedef ValueToColor::Colorf Color;

	struct SingleVertData
	{
		GLfloat x, y, z;
		GLfloat r, g, b;
	};
	struct VertData
	{
		SingleVertData v1, v2, v3, v4, v5, v6, v7, v8;
		inline void init(GLfloat xc, GLfloat yc, GLfloat zc,
						 GLfloat vol, Color &c)
		{
			GLfloat half_len = pow(vol, 0.33333333f) * 0.5;
			// v1
			v1.x = xc + half_len;
			v1.y = yc + half_len;
			v1.z = zc + half_len;
			v1.r = c.r;
			v1.g = c.g;
			v1.b = c.b;
			// v2
			v2.x = xc - half_len;
			v2.y = yc + half_len;
			v2.z = zc + half_len;
			v2.r = c.r;
			v2.g = c.g;
			v2.b = c.b;
			// v3
			v3.x = xc - half_len;
			v3.y = yc - half_len;
			v3.z = zc + half_len;
			v3.r = c.r;
			v3.g = c.g;
			v3.b = c.b;
			// v4
			v4.x = xc + half_len;
			v4.y = yc - half_len;
			v4.z = zc + half_len;
			v4.r = c.r;
			v4.g = c.g;
			v4.b = c.b;
			// v5
			v5.x = xc + half_len;
			v5.y = yc + half_len;
			v5.z = zc - half_len;
			v5.r = c.r;
			v5.g = c.g;
			v5.b = c.b;
			// v6
			v6.x = xc - half_len;
			v6.y = yc + half_len;
			v6.z = zc - half_len;
			v6.r = c.r;
			v6.g = c.g;
			v6.b = c.b;
			// v7
			v7.x = xc - half_len;
			v7.y = yc - half_len;
			v7.z = zc - half_len;
			v7.r = c.r;
			v7.g = c.g;
			v7.b = c.b;
			// v8
			v8.x = xc + half_len;
			v8.y = yc - half_len;
			v8.z = zc - half_len;
			v8.r = c.r;
			v8.g = c.g;
			v8.b = c.b;
		}
	};
	MemoryUtils::ItemArray<VertData> vert_data;

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
	
public:
	MultiColorCubeParticleSystem() {}
	~MultiColorCubeParticleSystem() { clear(); }
	inline void clear()
	{
		vert_data.clear();
		elem_data.clear();
	}

	inline GLvoid *get_vert_data() { return (GLvoid *)vert_data.get_mem(); }
	inline GLvoid *get_elem_data() { return (GLvoid *)elem_data.get_mem(); }
	
	// Assumptions:
	// Particle has memebers x, y, z, vol
	// pcl value offset means offsetof(Particle, data_member)
	template <typename Particle>
	int init_data(Particle *pcls, size_t pcl_num,
		size_t pcl_value_offset, ValueToColor &v2c)
	{
		if (pcl_num == 0)
			return -1;

		vert_data.reserve(pcl_num);
		VertData *verts = vert_data.get_mem();
		elem_data.reserve(pcl_num);
		ElemData *elems = elem_data.get_mem();

		size_t cur_id_off = 0;
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			Particle &pcl = pcls[pcl_id];
			double value = *(double *)(((char*)(&pcl)) + pcl_value_offset);
			VertData &v = verts[pcl_id];
			v.init((GLfloat)pcl.x, (GLfloat)pcl.y, (GLfloat)pcl.z,
				   (GLfloat)pcl.vol, v2c.get_color(value));
			ElemData &e = elem[pcl_id];
			e.offset(cur_id_off);
			cur_id_off += 8;
		}

		return 0;
	}
};

#endif