#ifndef __Phong_Cube_Particle_System_h__
#define __Phong_Cube_Particle_System_h__

#include <QOpenGLWidget> // for opengl type

#include "ItemArray.hpp"
#include "Geometry.h"
#include "ValueToColor.h"

// particle system rendered with the Phong model
class PhongCubeParticleSystem
{
protected:
	typedef ValueToColor::Colorf Color;

	struct SingleVertData
	{
		GLfloat x, y, z;    // coordinates
		GLfloat r, g, b;    // color
		GLfloat nx, ny, nz; // normal vector
	};
	struct VertData
	{
		SingleVertData v0,  v1,  v2,  v3;  // face 1
		SingleVertData v4,  v5,  v6,  v7;  // face 2
		SingleVertData v8,  v9,  v10, v11; // face 3
		SingleVertData v12, v13, v14, v15; // face 4
		SingleVertData v16, v17, v18, v19; // face 5
		SingleVertData v20, v21, v22, v23; // face 6
		inline void init(GLfloat xc, GLfloat yc, GLfloat zc,
						 GLfloat vol, Color &c)
		{
			GLfloat half_len = pow(vol, 0.33333333f) * 0.5f;
			// face 1
			// v0
			v0.x = xc + half_len;
			v0.y = yc - half_len;
			v0.z = zc - half_len;
			v0.r = c.r;
			v0.g = c.g;
			v0.b = c.b;
			v0.nx = 1.0f;
			v0.ny = 0.0f;
			v0.nz = 0.0f;
			// v1
			v1.x = xc + half_len;
			v1.y = yc + half_len;
			v1.z = zc - half_len;
			v1.r = c.r;
			v1.g = c.g;
			v1.b = c.b;
			v1.nx = 1.0f;
			v1.ny = 0.0f;
			v1.nz = 0.0f;
			// v2
			v2.x = xc + half_len;
			v2.y = yc + half_len;
			v2.z = zc + half_len;
			v2.r = c.r;
			v2.g = c.g;
			v2.b = c.b;
			v2.nx = 1.0f;
			v2.ny = 0.0f;
			v2.nz = 0.0f;
			// v3
			v3.x = xc + half_len;
			v3.y = yc - half_len;
			v3.z = zc + half_len;
			v3.r = c.r;
			v3.g = c.g;
			v3.b = c.b;
			v3.nx = 1.0f;
			v3.ny = 0.0f;
			v3.nz = 0.0f;
			// face 2
			// v4
			v4.x = xc + half_len;
			v4.y = yc + half_len;
			v4.z = zc - half_len;
			v4.r = c.r;
			v4.g = c.g;
			v4.b = c.b;
			v4.nx = 0.0f;
			v4.ny = 1.0f;
			v4.nz = 0.0f;
			// v5
			v5.x = xc - half_len;
			v5.y = yc + half_len;
			v5.z = zc - half_len;
			v5.r = c.r;
			v5.g = c.g;
			v5.b = c.b;
			v5.nx = 0.0f;
			v5.ny = 1.0f;
			v5.nz = 0.0f;
			// v6
			v6.x = xc - half_len;
			v6.y = yc + half_len;
			v6.z = zc + half_len;
			v6.r = c.r;
			v6.g = c.g;
			v6.b = c.b;
			v6.nx = 0.0f;
			v6.ny = 1.0f;
			v6.nz = 0.0f;
			// v7
			v7.x = xc + half_len;
			v7.y = yc + half_len;
			v7.z = zc + half_len;
			v7.r = c.r;
			v7.g = c.g;
			v7.b = c.b;
			v7.nx = 0.0f;
			v7.ny = 1.0f;
			v7.nz = 0.0f;
			// face 3
			// v8
			v8.x = xc - half_len;
			v8.y = yc - half_len;
			v8.z = zc + half_len;
			v8.r = c.r;
			v8.g = c.g;
			v8.b = c.b;
			v8.nx = 0.0f;
			v8.ny = 0.0f;
			v8.nz = 1.0f;
			// v9
			v9.x = xc + half_len;
			v9.y = yc - half_len;
			v9.z = zc + half_len;
			v9.r = c.r;
			v9.g = c.g;
			v9.b = c.b;
			v9.nx = 0.0f;
			v9.ny = 0.0f;
			v9.nz = 1.0f;
			// v10
			v10.x = xc + half_len;
			v10.y = yc + half_len;
			v10.z = zc + half_len;
			v10.r = c.r;
			v10.g = c.g;
			v10.b = c.b;
			v10.nx = 0.0f;
			v10.ny = 0.0f;
			v10.nz = 1.0f;
			// v11
			v11.x = xc - half_len;
			v11.y = yc + half_len;
			v11.z = zc + half_len;
			v11.r = c.r;
			v11.g = c.g;
			v11.b = c.b;
			v11.nx = 0.0f;
			v11.ny = 0.0f;
			v11.nz = 1.0f;
			// face 4
			// v12
			v12.x = xc - half_len;
			v12.y = yc + half_len;
			v12.z = zc - half_len;
			v12.r = c.r;
			v12.g = c.g;
			v12.b = c.b;
			v12.nx = -1.0f;
			v12.ny = 0.0f;
			v12.nz = 0.0f;
			// v13
			v13.x = xc - half_len;
			v13.y = yc - half_len;
			v13.z = zc - half_len;
			v13.r = c.r;
			v13.g = c.g;
			v13.b = c.b;
			v13.nx = -1.0f;
			v13.ny = 0.0f;
			v13.nz = 0.0f;
			// v14
			v14.x = xc - half_len;
			v14.y = yc - half_len;
			v14.z = zc + half_len;
			v14.r = c.r;
			v14.g = c.g;
			v14.b = c.b;
			v14.nx = -1.0f;
			v14.ny = 0.0f;
			v14.nz = 0.0f;
			// v15
			v15.x = xc - half_len;
			v15.y = yc + half_len;
			v15.z = zc + half_len;
			v15.r = c.r;
			v15.g = c.g;
			v15.b = c.b;
			v15.nx = -1.0f;
			v15.ny = 0.0f;
			v15.nz = 0.0f;
			// face 5
			// v16
			v16.x = xc - half_len;
			v16.y = yc - half_len;
			v16.z = zc - half_len;
			v16.r = c.r;
			v16.g = c.g;
			v16.b = c.b;
			v16.nx = 0.0f;
			v16.ny = -1.0f;
			v16.nz = 0.0f;
			// v17
			v17.x = xc + half_len;
			v17.y = yc - half_len;
			v17.z = zc - half_len;
			v17.r = c.r;
			v17.g = c.g;
			v17.b = c.b;
			v17.nx = 0.0f;
			v17.ny = -1.0f;
			v17.nz = 0.0f;
			// v18
			v18.x = xc + half_len;
			v18.y = yc - half_len;
			v18.z = zc + half_len;
			v18.r = c.r;
			v18.g = c.g;
			v18.b = c.b;
			v18.nx = 0.0f;
			v18.ny = -1.0f;
			v18.nz = 0.0f;
			// v19
			v19.x = xc - half_len;
			v19.y = yc - half_len;
			v19.z = zc + half_len;
			v19.r = c.r;
			v19.g = c.g;
			v19.b = c.b;
			v19.nx = 0.0f;
			v19.ny = -1.0f;
			v19.nz = 0.0f;
			// face 6
			// v20
			v20.x = xc - half_len;
			v20.y = yc - half_len;
			v20.z = zc - half_len;
			v20.r = c.r;
			v20.g = c.g;
			v20.b = c.b;
			v20.nx = 0.0f;
			v20.ny = 0.0f;
			v20.nz = -1.0f;
			// v21
			v21.x = xc - half_len;
			v21.y = yc + half_len;
			v21.z = zc - half_len;
			v21.r = c.r;
			v21.g = c.g;
			v21.b = c.b;
			v21.nx = 0.0f;
			v21.ny = 0.0f;
			v21.nz = -1.0f;
			// v22
			v22.x = xc + half_len;
			v22.y = yc + half_len;
			v22.z = zc - half_len;
			v22.r = c.r;
			v22.g = c.g;
			v22.b = c.b;
			v22.nx = 0.0f;
			v22.ny = 0.0f;
			v22.nz = -1.0f;
			// v23
			v23.x = xc + half_len;
			v23.y = yc - half_len;
			v23.z = zc - half_len;
			v23.r = c.r;
			v23.g = c.g;
			v23.b = c.b;
			v23.nx = 0.0f;
			v23.ny = 0.0f;
			v23.nz = -1.0f;
		}
	};
	MemoryUtils::ItemArray<VertData> vert_data;
	GLsizei vert_data_size;

	struct SingleElemData { GLuint n1, n2, n3; };
	struct ElemData
	{
		SingleElemData e01, e02; // face 1
		SingleElemData e11, e12; // face 2
		SingleElemData e21, e22; // face 3
		SingleElemData e31, e32; // face 4
		SingleElemData e41, e42; // face 5
		SingleElemData e51, e52; // face 6
		inline void offset(GLuint off)
		{
			// face 1
			e01.n1 = 0 + off;
			e01.n2 = 1 + off;
			e01.n3 = 3 + off;
			e02.n1 = 1 + off;
			e02.n2 = 2 + off;
			e02.n3 = 3 + off;
			// face 2
			e11.n1 = 4 + off;
			e11.n2 = 5 + off;
			e11.n3 = 6 + off;
			e12.n1 = 4 + off;
			e12.n2 = 6 + off;
			e12.n3 = 7 + off;
			// face 3
			e21.n1 = 8 + off;
			e21.n2 = 9 + off;
			e21.n3 = 11 + off;
			e22.n1 = 9 + off;
			e22.n2 = 10 + off;
			e22.n3 = 11 + off;
			// face 4
			e31.n1 = 12 + off;
			e31.n2 = 13 + off;
			e31.n3 = 14 + off;
			e32.n1 = 12 + off;
			e32.n2 = 14 + off;
			e32.n3 = 15 + off;
			// face 5
			e41.n1 = 16 + off;
			e41.n2 = 17 + off;
			e41.n3 = 19 + off;
			e42.n1 = 17 + off;
			e42.n2 = 18 + off;
			e42.n3 = 19 + off;
			// face 6
			e51.n1 = 20 + off;
			e51.n2 = 21 + off;
			e51.n3 = 23 + off;
			e52.n1 = 21 + off;
			e52.n2 = 22 + off;
			e52.n3 = 23 + off;
		}
	};
	MemoryUtils::ItemArray<ElemData> elem_data; // face id
	GLsizei elem_data_size;

	size_t pcl_size, pcl_num;
	size_t x_off, y_off, z_off;
	size_t vol_off;
	float vol_scale;
	size_t fld_off;

public:
	PhongCubeParticleSystem() :
		vert_data_size(0), elem_data_size(0), pcl_num(0) {}
	~PhongCubeParticleSystem() { clear(); }
	inline void clear()
	{
		vert_data.clear();
		elem_data.clear();
	}

	inline GLvoid* get_vert_data() { return (GLvoid*)vert_data.get_mem(); }
	inline GLsizei get_vert_data_size() { return vert_data_size; }
	inline GLvoid* get_elem_data() { return (GLvoid*)elem_data.get_mem(); }
	inline GLsizei get_elem_data_size() { return elem_data_size; }

	// mono color version
	int init_data(char *pcls_data, size_t _pcl_size, size_t _pcl_num,
		size_t _x_off, size_t _y_off, size_t _z_off,
		size_t _vol_off, float _vol_scale, Color &_color);
	int update_data(char *pcls_data, Color &_color);

	// assume particle has x, y, z, vol
	template <typename Particle>
	int init_data(Particle *pcls, size_t _pcl_num,
				  float _vol_scale, Color &_color)
	{
		pcl_num = _pcl_num;
		vol_scale = _vol_scale;

		vert_data_size = pcl_num * sizeof(VertData);
		vert_data.reserve(pcl_num);
		VertData* verts = vert_data.get_mem();

		elem_data_size = pcl_num * sizeof(ElemData);
		elem_data.reserve(pcl_num);
		ElemData* elems = elem_data.get_mem();

		size_t cur_id_off = 0;
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			Particle& pcl = pcls[pcl_id];

			VertData& v = verts[pcl_id];
			v.init((GLfloat)pcl.x, (GLfloat)pcl.y, (GLfloat)pcl.z,
				   (GLfloat)pcl.get_vol() * vol_scale, _color);

			ElemData& e = elems[pcl_id];
			e.offset(cur_id_off);
			cur_id_off += 24;
		}

		return 0;
	}

	// Assumptions:
	// 1. Particle has members x, y, z, vol and field_value
	// 2. pcl_value_offset means offsetof(Particle, data_member)
	template <typename FieldType>
	int init_data(
		char *pcls_data, size_t _pcl_size, size_t _pcl_num,
		size_t _x_off, size_t _y_off, size_t _z_off,
		size_t _vol_off, float _vol_scale,
		size_t _fld_off, ValueToColor& v2c
		)
	{
		if (!pcls_data || _pcl_num == 0)
			return -1;

		pcl_size = _pcl_size;
		pcl_num = _pcl_num;
		x_off = _x_off;
		y_off = _y_off;
		z_off = _z_off;
		vol_off = _vol_off;
		vol_scale = _vol_scale;
		fld_off = _fld_off;

		vert_data_size = pcl_num * sizeof(VertData);
		vert_data.reserve(pcl_num);
		VertData *verts = vert_data.get_mem();

		elem_data_size = pcl_num * sizeof(ElemData);
		elem_data.reserve(pcl_num);
		ElemData *elems = elem_data.get_mem();

		char *pcl_data = pcls_data;
		double pcl_x, pcl_y, pcl_z, pcl_vol, pcl_fld;
		size_t cur_id_off = 0;
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			pcl_x = *(double*)(pcl_data + x_off);
			pcl_y = *(double*)(pcl_data + y_off);
			pcl_z = *(double*)(pcl_data + z_off);
			pcl_vol = *(double*)(pcl_data + vol_off);
			pcl_fld = double(*(FieldType*)(pcl_data + fld_off));
			VertData &v = verts[pcl_id];
			v.init((GLfloat)pcl_x, (GLfloat)pcl_y, (GLfloat)pcl_z,
				   (GLfloat)pcl_vol * vol_scale,
					v2c.get_color(pcl_fld));
			ElemData &e = elems[pcl_id];
			e.offset(cur_id_off);
			cur_id_off += 24;
			pcl_data += pcl_size;
		}

		return 0;
	}

	template <typename FieldType>
	int update_data(char *pcls_data, ValueToColor& v2c)
	{
		if (!pcls_data || pcl_num == 0)
			return -1;

		VertData* verts = vert_data.get_mem();
		char *pcl_data = pcls_data;
		double pcl_x, pcl_y, pcl_z, pcl_vol, pcl_fld;
		for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
		{
			pcl_x = *(double*)(pcl_data + x_off);
			pcl_y = *(double*)(pcl_data + y_off);
			pcl_z = *(double*)(pcl_data + z_off);
			pcl_vol = *(double*)(pcl_data + vol_off);
			pcl_fld = double(*(FieldType*)(pcl_data + fld_off));
			VertData& v = verts[pcl_id];
			v.init((GLfloat)pcl_x, (GLfloat)pcl_y, (GLfloat)pcl_z,
				   (GLfloat)pcl_vol * vol_scale,
					v2c.get_color(pcl_fld));
			pcl_data += pcl_size;
		}

		return 0;
	}

	int init_points(Point3D* pcls, size_t pcl_num,
		float pcl_vol, Color& color);
};

#endif