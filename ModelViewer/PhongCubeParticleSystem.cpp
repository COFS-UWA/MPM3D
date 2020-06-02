#include "ModelViewer_pcp.h"

#include "PhongCubeParticleSystem.h"

// Mono color particles
int PhongCubeParticleSystem::init_data(
	char *pcls_data, size_t _pcl_size, size_t _pcl_num,
	size_t _x_off, size_t _y_off, size_t _z_off,
	size_t _vol_off, float _vol_scale,
	Color& _color
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

	vert_data_size = pcl_num * sizeof(VertData);
	vert_data.reserve(pcl_num);
	VertData* verts = vert_data.get_mem();

	elem_data_size = pcl_num * sizeof(ElemData);
	elem_data.reserve(pcl_num);
	ElemData* elems = elem_data.get_mem();

	char* pcl_data = pcls_data;
	double pcl_x, pcl_y, pcl_z, pcl_vol;
	size_t cur_id_off = 0;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		pcl_x = *(double*)(pcl_data + x_off);
		pcl_y = *(double*)(pcl_data + y_off);
		pcl_z = *(double*)(pcl_data + z_off);
		pcl_vol = *(double*)(pcl_data + vol_off);
		VertData& v = verts[pcl_id];
		v.init((GLfloat)pcl_x, (GLfloat)pcl_y, (GLfloat)pcl_z,
			   (GLfloat)pcl_vol * vol_scale, _color);
		ElemData& e = elems[pcl_id];
		e.offset(cur_id_off);
		cur_id_off += 24;
		pcl_data += pcl_size;
	}

	return 0;
}

int PhongCubeParticleSystem::update_data(
	char* pcls_data,
	Color& _color
	)
{
	if (!pcls_data || pcl_num == 0)
		return -1;

	VertData* verts = vert_data.get_mem();
	char* pcl_data = pcls_data;
	double pcl_x, pcl_y, pcl_z, pcl_vol;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		pcl_x = *(double*)(pcl_data + x_off);
		pcl_y = *(double*)(pcl_data + y_off);
		pcl_z = *(double*)(pcl_data + z_off);
		pcl_vol = *(double*)(pcl_data + vol_off);
		VertData& v = verts[pcl_id];
		v.init((GLfloat)pcl_x, (GLfloat)pcl_y, (GLfloat)pcl_z,
			   (GLfloat)pcl_vol * vol_scale, _color);
	}
	return 0;
}
