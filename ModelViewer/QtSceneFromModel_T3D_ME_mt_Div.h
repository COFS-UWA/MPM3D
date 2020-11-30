#ifndef __Qt_Scene_From_Model_T3D_ME_mt_Div_h__
#define __Qt_Scene_From_Model_T3D_ME_mt_Div_h__

#include <QOpenGLShaderProgram>

#include "DivisionSet.h"
#include "QtSceneFromModel_T3D_ME_mt.h"

template <class DivisionSet>
class QtSceneFromModel_T3D_ME_mt_Div :
	public QtSceneFromModel_T3D_ME_mt
{
protected:
	DivisionSet div_set;

public:
	QtSceneFromModel_T3D_ME_mt_Div(QOpenGLFunctions_3_3_Core& _gl) :
		QtSceneFromModel_T3D_ME_mt(_gl) {}
	~QtSceneFromModel_T3D_ME_mt_Div() {}

	inline DivisionSet &get_div_set() noexcept { return div_set; }

	int initialize(int wd, int ht)
	{
		gl.glEnable(GL_DEPTH_TEST);
		width = wd; height = ht;

		init_shaders();

		// init bg_mesh
		QVector3D gray(0.5f, 0.5f, 0.5f);
		bg_mesh_obj.init_from_elements<
			Model_T3D_ME_mt::Position,
			Model_T3D_ME_mt::ElemNodeIndex,
			DivisionSet>(
			model->get_node_pos(),
			model->get_node_num(),
			model->get_elem_node_index(),
			model->get_elem_num(),
			gray,
			div_set
		);

		// init pcls
		size_t pcl_num = model->get_pcl_num();
		double *pcl_x = new double[pcl_num * 4];
		double* pcl_y = pcl_x + pcl_num;
		double* pcl_z = pcl_y + pcl_num;
		double* pcl_vol = pcl_z + pcl_num;
		auto *p_pos = model->get_pcl_pos();
		auto *p_m = model->get_pcl_m();
		auto *p_d = model->get_pcl_density0();
		for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		{
			pcl_x[p_id] = p_pos[p_id].x;
			pcl_y[p_id] = p_pos[p_id].y;
			pcl_z[p_id] = p_pos[p_id].z;
			pcl_vol[p_id] = p_m[p_id] / p_d[p_id];
		}
		QVector3D moccasin(1.0f, 0.8941f, 0.7098f);
		pcls_obj.init<DivisionSet>(pcl_x, pcl_y, pcl_z, pcl_vol,
					  pcl_num, moccasin, 0.5f, div_set);
		delete[] pcl_x;

		// init pts
		init_pts_buffer();
		// init rigid object
		init_rigid_objects_buffer();
		return 0;
	}
};

#endif