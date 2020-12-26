#ifndef __Qt_Scene_From_Hdf5_T3D_ME_mt_Div_h__
#define __Qt_Scene_From_Hdf5_T3D_ME_mt_Div_h__

#include <QOpenGLShaderProgram>

#include "DivisionSet.h"
#include "QtSceneFromHdf5_T3D_ME_mt.h"

template <class DivisionSet>
class QtSceneFromHdf5_T3D_ME_mt_Div :
	public QtSceneFromHdf5_T3D_ME_mt
{
protected:
	DivisionSet div_set;

public:
	explicit QtSceneFromHdf5_T3D_ME_mt_Div(QOpenGLFunctions_3_3_Core& _gl) :
		QtSceneFromHdf5_T3D_ME_mt(_gl) {}
	~QtSceneFromHdf5_T3D_ME_mt_Div() {}

	inline DivisionSet& get_div_set() noexcept { return div_set; }
	
	int init_scene(int wd, int ht, size_t frame_id) override;
	void update_scene(size_t frame_id) override;
};

template <class DivisionSet>
int QtSceneFromHdf5_T3D_ME_mt_Div<DivisionSet>
	::init_scene(
	int wd, int ht,
	size_t frame_id
	)
{
	using namespace Model_hdf5_utilities;

	if (!res_file || !res_file->is_open())
		return -1;

	win_wd = wd; win_ht = ht;

	int res;
	ResultFile_hdf5& rf = *res_file;

	// write model data to view
	hid_t md_data_grp_id = rf.get_model_data_grp_id();
	hid_t bg_mesh_id = rf.open_group(md_data_grp_id, "BackgroundMesh");
	// read nodes
	size_t node_num;
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	hid_t node_dt_id = get_nd_3d_dt_id();
	Node3DData* nodes_data = new Node3DData[node_num];
	rf.read_dataset(bg_mesh_id, "NodeData", node_num, nodes_data, node_dt_id);
	H5Tclose(node_dt_id);
	// read elements
	size_t elem_num;
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	hid_t elem_dt_id = get_ed_3d_dt_id();
	Elem3DData* elems_data = new Elem3DData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementData", elem_num, elems_data, elem_dt_id);
	H5Tclose(elem_dt_id);
	// get bounding box
	Cube mh_bbox;
	mh_bbox.xl = nodes_data[0].x;
	mh_bbox.xu = mh_bbox.xl;
	mh_bbox.yl = nodes_data[0].y;
	mh_bbox.yu = mh_bbox.yl;
	mh_bbox.zl = nodes_data[0].z;
	mh_bbox.zu = mh_bbox.zl;
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		Node3DData& n = nodes_data[n_id];
		if (mh_bbox.xl > n.x)
			mh_bbox.xl = n.x;
		if (mh_bbox.xu < n.x)
			mh_bbox.xu = n.x;
		if (mh_bbox.yl > n.y)
			mh_bbox.yl = n.y;
		if (mh_bbox.yu < n.y)
			mh_bbox.yu = n.y;
		if (mh_bbox.zl > n.z)
			mh_bbox.zl = n.z;
		if (mh_bbox.zu < n.z)
			mh_bbox.zu = n.z;
	}
	// init bg_mesh buffer
	QVector3D gray(0.5f, 0.5f, 0.5f);
	res = bg_mesh_obj.init_from_elements<
		Node3DData, Elem3DData, DivisionSet>(
			nodes_data, node_num,
			elems_data, elem_num,
			gray, div_set
			);
	delete[] nodes_data;
	delete[] elems_data;
	rf.close_group(bg_mesh_id);
	if (res)
		return res;

	// init particle data
	need_mat_model_data = pfld->need_mat_model_data();
	res = data_loader.load_frame_data(frame_id, need_mat_model_data);
	if (res)
		return res;

	size_t pcl_num = data_loader.get_pcl_num();
	if (pcl_num)
	{
		pcl_fld_mem.reserve(pcl_num * 5);
		// x coord
		double* pcl_x_data = pcl_fld_mem.get_mem();
		x_fld.extract_pcl_fld_data(pcl_x_data);
		// y coord
		double* pcl_y_data = pcl_x_data + pcl_num;
		y_fld.extract_pcl_fld_data(pcl_y_data);
		// z coord
		double* pcl_z_data = pcl_y_data + pcl_num;
		z_fld.extract_pcl_fld_data(pcl_z_data);
		// vol
		double* pcl_vol_data = pcl_z_data + pcl_num;
		vol_fld.extract_pcl_fld_data(pcl_vol_data);
		// pcl fld
		double* pcl_fld_data = pcl_vol_data + pcl_num;
		pfld->extract_pcl_fld_data(pcl_fld_data);
		pcls_obj.init<DivisionSet>(
			pcl_num,
			pcl_x_data,
			pcl_y_data,
			pcl_z_data,
			pcl_vol_data,
			pcl_fld_data,
			0.5f, div_set
			);
	}

	init_rigid_objects_buffer(mh_bbox, rf);
	if (!init_color_map_texture())
		return -2;
	init_shaders(mh_bbox, wd, ht);
	return 0;
}

template <class DivisionSet>
void QtSceneFromHdf5_T3D_ME_mt_Div<DivisionSet>
	::update_scene(size_t frame_id)
{
	ResultFile_hdf5& rf = *res_file;

	data_loader.load_frame_data(frame_id, need_mat_model_data);
	size_t pcl_num = data_loader.get_pcl_num();
	if (pcl_num)
	{
		pcl_fld_mem.reserve(pcl_num * 5);
		// x coord
		double* pcl_x_data = pcl_fld_mem.get_mem();
		x_fld.extract_pcl_fld_data(pcl_x_data);
		// y coord
		double* pcl_y_data = pcl_x_data + pcl_num;
		y_fld.extract_pcl_fld_data(pcl_y_data);
		// z coord
		double* pcl_z_data = pcl_y_data + pcl_num;
		z_fld.extract_pcl_fld_data(pcl_z_data);
		// vol
		double* pcl_vol_data = pcl_z_data + pcl_num;
		vol_fld.extract_pcl_fld_data(pcl_vol_data);
		// pcl fld
		double* pcl_fld_data = pcl_vol_data + pcl_num;
		pfld->extract_pcl_fld_data(pcl_fld_data);
		pcls_obj.update<DivisionSet>(
			pcl_num,
			pcl_x_data,
			pcl_y_data,
			pcl_z_data,
			pcl_vol_data,
			pcl_fld_data,
			0.5f,
			div_set
		);
		//pcls_obj.update(
		//	pcl_num,
		//	pcl_x_data,
		//	pcl_y_data,
		//	pcl_z_data,
		//	pcl_vol_data,
		//	pcl_fld_data,
		//	0.5
		//	);
	}

	char frame_name[50];
	snprintf(frame_name, 50, "frame_%zu", frame_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	update_rigid_objects_buffer(frame_grp_id, rf);
	rf.close_group(frame_grp_id);
}

#endif