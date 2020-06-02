#include "ModelViewer_pcp.h"

#include "Model_hdf5_utilities.h"
#include "PospSingleFrameController.h"

PospSingleFrameController::PospSingleFrameController(MPM3DModelView& v) :
	MPM3DModelView::Controller(v),
	th_id(-1), frame_grp_id(-1), pcl_dt_id(-1) {}

PospSingleFrameController::~PospSingleFrameController() { close(); }

void PospSingleFrameController::close()
{
	ResultFile_hdf5& rf = *res_file;

	if (frame_grp_id >= 0)
	{
		rf.close_group(frame_grp_id);
		frame_grp_id = -1;
	}
	if (th_id >= 0)
	{
		rf.close_group(th_id);
		th_id = -1;
	}
	if (pcl_dt_id >= 0)
	{
		H5Tclose(pcl_dt_id);
		pcl_dt_id = -1;
	}
	res_file = nullptr;
}

int PospSingleFrameController::set_res_file(
	ResultFile_hdf5& rf,
	const char* th_na,
	size_t frame_id,
	const char* field_na
	)
{
	close();
	res_file = &rf;

	hid_t th_grp_id = rf.get_time_history_grp_id();

	// check if the time history exists
	if (!rf.has_group(th_grp_id, th_na))
	{
		close();
		return -1;
	}

	th_id = rf.open_group(th_grp_id, th_na);

	char frame[50];
	snprintf(frame, 50, "frame_%zu", frame_id);
	if (!rf.has_group(th_id, frame))
		return -2;

	frame_grp_id = rf.open_group(th_id, frame);
	hid_t pcl_grp_id = rf.open_group(frame_grp_id, "ParticleData");
	hid_t pcl_dset_id = rf.open_dataset(pcl_grp_id, "field");
	
	pcl_dt_id = H5Dget_type(pcl_dset_id);
	pcl_size = H5Tget_size(pcl_dt_id);
	int init_mem_num = 0;
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char* mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
		if (strcmp(mem_name, "x") == 0)
		{
			pcl_x_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, "y") == 0)
		{
			pcl_y_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, "z") == 0)
		{
			pcl_z_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, "vol") == 0)
		{
			pcl_vol_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			++init_mem_num;
		}
		else if (strcmp(mem_name, field_na) == 0)
		{
			field_name = field_na;
			pcl_fld_off = H5Tget_member_offset(pcl_dt_id, mem_id);
			pcl_fld_type = H5Tget_member_class(pcl_dt_id, mem_id);
			++init_mem_num;
		}
	}

	rf.close_dataset(pcl_dset_id);
	rf.close_group(pcl_grp_id);
	if (init_mem_num == 5)
		return 0;

	// no such data member
	return -3;
}


using namespace Model_hdf5_utilities;

int PospSingleFrameController::initialize_model_view_data()
{
	if (!res_file || !res_file->is_open())
		return -1;

	ResultFile_hdf5& rf = *res_file;

	// write model data to view
	hid_t md_data_grp_id = rf.get_model_data_grp_id();
	hid_t bg_mesh_id = rf.open_group(md_data_grp_id, "BackgroundMesh");
	size_t node_num, elem_num;
	// read nodes
	rf.read_attribute(bg_mesh_id, "node_num", node_num);
	hid_t node_dt_id = get_nd_3d_dt_id();
	Node3DData* nodes_data = new Node3DData[node_num];
	rf.read_dataset(bg_mesh_id, "NodeCoordinate", node_num, nodes_data, node_dt_id);
	H5Tclose(node_dt_id);
	// read elements
	rf.read_attribute(bg_mesh_id, "element_num", elem_num);
	hid_t elem_dt_id = get_ed_3d_dt_id();
	Elem3DData* elems_data = new Elem3DData[elem_num];
	rf.read_dataset(bg_mesh_id, "ElementTopology", elem_num, elems_data, elem_dt_id);
	H5Tclose(elem_dt_id);
	// write nodess and elems to window
	view->init_bg_mesh<Node3DData, Elem3DData>(
		nodes_data, node_num,
		elems_data, elem_num,
		QVector3D(0.5f, 0.5f, 0.5f)
		);
	delete[] nodes_data;
	delete[] elems_data;
	rf.close_group(bg_mesh_id);

	// write particle data to view
	hid_t frame_id = rf.open_group(th_id, "frame_0");
	hid_t pcl_data_id = rf.open_group(frame_id, "ParticleData");
	// pcl num
	rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	// pcl data
	char* pcls_data = new char[pcl_size * pcl_num];
	rf.read_dataset(pcl_data_id, "field", pcl_num, (void*)pcls_data, pcl_dt_id);
	rf.close_group(pcl_data_id);
	rf.close_group(frame_id);
	// write pcl data to view
	// support field data of different type
	//if (data_type == H5T_NATIVE_DOUBLE)
	view->init_multicolor_pcl_data<double>(
			pcls_data, pcl_size, pcl_num,
			pcl_x_off, pcl_y_off, pcl_z_off,
			pcl_vol_off, 0.125,
			pcl_fld_off, view->get_color_scale());
	delete[] pcls_data;

	return 0;
}
