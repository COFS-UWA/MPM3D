#ifndef __Posp_Single_Frame_Controller_h__
#define __Posp_Single_Frame_Controller_h__

#include <QObject>

#include "ItemArray.hpp"
#include "ResultFile_hdf5.h"
#include "MPM3DModelView.h"

// load result from hdf5
// display single frame of the hdf5 file
class PospSingleFrameController : public QObject,
	public MPM3DModelView::Controller
{
	Q_OBJECT

protected:
	ResultFile_hdf5* res_file;
	hid_t th_id;
	hid_t frame_grp_id;
	size_t pcl_num;
	hid_t pcl_dt_id;
	size_t pcl_size;
	size_t pcl_x_off;
	size_t pcl_y_off;
	size_t pcl_z_off;
	size_t pcl_vol_off;
	std::string field_name;
	size_t pcl_fld_off;
	hid_t pcl_fld_type;
	
	// add png output here...

public:
	PospSingleFrameController(MPM3DModelView& v);
	~PospSingleFrameController();
	void close();

	int set_res_file(
		ResultFile_hdf5& rf,
		const char* th_na,
		size_t frame_id,
		const char* field_na
		);

protected:
	int initialize_model_view_data() override;
};

#endif