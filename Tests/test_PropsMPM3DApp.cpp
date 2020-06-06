#include "Tests_pcp.h"

#include "test_model_view.h"

#include "utils.h"
#include "ResultFile_hdf5.h"
#include "PospMPM3DApp.h"

void test_PospMPM3DApp(int argc, char **argv)
{
	PospMPM3DApp app(argc, argv, PospMPM3DApp::Animation);
	app.set_view_dir(10.0f, 30.0f);
	app.set_light_dir(10.0f, 30.0f);
	
	app.set_ani_time(5.0);
	app.set_gif_name("bar_com.gif");

	app.init_color_scale(-8.0e-6, 0.0, 
		ColorScaleExamples::get_color_scale(),
		ColorScaleExamples::get_color_num());

	ResultFile_hdf5 rf;
	rf.open("t3d_me_s_1d_compression.h5");
	int res = app.set_res_file(rf, "compression", "s33",
							   MPM3DModelView::BallShape);

	app.start();
}


// test the set_res_file_function
#include "Model_T3D_ME_s_hdf5_utilities.h"

void test_PospMPM3DApp_set_res_file()
{
	//using Model_T3D_ME_s_hdf5_utilities::ParticleData;
	//size_t psize = sizeof(ParticleData);
	//size_t x_off = offsetof(ParticleData, x);
	//size_t y_off = offsetof(ParticleData, y);
	//size_t z_off = offsetof(ParticleData, z);
	//size_t vol_off = offsetof(ParticleData, vol);
	//size_t vx_off = offsetof(ParticleData, vx);
}
