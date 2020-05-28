#include "Tests_pcp.h"

#include "test_model_view.h"

#include "ResultFile_hdf5.h"
#include "PospMPM3DApp.h"

void test_PospMPM3DApp(int argc, char **argv)
{
	PospMPM3DApp app(argc, argv, PospMPM3DApp::Animation);
	app.set_view_dir(10.0f, 30.0f);
	app.set_ani_time(5.0);

	app.set_gif_name("bar_com.gif");

	PospMPM3DApp::Colori colors[] = {
		{ 0,   0,   255 },
		{ 0,   93,  255 },
		{ 0,   185, 255 },
		{ 0,   255, 232 },
		{ 0,   255, 139 },
		{ 0,   255, 46 },
		{ 46,  255, 0 },
		{ 139, 255, 0 },
		{ 232, 255, 0 },
		{ 255, 185, 0 },
		{ 255, 93,  0 },
		{ 255, 0,   0 }
	};
	app.init_color_scale(-8.0e-6, 0.0, colors, sizeof(colors) / sizeof(colors[0]));

	ResultFile_hdf5 rf;
	rf.open("t3d_me_s_1d_compression.h5");
	int res = app.set_res_file(rf, "compression", "s33");

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