#include "Tests_pcp.h"

#include "test_model_view.h"

#include "utils.h"
#include "ResultFile_hdf5.h"
#include "PospMPM3DApp.h"

void test_PospSingleFrame_display(int argc, char **argv)
{
	PospMPM3DApp app(argc, argv, PospMPM3DApp::SingleFrame);
	app.set_view_dir(10.0f, 10.0f);
	app.set_light_dir(10.0f, 30.0f);

	app.init_color_scale(-8.0e-6, 0.0,
		ColorScaleExamples::get_color_scale(),
		ColorScaleExamples::get_color_num());

	ResultFile_hdf5 rf;
	rf.open("t3d_me_s_1d_compression.h5");
	int res = app.set_res_file(rf, "compression", 0, "s33");

	app.start();
}