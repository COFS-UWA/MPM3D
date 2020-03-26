#include "QtPostProcessor_pcp.h"

#include <iostream>
#include "hdf5_utils.h"

#include "MainWindow.h"

herr_t get_time_history_name_func(
	hid_t loc_id,
	const char *name,
	const H5L_info_t *info,
	void *operator_data
	);

MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent),
	file_id(-1), md_id(-1), ths_id(-1),
	th_num(1), cur_th_id(0), cur_th_grp_id(-1)
{
	setObjectName(QString::fromUtf8("main_window"));
	resize(1000, 600);
	setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));

	// Menu bar
	menu_bar = new QMenuBar(this);
	menu_bar->setObjectName(QString::fromUtf8("menu_bar"));
	setMenuBar(menu_bar);

	// files menu
	menu_files = new QMenu(menu_bar);
	menu_files->setObjectName(QString::fromUtf8("menu_files"));
	menu_files->setTitle(QCoreApplication::translate("MainWindow", "Files", nullptr));
	menu_bar->addMenu(menu_files);

	action_open = new QAction(menu_files);
	action_open->setObjectName(QString::fromUtf8("action_open"));
	action_open->setText(tr("Open"));
	menu_files->addAction(action_open);
	
	action_close = new QAction(menu_files);
	action_close->setObjectName(QString::fromUtf8("action_close"));
	action_close->setText(tr("Close"));
	menu_files->addAction(action_close);

	// view menu
	menu_view = new QMenu(menu_bar);
	menu_view->setObjectName(QString::fromUtf8("menu_view"));
	menu_view->setTitle(tr("View"));
	menu_bar->addMenu(menu_view);

	action_show_bg_mesh = new QAction(menu_view);
	action_show_bg_mesh->setObjectName(QString::fromUtf8("action_close"));
	action_show_bg_mesh->setText(tr("Show mesh"));
	action_show_bg_mesh->setCheckable(true);
	menu_view->addAction(action_show_bg_mesh);

	// central widget
	central_widget = new QWidget(this);
	setCentralWidget(central_widget);

	// status bar
	status_bar = new QStatusBar(this);
	status_bar->setObjectName(QString::fromUtf8("status_bar"));
	setStatusBar(status_bar);

	// central widget
	main_layout = new QHBoxLayout(central_widget);
	main_layout->setObjectName(QString::fromUtf8("main_layout"));

	// left layout
	left_box_widget = new QWidget(central_widget);
	left_box_widget->setFixedWidth(350);
	main_layout->addWidget(left_box_widget);
	left_box = new QVBoxLayout;
	left_box_widget->setLayout(left_box);
	left_box->setObjectName(QString::fromUtf8("left_box"));
	left_box->setContentsMargins(0, 0, 0, 0);

	// right box
	right_box_widget = new QWidget(central_widget);
	main_layout->addWidget(right_box_widget);
	right_box = new QHBoxLayout;
	right_box_widget->setLayout(right_box);
	right_box->setObjectName(QString::fromUtf8("right_box"));
	right_box->setContentsMargins(0, 0, 0, 0);

	// left box
	output_combo = new QComboBox(this);
	output_combo->setObjectName(QString::fromUtf8("output_combo"));
	left_box->addWidget(output_combo);

	frame_list = new QListWidget(this);
	frame_list->setObjectName(QString::fromUtf8("frame_list"));
	left_box->addWidget(frame_list, 1);

	field_combo = new QComboBox(this);
	field_combo->setObjectName(QString::fromUtf8("field_combo"));
	left_box->addWidget(field_combo);

	model_info = new QLabel(this);
	model_info->setObjectName(QString::fromUtf8("model_info"));
	left_box->addWidget(model_info, 1);

	// right box
	openGL_widget = new GLWindow(this);
	openGL_widget->setObjectName(QString::fromUtf8("openGL_widget"));
	right_box->addWidget(openGL_widget, 1);

	view_angle_slider = new QSlider(this);
	view_angle_slider->setObjectName(QString::fromUtf8("view_angle_slider"));
	view_angle_slider->setMinimum(-180);
	view_angle_slider->setMaximum( 180);
	view_angle_slider->setValue(0);
	right_box->addWidget(view_angle_slider);

	// open hdf5 result file
	connect(action_open, SIGNAL(triggered()), this, SLOT(on_open_file()));
	// close hdf5 result file
	connect(action_close, SIGNAL(triggered()), this, SLOT(on_close_file()));
	// choose time history output
	connect(output_combo, SIGNAL(currentIndexChanged(int)), this, SLOT(on_select_time_history(int)));
	// 
	connect(action_show_bg_mesh, SIGNAL(triggered()), this, SLOT(on_show_bg_mesh()));
}

MainWindow::~MainWindow()
{
	on_close_file();
}

void MainWindow::on_open_file()
{
	on_close_file();

	//QString filename = QFileDialog::getOpenFileName(
	//	this,
	//	tr("Open result file"),
	//	QDir::currentPath(),
	//	tr("hdf5 files(*hdf5 *h5)")
	//	);
	
	QString filename("..\\..\\Asset\\t2d_chm_s_geostatic_hdf5.hdf5");

	if (filename == QString(""))
		return;

	std::string fn_str = filename.toStdString();
	file_id = H5Fopen(fn_str.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
		return;
	
	md_id = H5Gopen(file_id, "ModelData", H5P_DEFAULT);
	ths_id = H5Gopen(file_id, "TimeHistory", H5P_DEFAULT);

	// get name of all time history output
	QStringList th_name_list;
	th_name_list.append("Select time history...");
	H5Literate(
		ths_id,
		H5_INDEX_CRT_ORDER, //must H5Pset_link_creation_order() when created group 
		H5_ITER_NATIVE,
		nullptr,
		get_time_history_name_func,
		(void *)&th_name_list
		);
	output_combo->addItems(th_name_list);
	th_num = th_name_list.length();
	cur_th_id = 0;
}

herr_t get_time_history_name_func(
	hid_t loc_id,
	const char *name,
	const H5L_info_t *info,
	void *ext_data
	)
{
	QStringList &th_name_list = *(QStringList *)ext_data;

	H5O_info_t infobuf;
	H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
	if (infobuf.type == H5O_TYPE_GROUP)
	{
		th_name_list.push_back(name);
	}

	//std::cout << name << "\n";
	return 0;
}


void MainWindow::on_close_file()
{
	// clear output combobox
	th_num = 1;
	output_combo->clear();
	field_combo->clear();
	cur_th_id = 0;
	frame_list->clear();
	close_hdf5_group(cur_th_grp_id);

	close_hdf5_group(md_id);
	close_hdf5_group(ths_id);
	close_hdf5_file(file_id);
}


herr_t get_time_history_frame_func(
	hid_t loc_id,
	const char *name,
	const H5L_info_t *info,
	void *ext_data
	);

void MainWindow::on_select_time_history(int cur_id)
{
	//int new_th_id = output_combo->currentIndex();
	int new_th_id = cur_id;
	if (new_th_id != cur_th_id)
	{
		frame_list->clear();

		if (cur_th_grp_id >= 0)
		{
			H5Gclose(cur_th_grp_id);
			cur_th_grp_id = -1;
		}

		if (new_th_id != 0)
		{
			std::string th_name(output_combo->currentText().toStdString());
			cur_th_grp_id = H5Gopen(ths_id, th_name.c_str(), H5P_DEFAULT);

			// get time frame infos
			QStringList th_frame_list;
			H5Literate(
				cur_th_grp_id,
				H5_INDEX_CRT_ORDER, //must H5Pset_link_creation_order() when created group 
				H5_ITER_NATIVE,
				nullptr,
				get_time_history_frame_func,
				(void *)&th_frame_list
				);
			frame_list->addItems(th_frame_list);

			//// get type infos
			//QStringList mem_names;
			////cur_pcl_dset_id = H5Dopen(, "ParticleData", H5P_DEFAULT);
			//int mem_num = H5Tget_nmembers();
			//for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
			//{
			//	mem_names.append(H5Tget_member_name(, mem_id));
			//}
			//field_combo->addItems(mem_names);
		}

		cur_th_id = new_th_id;
	}
}

herr_t get_time_history_frame_func(
	hid_t loc_id,
	const char *name,
	const H5L_info_t *info,
	void *ext_data
	)
{
	QStringList &th_frame_list = *(QStringList *)ext_data;

	H5O_info_t infobuf;
	H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
	if (infobuf.type == H5O_TYPE_GROUP)
	{
		th_frame_list.push_back(name);
	}
	
	return 0;
}


void MainWindow::on_show_bg_mesh()
{
	if (action_show_bg_mesh->isChecked())
	{
		std::cout << "is checked\n";
	}
	std::cout << "show mesh\n";
}

