#ifndef __Main_Window_H__
#define __Main_Window_H__

#include <QApplication>
#include <QMainWindow>

#include <QMenuBar>
#include <QMenu>
#include <QAction>

#include <QStatusBar>

#include <QHBoxLayout>
#include <QVBoxLayout>

#include <QComboBox>
#include <QListWidget>
#include <QLabel>

#include <QSlider>

#include "hdf5.h"

#include "GLWindow.h"

class MainWindow : public QMainWindow
{
	Q_OBJECT
public:
	// menu bar
	QMenuBar *menu_bar;
	// files menu
	QMenu *menu_files;
	QAction *action_open;
	QAction *action_close;
	// view menu
	QMenu *menu_view;
	QAction *action_show_bg_mesh;

	// status bar
	QStatusBar *status_bar;

	// central widget
	QWidget *central_widget;
	QHBoxLayout *main_layout;

	QWidget *left_box_widget;
	QVBoxLayout *left_box;

	QWidget *right_box_widget;
	QHBoxLayout *right_box;
	
	QComboBox *output_combo;
	QListWidget *frame_list;
	QComboBox *field_combo;
	QLabel *model_info;

	GLWindow *openGL_widget;
	QSlider *view_angle_slider;

public:
	explicit MainWindow(QWidget *parent = Q_NULLPTR);
	~MainWindow();

protected:
	hid_t file_id;
	hid_t md_id;
	hid_t ths_id;

private slots:
	void on_open_file();
	void on_close_file();

	void on_select_time_history(int cur_id);

	void on_show_bg_mesh();

protected:
	size_t th_num, cur_th_id;
	hid_t cur_th_grp_id;
	hid_t cur_pcl_dset_id;
	hid_t cur_pcl_dset_type_id;
	
};

#endif