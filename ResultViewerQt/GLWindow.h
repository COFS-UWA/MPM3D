#ifndef __GL_Window_H__
#define __GL_Window_H__

#include <QMouseEvent>
#include <QWheelEvent>
#include <QKeyEvent>
#include <QElapsedTimer>
#include <QTimer>
#include <QOpenGLWidget>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core> 

class GLWindow : public QOpenGLWidget,
	public QOpenGLFunctions_3_3_Core
{
	Q_OBJECT
public:
	explicit GLWindow(QWidget *parent = Q_NULLPTR);
	~GLWindow();
	void initializeGL();
	void paintGL();
	void resizeGL(int w, int h);
	
private slots:
	void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
	void mouseReleaseEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
	void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
	void wheelEvent(QWheelEvent *event) Q_DECL_OVERRIDE;
	//void on_rotate(void);

	void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
	void keyReleaseEvent(QKeyEvent *event) Q_DECL_OVERRIDE;

/* =============== model ============== */
protected:
	QOpenGLShaderProgram shader;
	GLuint vao, vbo, ebo;

/* ======== for frame updating ======== */
protected:
	QElapsedTimer frame_time;
	int last_frame_time;
	QTimer frame_timer;

	bool should_render_frame;
	float min_frame_inv; // ms
	// frame interval adjustment
	float adj_frame_inv; // ms
	
protected slots:
	void on_frame_timer_timeout();

/* ====== for movement of camera ====== */
protected: // for movement of camera
	QMatrix4x4 view_mat;
	QMatrix4x4 proj_mat;

	bool is_roaming;
	QCursor cursor;

	// camera pos
	float cam_x, cam_y, cam_z;

	// camera direction
	float cam_yaw, cam_pitch, cam_roll;
	bool need_update_loc_dir;
	QVector3D loc_x_dir, loc_y_dir, loc_z_dir;

	float ang_sensitivity;
	int velocity_index;
	float velocity0;
	float velocity;

	QElapsedTimer cam_time;
	int time_start_moving_left;
	int time_start_moving_right;
	int time_start_moving_up;
	int time_start_moving_down;

protected:
	inline void update_loc_dir()
	{
		if (need_update_loc_dir)
		{
			float cos_pitch = cos(cam_pitch);
			loc_x_dir.setX(cos_pitch * cos(cam_yaw));
			loc_x_dir.setY(cos_pitch * sin(cam_yaw));
			loc_x_dir.setZ(sin(cam_pitch));
			loc_y_dir.setX(-loc_x_dir.y());
			loc_y_dir.setY(loc_x_dir.x());
			loc_y_dir.setZ(0.0f);
			need_update_loc_dir = false;
		}
	}
};

#endif