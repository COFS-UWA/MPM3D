#ifndef __Qt_Scene_From_Hdf5_h__
#define __Qt_Scene_From_Hdf5_h__

#include <QOpenGLFunctions_3_3_Core>

class QtSceneFromHdf5
{
protected:
	QOpenGLFunctions_3_3_Core& gl;

public:
	QtSceneFromHdf5(QOpenGLFunctions_3_3_Core& _gl) : gl(_gl) {}
	~QtSceneFromHdf5() {}

	virtual size_t get_frame_num() { return 0; }
	virtual double get_frame_time(size_t frame_id) { return 0.0; }

	virtual int init_scene(int wd, int ht, size_t frame_id) { return 0; }
	virtual void update_scene(size_t frame_id) {}
	virtual void draw() {}
	virtual void resize(int wd, int ht) {}

};

#endif