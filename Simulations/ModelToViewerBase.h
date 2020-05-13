#ifndef __Model_To_Viewer_3D_h__
#define __Model_To_Viewer_3D_h__

class Model;
class GLModelWindow;

class ModelToViewerBase
{
protected:
	Model *model;
	GLModelWindow *md_win;

public:
	ModelToViewerBase() : model(nullptr), md_win(nullptr) {}
	ModelToViewerBase(Model &_md, GLModelWindow &_md_win) : model(&_md), md_win(&_md_win) {}
	inline void set_model(Model &_md) { model = &_md; }
	inline void set_win(GLModelWindow &_md_win) { md_win = &_md_win; }
	
	virtual int init_win_data() = 0;
};

#endif