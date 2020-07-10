#ifndef __Qt_Scene_From_Model_h__
#define __Qt_Scene_From_Model_h__

#include <QOpenGLFunctions_3_3_Core>

class QtSceneFromModel
{
protected:
	QOpenGLFunctions_3_3_Core& gl;

public:
	QtSceneFromModel(QOpenGLFunctions_3_3_Core& _gl) : gl(_gl) {}
	~QtSceneFromModel() {}

	virtual int initialize(int wd, int ht) { return 0; }
	virtual void draw() {}
	virtual void resize(int wd, int ht) {}
};

#endif