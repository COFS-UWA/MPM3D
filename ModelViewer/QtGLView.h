#ifndef __Qt_GL_View_h__
#define __Qt_GL_View_h__

#include <QMatrix4x4>
#include <QOpenGLWidget>
#include <QOpenGLFunctions_3_3_Core>

class QtController;

class QtGLView : public QOpenGLWidget,
	public QOpenGLFunctions_3_3_Core
{
	Q_OBJECT
protected:
	QtController *controller;

private:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int wd, int ht) override;

	bool fully_loaded;
	unsigned int resize_num;
	typedef void (QtGLView::* DrawFunc)();
	DrawFunc draw_func;
	typedef void (QtGLView::* ResizeFunc)(int wd, int ht);
	ResizeFunc resize_func;

	void draw_view_ini();
	void resize_view_ini(int wd, int ht);

public:
	explicit QtGLView(QWidget* parent = Q_NULLPTR);
	virtual ~QtGLView();

	inline void set_controller(QtController &con) { controller = &con; }
	inline bool is_fully_loaded() { return fully_loaded; }

	virtual void draw_view();
	virtual void resize_view(int wd, int ht);
};

#endif