#ifndef __Prep_Single_Frame_Controller_Base_h__
#define __Prep_Single_Frame_Controller_Base_h__

#include <QObject>

#include "ItemArray.hpp"
#include "ResultFile_hdf5.h"
#include "MPM3DModelView.h"

class PrepSingleFrameControllerBase : public QObject,
	public MPM3DModelView::Controller
{
	Q_OBJECT
protected:
	Point3D *points;
	size_t point_num;

public:
	PrepSingleFrameControllerBase(MPM3DModelView& v) :
		MPM3DModelView::Controller(v),
		points(nullptr), point_num(0) {}
	~PrepSingleFrameControllerBase() { close(); }
	void close() {}

	inline void set_points(Point3D *pts, size_t pt_num) { points = pts, point_num = pt_num; }

protected:
	int initialize_model_view_data() override
	{
		int res;
		if (points && point_num != 0)
		{
			QVector3D red(1.0f, 0.0f, 0.0f);
			if ((res = view->init_point_data(points, point_num, 10.0f, red)) != 0)
				return res;
		}
		return 0;
	}
};

#endif