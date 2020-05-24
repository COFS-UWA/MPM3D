#ifndef __Time_History_Model_View_h__
#define __Time_History_Model_View_h__

#include "TimeHistory.h"

#include "Geometry.h"
#include "PrepMPM3DApp.h"

class TimeHistory_ModelView : public TimeHistory
{
protected:
	Point3D* points;
	size_t point_num;

	MPM3DModelView* view;

public:
	TimeHistory_ModelView(Point3D *pts, size_t pt_num, MPM3DModelView *v,
		const char* _name, const char* _type, TimeHistoryFunc _output_func) :
		TimeHistory(_name, _type, _output_func),
		points(pts), point_num(pt_num), view(v)
	{
		view->th_mv = this;
	}

	~TimeHistory_ModelView() {}

	inline void set_points(Point3D* pts, size_t pt_num) { points = pts, point_num = pt_num; }
	
	inline void set_view(MPM3DModelView& v)
	{
		view = &v;
		view->th_mv = this;
	}

	inline void set_view(PrepMPM3DApp& v)
	{
		view = &v.get_model_view();
		view->th_mv = this;
	}

	virtual int initialize_model_view_data();
};

#endif