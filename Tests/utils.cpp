#include "Tests_pcp.h"

#include "utils.h"

ValueToColor::Colori ColorScaleExamples::abaqus_color_scale[] = {
	{ 0,   0,   255 },
	{ 0,   93,  255 },
	{ 0,   185, 255 },
	{ 0,   255, 232 },
	{ 0,   255, 139 },
	{ 0,   255, 46  },
	{ 46,  255, 0   },
	{ 139, 255, 0   },
	{ 232, 255, 0   },
	{ 255, 185, 0   },
	{ 255, 93,  0   },
	{ 255, 0,   0   }
};

size_t ColorScaleExamples::abaqus_color_scale_num 
	= sizeof(ColorScaleExamples::abaqus_color_scale)
	/ sizeof(ColorScaleExamples::abaqus_color_scale[0]);


// ME Model
void init_vx_bcs_display(Model_T3D_ME_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_ME_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
	{
		Model_T3D_ME_s::Node& n = nodes[md.vxs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vy_bcs_display(Model_T3D_ME_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_ME_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
	{
		Model_T3D_ME_s::Node& n = nodes[md.vys[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vz_bcs_display(Model_T3D_ME_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_ME_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vz_num; ++v_id)
	{
		Model_T3D_ME_s::Node& n = nodes[md.vzs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

// CHM Model
void init_vsx_bcs_display(Model_T3D_CHM_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_CHM_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
	{
		Model_T3D_CHM_s::Node& n = nodes[md.vsxs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vsy_bcs_display(Model_T3D_CHM_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_CHM_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
	{
		Model_T3D_CHM_s::Node& n = nodes[md.vsys[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vsz_bcs_display(Model_T3D_CHM_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_CHM_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vsz_num; ++v_id)
	{
		Model_T3D_CHM_s::Node& n = nodes[md.vszs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vfx_bcs_display(Model_T3D_CHM_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_CHM_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
	{
		Model_T3D_CHM_s::Node& n = nodes[md.vfxs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vfy_bcs_display(Model_T3D_CHM_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_CHM_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
	{
		Model_T3D_CHM_s::Node& n = nodes[md.vfys[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vfz_bcs_display(Model_T3D_CHM_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_T3D_CHM_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vfz_num; ++v_id)
	{
		Model_T3D_CHM_s::Node& n = nodes[md.vfzs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

// FEM ME
void init_tz_face_bcs_display(Model_FEM_T3D_ME_s& md, Point3DArray& ptlist)
{
	typedef Model_FEM_T3D_ME_s::Element Element;
	typedef Model_FEM_T3D_ME_s::Node Node;
	Point3D pt;
	Element* elems = md.get_elems();
	Node* nodes = md.get_nodes();
	for (size_t t_id = 0; t_id < md.tz_num; ++t_id)
	{
		Element& e = elems[md.tzs[t_id].elem_id];
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		Node& n4 = nodes[e.n4];
		switch (md.tzs[t_id].face_id)
		{
		case 0:
			pt.x = n1.x;
			pt.y = n1.y;
			pt.z = n1.z;
			ptlist.add(pt);
			pt.x = n2.x;
			pt.y = n2.y;
			pt.z = n2.z;
			ptlist.add(pt);
			pt.x = n3.x;
			pt.y = n3.y;
			pt.z = n3.z;
			ptlist.add(pt);
			break;
		case 1:
			pt.x = n1.x;
			pt.y = n1.y;
			pt.z = n1.z;
			ptlist.add(pt);
			pt.x = n4.x;
			pt.y = n4.y;
			pt.z = n4.z;
			ptlist.add(pt);
			pt.x = n2.x;
			pt.y = n2.y;
			pt.z = n2.z;
			ptlist.add(pt);
			break;
		case 2:
			pt.x = n1.x;
			pt.y = n1.y;
			pt.z = n1.z;
			ptlist.add(pt);
			pt.x = n3.x;
			pt.y = n3.y;
			pt.z = n3.z;
			ptlist.add(pt);
			pt.x = n4.x;
			pt.y = n4.y;
			pt.z = n4.z;
			ptlist.add(pt);
			break;
		case 3:
			pt.x = n2.x;
			pt.y = n2.y;
			pt.z = n2.z;
			ptlist.add(pt);
			pt.x = n4.x;
			pt.y = n4.y;
			pt.z = n4.z;
			ptlist.add(pt);
			pt.x = n3.x;
			pt.y = n3.y;
			pt.z = n3.z;
			break;
		default:
			break;
		}
	}


}

void init_vx_bcs_display(Model_FEM_T3D_ME_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_FEM_T3D_ME_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
	{
		Model_FEM_T3D_ME_s::Node& n = nodes[md.vxs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vy_bcs_display(Model_FEM_T3D_ME_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_FEM_T3D_ME_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
	{
		Model_FEM_T3D_ME_s::Node& n = nodes[md.vys[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}

void init_vz_bcs_display(Model_FEM_T3D_ME_s& md, Point3DArray& ptlist)
{
	Point3D pt;
	Model_FEM_T3D_ME_s::Node* nodes = md.get_nodes();
	for (size_t v_id = 0; v_id < md.vz_num; ++v_id)
	{
		Model_FEM_T3D_ME_s::Node& n = nodes[md.vzs[v_id].node_id];
		pt.x = n.x;
		pt.y = n.y;
		pt.z = n.z;
		ptlist.add(pt);
	}
}
