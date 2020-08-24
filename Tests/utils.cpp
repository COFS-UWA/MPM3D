#include "Tests_pcp.h"

#include "utils.h"

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
