#ifndef __Model_To_Viewer_3D_hpp__
#define __Model_To_Viewer_3D_hpp__

#include "ModelToViewerBase.h"
#include "ModelWindow.h"

template <typename ModelType>
class ModelToViewer3D : public ModelToViewerBase
{
public:
	ModelToViewer3D() {}
	ModelToViewer3D(Model &_md, ModelWindow &_md_win) : 
		ModelToViewerBase(_md, _md_win.get_gl_win()) {}
	
	int init_win_data() override
	{
		ModelType *md = static_cast<ModelType *>(model);
		// init mesh data
		md_win->glGenVertexArrays(1, &md_win->vao_mh);
		md_win->glBindVertexArray(md_win->vao_mh);
		// node data
		size_t node_num = md->get_node_num();
		typename ModelType::Node *nodes = md->get_nodes();
		GLfloat *node_data = new GLfloat[node_num * 3];
		GLfloat *pnd = node_data;
		for (size_t n_id = 0; n_id < node_num; ++n_id)
		{
			typename ModelType::Node &n = nodes[n_id];
			*pnd = n.x;
			++pnd;
			*pnd = n.y;
			++pnd;
			*pnd = n.z;
			++pnd;
		}
		md_win->glGenBuffers(1, &md_win->vbo_mh);
		md_win->glBindBuffer(GL_ARRAY_BUFFER, md_win->vbo_mh);
		md_win->glBufferData(
			GL_ARRAY_BUFFER,
			node_num * 3 * sizeof(GLfloat),
			node_data,
			GL_STATIC_DRAW
			);
		delete[] node_data;
		// coordinates
		md_win->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
		md_win->glEnableVertexAttribArray(0);
		// element edge data
		size_t edge_num = md->get_edge_num();
		md_win->line_num = edge_num;
		typename ModelType::Edge *edges = md->get_edges();
		GLuint *edge_data = new GLuint[edge_num * 2];
		GLuint *ped = edge_data;
		for (size_t e_id = 0; e_id < edge_num; ++e_id)
		{
			typename ModelType::Edge &e = edges[e_id];
			*ped = e.n1;
			++ped;
			*ped = e.n2;
			++ped;
		}
		md_win->glGenBuffers(1, &md_win->ebo_mh);
		md_win->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, md_win->ebo_mh);
		md_win->glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			edge_num * 2 * sizeof(GLuint),
			edge_data,
			GL_STATIC_DRAW
			);
		delete[] edge_data;
		// set view angle
		Point3D mh_cen = md->get_centre();
		md_win->mh_centre.setX(mh_cen.x);
		md_win->mh_centre.setY(mh_cen.y);
		md_win->mh_centre.setZ(mh_cen.z);
		Cube bbox = md->get_bounding_box();
		float bbox_dx = bbox.xu - bbox.xl;
		float bbox_dy = bbox.yu - bbox.yl;
		float bbox_dz = bbox.zu - bbox.zl;
		md_win->mh_radius = sqrt(bbox_dx*bbox_dx + bbox_dy*bbox_dy + bbox_dz*bbox_dz) * 0.5;
		md_win->mh_radius *= 1.05; // expand a little bit
		return 0;
	}
};

#endif