#include "ModelViewer_pcp.h"

#include "ball_mesh_data.h"
#include "QtMultiColorBallGLObject.h"

QtMultiColorBallGLObject::QtMultiColorBallGLObject(QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), vao(0), vbo_cs(0), veo_cs(0), vbo_pts(0),
    c_elem_node_num(0), pt_num(0) {}

QtMultiColorBallGLObject::~QtMultiColorBallGLObject() { clear(); }

void QtMultiColorBallGLObject::clear()
{
    if (veo_cs)
    {
        gl.glDeleteBuffers(1, &veo_cs);
        veo_cs = 0;
    }
    if (vbo_cs)
    {
        gl.glDeleteBuffers(1, &vbo_cs);
        vbo_cs = 0;
    }
    if (vbo_pts)
    {
        gl.glDeleteBuffers(1, &vbo_pts);
        vbo_pts = 0;
    }
    if (vao)
    {
        gl.glDeleteVertexArrays(1, &vao);
        vao = 0;
    }
}

void QtMultiColorBallGLObject::draw(QOpenGLShaderProgram& shader)
{
    (void)shader;

    gl.glBindVertexArray(vao);
    gl.glDrawElementsInstanced(
        GL_TRIANGLES,
        c_elem_node_num,
        GL_UNSIGNED_INT,
        (GLvoid*)0,
        pt_num
        );
}

int QtMultiColorBallGLObject::init_gl_buffer(
    PointData* pds,
    size_t pd_num
    )
{
    clear();

    pt_num = pd_num;

    gl.glGenVertexArrays(1, &vao);
    if (vao == 0)
        return -1;
    gl.glBindVertexArray(vao);
    
    int res = init_ball_data();
    if (res)
        return res;

    gl.glGenBuffers(1, &vbo_pts);
    if (vbo_pts == 0)
        return -1;
    gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_pts);
    gl.glBufferData(GL_ARRAY_BUFFER,
        pd_num * sizeof(PointData),
        (GLvoid *)pds,
        GL_STREAM_DRAW
        );

    // pt_type
    gl.glVertexAttribIPointer(1,
        1, GL_UNSIGNED_INT,
        sizeof(PointData),
        (GLvoid*)offsetof(PointData, type)
        );
    gl.glEnableVertexAttribArray(1);
    gl.glVertexAttribDivisor(1, 1);
    // pt_pos
    gl.glVertexAttribPointer(2,
        3, GL_FLOAT, GL_FALSE,
        sizeof(PointData),
        (GLvoid *)offsetof(PointData, x)
        );
    gl.glEnableVertexAttribArray(2);
    gl.glVertexAttribDivisor(2, 1);
    // pt_radius
    gl.glVertexAttribPointer(3,
        1, GL_FLOAT, GL_FALSE,
        sizeof(PointData),
        (GLvoid*)offsetof(PointData, radius)
        );
    gl.glEnableVertexAttribArray(3);
    gl.glVertexAttribDivisor(3, 1);
    // point value
    gl.glVertexAttribPointer(4,
        1, GL_FLOAT, GL_FALSE,
        sizeof(PointData),
        (GLvoid*)offsetof(PointData, fld)
        );
    gl.glEnableVertexAttribArray(4);
    gl.glVertexAttribDivisor(4, 1);

    return 0;
}

int QtMultiColorBallGLObject::update_gl_buffer(
    PointData* pds,
    size_t pd_num
    )
{
    if (!vao || !vbo_pts)
        return -1;

    gl.glBindVertexArray(vao);

    gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_pts);
    gl.glBufferData(GL_ARRAY_BUFFER,
        pd_num * sizeof(PointData),
        (GLvoid *)pds,
        GL_STREAM_DRAW
        );

    return 0;
}

int QtMultiColorBallGLObject::init_ball_data()
{
    gl.glGenBuffers(1, &vbo_cs);
    if (vbo_cs == 0)
        return -1;
    gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_cs);
    gl.glBufferData(GL_ARRAY_BUFFER,
        sizeof(ball_nodes),
        ball_nodes,
        GL_STATIC_DRAW
        );

    gl.glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    gl.glEnableVertexAttribArray(0);

    gl.glGenBuffers(1, &veo_cs);
    if (veo_cs == 0)
        return -1;
    gl.glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, veo_cs);
    gl.glBufferData(GL_ELEMENT_ARRAY_BUFFER,
        sizeof(ball_elems),
        ball_elems,
        GL_STATIC_DRAW
        );

    c_elem_node_num = sizeof(ball_elems) / sizeof(ball_elems[0]);
    
    return 0;
}

int QtMultiColorBallGLObject::init(
    size_t pcl_num,
    double* pcl_x_data,
    double* pcl_y_data,
    double* pcl_z_data,
    double* pcl_vol_data,
    double *pcl_fld_data,
    float radius_scale
    )
{
    pt_data_mem.reserve(pcl_num);
    PointData* pt_data = pt_data_mem.get_mem();
    for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
    {
        PointData& pd = pt_data[pcl_id];
        pd.type = 1; // multi color
        pd.x = GLfloat(pcl_x_data[pcl_id]);
        pd.y = GLfloat(pcl_y_data[pcl_id]);
        pd.z = GLfloat(pcl_z_data[pcl_id]);
        pd.radius = GLfloat(pow(3.0 * pcl_vol_data[pcl_id] / (4.0 * 3.14159265359), 0.33333)) * radius_scale;
        pd.fld = GLfloat(pcl_fld_data[pcl_id]);
        int fefe = 0;
    }
    int res = init_gl_buffer(pt_data, pcl_num);
    return res;
}

int QtMultiColorBallGLObject::update(
    size_t pcl_num,
    double* pcl_x_data,
    double* pcl_y_data,
    double* pcl_z_data,
    double* pcl_vol_data,
    double *pcl_fld_data,
    float radius_scale
    )
{
    pt_data_mem.reserve(pcl_num);
    PointData* pt_data = pt_data_mem.get_mem();
    for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
    {
        PointData& pd = pt_data[pcl_id];
        pd.type = 1; // multi color
        pd.x = GLfloat(pcl_x_data[pcl_id]);
        pd.y = GLfloat(pcl_y_data[pcl_id]);
        pd.z = GLfloat(pcl_z_data[pcl_id]);
        pd.radius = GLfloat(pow(3.0 * pcl_vol_data[pcl_id] / (4.0 * 3.14159265359), 0.33333)) * radius_scale;
        pd.fld = GLfloat(pcl_fld_data[pcl_id]);
    }
    int res = update_gl_buffer(pt_data, pcl_num);
    return res;
}

