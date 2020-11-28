#include "ModelViewer_pcp.h"

#include "QtRigidCubeObject.h"

QtRigidCubeObject::QtRigidCubeObject(
    QOpenGLFunctions_3_3_Core& _gl) :
    gl(_gl), vao(0), vbo(0), vbo_index_num(0) {}

QtRigidCubeObject::~QtRigidCubeObject() {}

void QtRigidCubeObject::clear()
{
    if (vbo)
    {
        gl.glDeleteBuffers(1, &vbo);
        vbo = 0;
    }
    if (vao)
    {
        gl.glDeleteVertexArrays(1, &vao);
        vao = 0;
    }
    vbo_index_num = 0;
}

static const unsigned int cube_node_num = 8;

static const float cube_nodes[] = {
    -1.0f, -1.0f, -1.0f,
     1.0f, -1.0f, -1.0f,
     1.0f,  1.0f, -1.0f,
    -1.0f,  1.0f, -1.0f,
    -1.0f, -1.0f,  1.0f,
     1.0f, -1.0f,  1.0f,
     1.0f,  1.0f,  1.0f,
    -1.0f,  1.0f,  1.0f
};

static const unsigned int cube_elem_num = 12;

static const unsigned int cube_elems[] = {
    1, 2, 6,
    1, 6, 5,
    3, 4, 7,
    0, 4, 3,
    0, 1, 4,
    1, 5, 4,
    2, 3, 6,
    3, 7, 6,
    3, 1, 0,
    3, 2, 1,
    4, 6, 7,
    4, 5, 6
};

int QtRigidCubeObject::init(
    double _x,
    double _y,
    double _z,
    double _hx,
    double _hy,
    double _hz,
    QVector3D& c
    )
{
    color = c;

    gl.glGenVertexArrays(1, &vao);
    gl.glBindVertexArray(vao);

    double hhx = 0.5 * _hx;
    double hhy = 0.5 * _hy;
    double hhz = 0.5 * _hz;
    const float* pcn = cube_nodes;
    for (size_t n_id = 0; n_id < cube_node_num; ++n_id)
    {
        NodeCoord& nc = node_coords[n_id];
        nc.x = pcn[0] * hhx;
        nc.y = pcn[1] * hhy;
        nc.z = pcn[2] * hhz;
        pcn += 3;
    }

    vbo_index_num = 3 * cube_elem_num;

    const unsigned int* pce = cube_elems;
    NodeData* pnd = node_datas;
    GLfloat v12_x, v12_y, v12_z, v13_x, v13_y, v13_z;
    GLfloat v_norm, _nx, _ny, _nz;
    for (size_t e_id = 0; e_id < cube_elem_num; ++e_id)
    {
        NodeCoord& n1 = node_coords[pce[0]];
        NodeCoord& n2 = node_coords[pce[1]];
        NodeCoord& n3 = node_coords[pce[2]];
        v12_x = n2.x - n1.x;
        v12_y = n2.y - n1.y;
        v12_z = n2.z - n1.z;
        v13_x = n3.x - n1.x;
        v13_y = n3.y - n1.y;
        v13_z = n3.z - n1.z;
        _nx = v12_y * v13_z - v12_z * v13_y;
        _ny = v12_z * v13_x - v12_x * v13_z;
        _nz = v12_x * v13_y - v12_y * v13_x;
        v_norm = sqrt(_nx * _nx + _ny * _ny + _nz * _nz);
        _nx /= v_norm;
        _ny /= v_norm;
        _nz /= v_norm;

        NodeData& nd1 = pnd[0];
        nd1.type = 1;
        nd1.x = n1.x;
        nd1.y = n1.y;
        nd1.z = n1.z;
        nd1.nx = _nx;
        nd1.ny = _ny;
        nd1.nz = _nz;
        NodeData& nd2 = pnd[1];
        nd2.type = 1;
        nd2.x = n2.x;
        nd2.y = n2.y;
        nd2.z = n2.z;
        nd2.nx = _nx;
        nd2.ny = _ny;
        nd2.nz = _nz;
        NodeData& nd3 = pnd[2];
        nd3.type = 1;
        nd3.x = n3.x;
        nd3.y = n3.y;
        nd3.z = n3.z;
        nd3.nx = _nx;
        nd3.ny = _ny;
        nd3.nz = _nz;

        pce += 3;
        pnd += 3;
    }

    gl.glGenBuffers(1, &vbo);
    gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);
    gl.glBufferData(
        GL_ARRAY_BUFFER,
        vbo_index_num * sizeof(NodeData),
        node_datas,
        GL_STREAM_DRAW
    );

    // v_type
    gl.glVertexAttribIPointer(0,
        1, GL_UNSIGNED_INT,
        sizeof(NodeData),
        (GLvoid*)offsetof(NodeData, type)
    );
    gl.glEnableVertexAttribArray(0);
    // v_pos
    gl.glVertexAttribPointer(1,
        3, GL_FLOAT, GL_FALSE,
        sizeof(NodeData),
        (GLvoid*)offsetof(NodeData, x)
    );
    gl.glEnableVertexAttribArray(1);
    // v_norm
    gl.glVertexAttribPointer(2,
        3, GL_FLOAT, GL_FALSE,
        sizeof(NodeData),
        (GLvoid*)offsetof(NodeData, nx)
    );
    gl.glEnableVertexAttribArray(2);

    form_model_mat(Point3D(_x, _y, _z));
    return 0;
}

int QtRigidCubeObject::update(
    double x,
    double y,
    double z
)
{
    form_model_mat(Point3D(x, y, z));
    return 0;
}

void QtRigidCubeObject::draw(QOpenGLShaderProgram& shader)
{
    shader.bind();
    shader.setUniformValue("model_mat", model_mat);
    shader.setUniformValue("g_color", color);

    gl.glBindVertexArray(vao);

    //gl.glFrontFace(GL_CW);
    //gl.glCullFace(GL_BACK);
    gl.glDrawArrays(GL_TRIANGLES, 0, vbo_index_num);
}

inline void QtRigidCubeObject::form_model_mat(const Point3D& cen)
{
    float(*mm_data)[4] = reinterpret_cast<float(*)[4]>(model_mat.data());
    // opengl is column major
    // c++ is row major
    mm_data[0][0] = 1.0f;
    mm_data[0][1] = 0.0f;
    mm_data[0][2] = 0.0f;
    mm_data[0][3] = 0.0f;
    mm_data[1][0] = 0.0f;
    mm_data[1][1] = 1.0f;
    mm_data[1][2] = 0.0f;
    mm_data[1][3] = 0.0f;
    mm_data[2][0] = 0.0f;
    mm_data[2][1] = 0.0f;
    mm_data[2][2] = 1.0f;
    mm_data[2][3] = 0.0f;
    mm_data[3][0] = float(cen.x);
    mm_data[3][1] = float(cen.y);
    mm_data[3][2] = float(cen.z);
    mm_data[3][3] = 1.0f;
}
