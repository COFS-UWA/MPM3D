#include "ModelViewer_pcp.h"

#include "QtRigidRectObject.h"

QtRigidRectObject::QtRigidRectObject(
	QOpenGLFunctions_3_3_Core& _gl) :
	gl(_gl), vao(0), vbo(0) {}

QtRigidRectObject::~QtRigidRectObject() { clear(); }

void QtRigidRectObject::clear()
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
}

int QtRigidRectObject::init(
    double x, double y, double ang,
    double hx, double hy,
    QVector3D& c, GLfloat line_wd
	)
{
	color = c;
	line_width = line_wd;

	gl.glGenVertexArrays(1, &vao);
	gl.glBindVertexArray(vao);

    double sin_ang = sin(ang);
    double cos_ang = cos(ang);
    hhx = 0.5 * hx;
    hhy = 0.5 * hy;
    NodeData& nd0 = node_datas[0];
    nd0.type = 0;
    nd0.x = cos_ang * hhx + -sin_ang * hhy + x;
    nd0.y = sin_ang * hhx +  cos_ang * hhy + y;
    NodeData& nd1 = node_datas[1];
    nd1.type = 0;
    nd1.x = cos_ang * -hhx + -sin_ang * hhy + x;
    nd1.y = sin_ang * -hhx +  cos_ang * hhy + y;
    NodeData& nd2 = node_datas[2];
    nd2.type = 0;
    nd2.x = cos_ang * -hhx + -sin_ang * -hhy + x;
    nd2.y = sin_ang * -hhx +  cos_ang * -hhy + y;
    NodeData& nd3 = node_datas[3];
    nd3.type = 0;
    nd3.x = cos_ang * hhx + -sin_ang * -hhy + x;
    nd3.y = sin_ang * hhx +  cos_ang * -hhy + y;

	gl.glGenBuffers(1, &vbo);
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);
    gl.glBufferData(
        GL_ARRAY_BUFFER,
        4 * sizeof(NodeData),
        node_datas,
        GL_STREAM_DRAW
        );

    // v_type
    gl.glVertexAttribIPointer(0,
        1, GL_INT,
        sizeof(NodeData),
        (GLvoid*)offsetof(NodeData, type)
        );
    gl.glEnableVertexAttribArray(0);
    // v_pos
    gl.glVertexAttribPointer(1,
        2, GL_FLOAT, GL_FALSE,
        sizeof(NodeData),
        (GLvoid*)offsetof(NodeData, x)
        );
    gl.glEnableVertexAttribArray(1);
    // v_color (not used)
    gl.glVertexAttribPointer(2,
        3, GL_FLOAT, GL_FALSE,
        0, (GLvoid*)0
        );
    gl.glEnableVertexAttribArray(2);

	return 0;
}

int QtRigidRectObject::update(double x, double y, double ang)
{
    if (vao == 0)
        return 0;

    double sin_ang = sin(ang);
    double cos_ang = cos(ang);
    NodeData& nd0 = node_datas[0];
    nd0.type = 0;
    nd0.x = cos_ang * hhx + -sin_ang * hhy + x;
    nd0.y = sin_ang * hhx + cos_ang * hhy + y;
    NodeData& nd1 = node_datas[1];
    nd1.type = 0;
    nd1.x = cos_ang * -hhx + -sin_ang * hhy + x;
    nd1.y = sin_ang * -hhx + cos_ang * hhy + y;
    NodeData& nd2 = node_datas[2];
    nd2.type = 0;
    nd2.x = cos_ang * -hhx + -sin_ang * -hhy + x;
    nd2.y = sin_ang * -hhx + cos_ang * -hhy + y;
    NodeData& nd3 = node_datas[3];
    nd3.type = 0;
    nd3.x = cos_ang * hhx + -sin_ang * -hhy + x;
    nd3.y = sin_ang * hhx + cos_ang * -hhy + y;

    gl.glBindVertexArray(vao);
    gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);
    gl.glBufferSubData(GL_ARRAY_BUFFER,
        0, 4 * sizeof(NodeData),
        (GLvoid *)node_datas
        );

    return 0;
}

void QtRigidRectObject::draw(QOpenGLShaderProgram& shader)
{
	shader.bind();
	shader.setUniformValue("g_color", color);

	gl.glLineWidth(line_width);
	gl.glBindVertexArray(vao);
    gl.glDrawArrays(GL_LINE_LOOP, 0, 4);
}