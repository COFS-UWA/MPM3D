#include "ModelViewer_pcp.h"

#include "QtRigidCylinderObject.h"

QtRigidCylinderObject::QtRigidCylinderObject(
    QOpenGLFunctions_3_3_Core& _gl) :
    gl(_gl), vao(0), vbo(0), vbo_index_num(0) {}

QtRigidCylinderObject::~QtRigidCylinderObject() {}

void QtRigidCylinderObject::clear()
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

static const unsigned int cylinder_node_num = 242;

static const float cylinder_nodes[] = {
    0.000000f, 0.000000f, 0.500000f, // node 0
    0.000000f, 0.000000f, -0.500000f, // node 1
    1.000000f, 0.000000f, 0.500000f, // node 2
    0.998630f, 0.052336f, 0.500000f, // node 3
    0.994522f, 0.104528f, 0.500000f, // node 4
    0.987688f, 0.156434f, 0.500000f, // node 5
    0.978148f, 0.207912f, 0.500000f, // node 6
    0.965926f, 0.258819f, 0.500000f, // node 7
    0.951057f, 0.309017f, 0.500000f, // node 8
    0.933580f, 0.358368f, 0.500000f, // node 9
    0.913545f, 0.406737f, 0.500000f, // node 10
    0.891007f, 0.453990f, 0.500000f, // node 11
    0.866025f, 0.500000f, 0.500000f, // node 12
    0.838671f, 0.544639f, 0.500000f, // node 13
    0.809017f, 0.587785f, 0.500000f, // node 14
    0.777146f, 0.629320f, 0.500000f, // node 15
    0.743145f, 0.669131f, 0.500000f, // node 16
    0.707107f, 0.707107f, 0.500000f, // node 17
    0.669131f, 0.743145f, 0.500000f, // node 18
    0.629320f, 0.777146f, 0.500000f, // node 19
    0.587785f, 0.809017f, 0.500000f, // node 20
    0.544639f, 0.838671f, 0.500000f, // node 21
    0.500000f, 0.866025f, 0.500000f, // node 22
    0.453990f, 0.891007f, 0.500000f, // node 23
    0.406737f, 0.913545f, 0.500000f, // node 24
    0.358368f, 0.933580f, 0.500000f, // node 25
    0.309017f, 0.951057f, 0.500000f, // node 26
    0.258819f, 0.965926f, 0.500000f, // node 27
    0.207912f, 0.978148f, 0.500000f, // node 28
    0.156434f, 0.987688f, 0.500000f, // node 29
    0.104528f, 0.994522f, 0.500000f, // node 30
    0.052336f, 0.998630f, 0.500000f, // node 31
    0.000000f, 1.000000f, 0.500000f, // node 32
    -0.052336f, 0.998630f, 0.500000f, // node 33
    -0.104528f, 0.994522f, 0.500000f, // node 34
    -0.156434f, 0.987688f, 0.500000f, // node 35
    -0.207912f, 0.978148f, 0.500000f, // node 36
    -0.258819f, 0.965926f, 0.500000f, // node 37
    -0.309017f, 0.951057f, 0.500000f, // node 38
    -0.358368f, 0.933580f, 0.500000f, // node 39
    -0.406737f, 0.913545f, 0.500000f, // node 40
    -0.453990f, 0.891007f, 0.500000f, // node 41
    -0.500000f, 0.866025f, 0.500000f, // node 42
    -0.544639f, 0.838671f, 0.500000f, // node 43
    -0.587785f, 0.809017f, 0.500000f, // node 44
    -0.629320f, 0.777146f, 0.500000f, // node 45
    -0.669131f, 0.743145f, 0.500000f, // node 46
    -0.707107f, 0.707107f, 0.500000f, // node 47
    -0.743145f, 0.669131f, 0.500000f, // node 48
    -0.777146f, 0.629320f, 0.500000f, // node 49
    -0.809017f, 0.587785f, 0.500000f, // node 50
    -0.838671f, 0.544639f, 0.500000f, // node 51
    -0.866025f, 0.500000f, 0.500000f, // node 52
    -0.891007f, 0.453990f, 0.500000f, // node 53
    -0.913545f, 0.406737f, 0.500000f, // node 54
    -0.933580f, 0.358368f, 0.500000f, // node 55
    -0.951057f, 0.309017f, 0.500000f, // node 56
    -0.965926f, 0.258819f, 0.500000f, // node 57
    -0.978148f, 0.207912f, 0.500000f, // node 58
    -0.987688f, 0.156434f, 0.500000f, // node 59
    -0.994522f, 0.104528f, 0.500000f, // node 60
    -0.998630f, 0.052336f, 0.500000f, // node 61
    -1.000000f, 0.000000f, 0.500000f, // node 62
    -0.998630f, -0.052336f, 0.500000f, // node 63
    -0.994522f, -0.104528f, 0.500000f, // node 64
    -0.987688f, -0.156434f, 0.500000f, // node 65
    -0.978148f, -0.207912f, 0.500000f, // node 66
    -0.965926f, -0.258819f, 0.500000f, // node 67
    -0.951057f, -0.309017f, 0.500000f, // node 68
    -0.933580f, -0.358368f, 0.500000f, // node 69
    -0.913545f, -0.406737f, 0.500000f, // node 70
    -0.891007f, -0.453990f, 0.500000f, // node 71
    -0.866025f, -0.500000f, 0.500000f, // node 72
    -0.838671f, -0.544639f, 0.500000f, // node 73
    -0.809017f, -0.587785f, 0.500000f, // node 74
    -0.777146f, -0.629320f, 0.500000f, // node 75
    -0.743145f, -0.669131f, 0.500000f, // node 76
    -0.707107f, -0.707107f, 0.500000f, // node 77
    -0.669131f, -0.743145f, 0.500000f, // node 78
    -0.629320f, -0.777146f, 0.500000f, // node 79
    -0.587785f, -0.809017f, 0.500000f, // node 80
    -0.544639f, -0.838671f, 0.500000f, // node 81
    -0.500000f, -0.866025f, 0.500000f, // node 82
    -0.453990f, -0.891007f, 0.500000f, // node 83
    -0.406737f, -0.913545f, 0.500000f, // node 84
    -0.358368f, -0.933580f, 0.500000f, // node 85
    -0.309017f, -0.951057f, 0.500000f, // node 86
    -0.258819f, -0.965926f, 0.500000f, // node 87
    -0.207912f, -0.978148f, 0.500000f, // node 88
    -0.156434f, -0.987688f, 0.500000f, // node 89
    -0.104528f, -0.994522f, 0.500000f, // node 90
    -0.052336f, -0.998630f, 0.500000f, // node 91
    -0.000000f, -1.000000f, 0.500000f, // node 92
    0.052336f, -0.998630f, 0.500000f, // node 93
    0.104528f, -0.994522f, 0.500000f, // node 94
    0.156434f, -0.987688f, 0.500000f, // node 95
    0.207912f, -0.978148f, 0.500000f, // node 96
    0.258819f, -0.965926f, 0.500000f, // node 97
    0.309017f, -0.951057f, 0.500000f, // node 98
    0.358368f, -0.933580f, 0.500000f, // node 99
    0.406737f, -0.913545f, 0.500000f, // node 100
    0.453990f, -0.891007f, 0.500000f, // node 101
    0.500000f, -0.866025f, 0.500000f, // node 102
    0.544639f, -0.838671f, 0.500000f, // node 103
    0.587785f, -0.809017f, 0.500000f, // node 104
    0.629320f, -0.777146f, 0.500000f, // node 105
    0.669131f, -0.743145f, 0.500000f, // node 106
    0.707107f, -0.707107f, 0.500000f, // node 107
    0.743145f, -0.669131f, 0.500000f, // node 108
    0.777146f, -0.629320f, 0.500000f, // node 109
    0.809017f, -0.587785f, 0.500000f, // node 110
    0.838671f, -0.544639f, 0.500000f, // node 111
    0.866025f, -0.500000f, 0.500000f, // node 112
    0.891007f, -0.453990f, 0.500000f, // node 113
    0.913545f, -0.406737f, 0.500000f, // node 114
    0.933580f, -0.358368f, 0.500000f, // node 115
    0.951057f, -0.309017f, 0.500000f, // node 116
    0.965926f, -0.258819f, 0.500000f, // node 117
    0.978148f, -0.207912f, 0.500000f, // node 118
    0.987688f, -0.156434f, 0.500000f, // node 119
    0.994522f, -0.104528f, 0.500000f, // node 120
    0.998630f, -0.052336f, 0.500000f, // node 121
    1.000000f, 0.000000f, -0.500000f, // node 122
    0.998630f, 0.052336f, -0.500000f, // node 123
    0.994522f, 0.104528f, -0.500000f, // node 124
    0.987688f, 0.156434f, -0.500000f, // node 125
    0.978148f, 0.207912f, -0.500000f, // node 126
    0.965926f, 0.258819f, -0.500000f, // node 127
    0.951057f, 0.309017f, -0.500000f, // node 128
    0.933580f, 0.358368f, -0.500000f, // node 129
    0.913545f, 0.406737f, -0.500000f, // node 130
    0.891007f, 0.453990f, -0.500000f, // node 131
    0.866025f, 0.500000f, -0.500000f, // node 132
    0.838671f, 0.544639f, -0.500000f, // node 133
    0.809017f, 0.587785f, -0.500000f, // node 134
    0.777146f, 0.629320f, -0.500000f, // node 135
    0.743145f, 0.669131f, -0.500000f, // node 136
    0.707107f, 0.707107f, -0.500000f, // node 137
    0.669131f, 0.743145f, -0.500000f, // node 138
    0.629320f, 0.777146f, -0.500000f, // node 139
    0.587785f, 0.809017f, -0.500000f, // node 140
    0.544639f, 0.838671f, -0.500000f, // node 141
    0.500000f, 0.866025f, -0.500000f, // node 142
    0.453990f, 0.891007f, -0.500000f, // node 143
    0.406737f, 0.913545f, -0.500000f, // node 144
    0.358368f, 0.933580f, -0.500000f, // node 145
    0.309017f, 0.951057f, -0.500000f, // node 146
    0.258819f, 0.965926f, -0.500000f, // node 147
    0.207912f, 0.978148f, -0.500000f, // node 148
    0.156434f, 0.987688f, -0.500000f, // node 149
    0.104528f, 0.994522f, -0.500000f, // node 150
    0.052336f, 0.998630f, -0.500000f, // node 151
    0.000000f, 1.000000f, -0.500000f, // node 152
    -0.052336f, 0.998630f, -0.500000f, // node 153
    -0.104528f, 0.994522f, -0.500000f, // node 154
    -0.156434f, 0.987688f, -0.500000f, // node 155
    -0.207912f, 0.978148f, -0.500000f, // node 156
    -0.258819f, 0.965926f, -0.500000f, // node 157
    -0.309017f, 0.951057f, -0.500000f, // node 158
    -0.358368f, 0.933580f, -0.500000f, // node 159
    -0.406737f, 0.913545f, -0.500000f, // node 160
    -0.453990f, 0.891007f, -0.500000f, // node 161
    -0.500000f, 0.866025f, -0.500000f, // node 162
    -0.544639f, 0.838671f, -0.500000f, // node 163
    -0.587785f, 0.809017f, -0.500000f, // node 164
    -0.629320f, 0.777146f, -0.500000f, // node 165
    -0.669131f, 0.743145f, -0.500000f, // node 166
    -0.707107f, 0.707107f, -0.500000f, // node 167
    -0.743145f, 0.669131f, -0.500000f, // node 168
    -0.777146f, 0.629320f, -0.500000f, // node 169
    -0.809017f, 0.587785f, -0.500000f, // node 170
    -0.838671f, 0.544639f, -0.500000f, // node 171
    -0.866025f, 0.500000f, -0.500000f, // node 172
    -0.891007f, 0.453990f, -0.500000f, // node 173
    -0.913545f, 0.406737f, -0.500000f, // node 174
    -0.933580f, 0.358368f, -0.500000f, // node 175
    -0.951057f, 0.309017f, -0.500000f, // node 176
    -0.965926f, 0.258819f, -0.500000f, // node 177
    -0.978148f, 0.207912f, -0.500000f, // node 178
    -0.987688f, 0.156434f, -0.500000f, // node 179
    -0.994522f, 0.104528f, -0.500000f, // node 180
    -0.998630f, 0.052336f, -0.500000f, // node 181
    -1.000000f, 0.000000f, -0.500000f, // node 182
    -0.998630f, -0.052336f, -0.500000f, // node 183
    -0.994522f, -0.104528f, -0.500000f, // node 184
    -0.987688f, -0.156434f, -0.500000f, // node 185
    -0.978148f, -0.207912f, -0.500000f, // node 186
    -0.965926f, -0.258819f, -0.500000f, // node 187
    -0.951057f, -0.309017f, -0.500000f, // node 188
    -0.933580f, -0.358368f, -0.500000f, // node 189
    -0.913545f, -0.406737f, -0.500000f, // node 190
    -0.891007f, -0.453990f, -0.500000f, // node 191
    -0.866025f, -0.500000f, -0.500000f, // node 192
    -0.838671f, -0.544639f, -0.500000f, // node 193
    -0.809017f, -0.587785f, -0.500000f, // node 194
    -0.777146f, -0.629320f, -0.500000f, // node 195
    -0.743145f, -0.669131f, -0.500000f, // node 196
    -0.707107f, -0.707107f, -0.500000f, // node 197
    -0.669131f, -0.743145f, -0.500000f, // node 198
    -0.629320f, -0.777146f, -0.500000f, // node 199
    -0.587785f, -0.809017f, -0.500000f, // node 200
    -0.544639f, -0.838671f, -0.500000f, // node 201
    -0.500000f, -0.866025f, -0.500000f, // node 202
    -0.453990f, -0.891007f, -0.500000f, // node 203
    -0.406737f, -0.913545f, -0.500000f, // node 204
    -0.358368f, -0.933580f, -0.500000f, // node 205
    -0.309017f, -0.951057f, -0.500000f, // node 206
    -0.258819f, -0.965926f, -0.500000f, // node 207
    -0.207912f, -0.978148f, -0.500000f, // node 208
    -0.156434f, -0.987688f, -0.500000f, // node 209
    -0.104528f, -0.994522f, -0.500000f, // node 210
    -0.052336f, -0.998630f, -0.500000f, // node 211
    -0.000000f, -1.000000f, -0.500000f, // node 212
    0.052336f, -0.998630f, -0.500000f, // node 213
    0.104528f, -0.994522f, -0.500000f, // node 214
    0.156434f, -0.987688f, -0.500000f, // node 215
    0.207912f, -0.978148f, -0.500000f, // node 216
    0.258819f, -0.965926f, -0.500000f, // node 217
    0.309017f, -0.951057f, -0.500000f, // node 218
    0.358368f, -0.933580f, -0.500000f, // node 219
    0.406737f, -0.913545f, -0.500000f, // node 220
    0.453990f, -0.891007f, -0.500000f, // node 221
    0.500000f, -0.866025f, -0.500000f, // node 222
    0.544639f, -0.838671f, -0.500000f, // node 223
    0.587785f, -0.809017f, -0.500000f, // node 224
    0.629320f, -0.777146f, -0.500000f, // node 225
    0.669131f, -0.743145f, -0.500000f, // node 226
    0.707107f, -0.707107f, -0.500000f, // node 227
    0.743145f, -0.669131f, -0.500000f, // node 228
    0.777146f, -0.629320f, -0.500000f, // node 229
    0.809017f, -0.587785f, -0.500000f, // node 230
    0.838671f, -0.544639f, -0.500000f, // node 231
    0.866025f, -0.500000f, -0.500000f, // node 232
    0.891007f, -0.453990f, -0.500000f, // node 233
    0.913545f, -0.406737f, -0.500000f, // node 234
    0.933580f, -0.358368f, -0.500000f, // node 235
    0.951057f, -0.309017f, -0.500000f, // node 236
    0.965926f, -0.258819f, -0.500000f, // node 237
    0.978148f, -0.207912f, -0.500000f, // node 238
    0.987688f, -0.156434f, -0.500000f, // node 239
    0.994522f, -0.104528f, -0.500000f, // node 240
    0.998630f, -0.052336f, -0.500000f // node 241
};

static const unsigned int cylinder_elem_num = 480;

static const unsigned int cylinder_elems[] = {
    0, 2, 3, // elem 0
    0, 3, 4, // elem 1
    0, 4, 5, // elem 2
    0, 5, 6, // elem 3
    0, 6, 7, // elem 4
    0, 7, 8, // elem 5
    0, 8, 9, // elem 6
    0, 9, 10, // elem 7
    0, 10, 11, // elem 8
    0, 11, 12, // elem 9
    0, 12, 13, // elem 10
    0, 13, 14, // elem 11
    0, 14, 15, // elem 12
    0, 15, 16, // elem 13
    0, 16, 17, // elem 14
    0, 17, 18, // elem 15
    0, 18, 19, // elem 16
    0, 19, 20, // elem 17
    0, 20, 21, // elem 18
    0, 21, 22, // elem 19
    0, 22, 23, // elem 20
    0, 23, 24, // elem 21
    0, 24, 25, // elem 22
    0, 25, 26, // elem 23
    0, 26, 27, // elem 24
    0, 27, 28, // elem 25
    0, 28, 29, // elem 26
    0, 29, 30, // elem 27
    0, 30, 31, // elem 28
    0, 31, 32, // elem 29
    0, 32, 33, // elem 30
    0, 33, 34, // elem 31
    0, 34, 35, // elem 32
    0, 35, 36, // elem 33
    0, 36, 37, // elem 34
    0, 37, 38, // elem 35
    0, 38, 39, // elem 36
    0, 39, 40, // elem 37
    0, 40, 41, // elem 38
    0, 41, 42, // elem 39
    0, 42, 43, // elem 40
    0, 43, 44, // elem 41
    0, 44, 45, // elem 42
    0, 45, 46, // elem 43
    0, 46, 47, // elem 44
    0, 47, 48, // elem 45
    0, 48, 49, // elem 46
    0, 49, 50, // elem 47
    0, 50, 51, // elem 48
    0, 51, 52, // elem 49
    0, 52, 53, // elem 50
    0, 53, 54, // elem 51
    0, 54, 55, // elem 52
    0, 55, 56, // elem 53
    0, 56, 57, // elem 54
    0, 57, 58, // elem 55
    0, 58, 59, // elem 56
    0, 59, 60, // elem 57
    0, 60, 61, // elem 58
    0, 61, 62, // elem 59
    0, 62, 63, // elem 60
    0, 63, 64, // elem 61
    0, 64, 65, // elem 62
    0, 65, 66, // elem 63
    0, 66, 67, // elem 64
    0, 67, 68, // elem 65
    0, 68, 69, // elem 66
    0, 69, 70, // elem 67
    0, 70, 71, // elem 68
    0, 71, 72, // elem 69
    0, 72, 73, // elem 70
    0, 73, 74, // elem 71
    0, 74, 75, // elem 72
    0, 75, 76, // elem 73
    0, 76, 77, // elem 74
    0, 77, 78, // elem 75
    0, 78, 79, // elem 76
    0, 79, 80, // elem 77
    0, 80, 81, // elem 78
    0, 81, 82, // elem 79
    0, 82, 83, // elem 80
    0, 83, 84, // elem 81
    0, 84, 85, // elem 82
    0, 85, 86, // elem 83
    0, 86, 87, // elem 84
    0, 87, 88, // elem 85
    0, 88, 89, // elem 86
    0, 89, 90, // elem 87
    0, 90, 91, // elem 88
    0, 91, 92, // elem 89
    0, 92, 93, // elem 90
    0, 93, 94, // elem 91
    0, 94, 95, // elem 92
    0, 95, 96, // elem 93
    0, 96, 97, // elem 94
    0, 97, 98, // elem 95
    0, 98, 99, // elem 96
    0, 99, 100, // elem 97
    0, 100, 101, // elem 98
    0, 101, 102, // elem 99
    0, 102, 103, // elem 100
    0, 103, 104, // elem 101
    0, 104, 105, // elem 102
    0, 105, 106, // elem 103
    0, 106, 107, // elem 104
    0, 107, 108, // elem 105
    0, 108, 109, // elem 106
    0, 109, 110, // elem 107
    0, 110, 111, // elem 108
    0, 111, 112, // elem 109
    0, 112, 113, // elem 110
    0, 113, 114, // elem 111
    0, 114, 115, // elem 112
    0, 115, 116, // elem 113
    0, 116, 117, // elem 114
    0, 117, 118, // elem 115
    0, 118, 119, // elem 116
    0, 119, 120, // elem 117
    0, 120, 121, // elem 118
    0, 121, 2, // elem 119
    1, 123, 122, // elem 120
    1, 124, 123, // elem 121
    1, 125, 124, // elem 122
    1, 126, 125, // elem 123
    1, 127, 126, // elem 124
    1, 128, 127, // elem 125
    1, 129, 128, // elem 126
    1, 130, 129, // elem 127
    1, 131, 130, // elem 128
    1, 132, 131, // elem 129
    1, 133, 132, // elem 130
    1, 134, 133, // elem 131
    1, 135, 134, // elem 132
    1, 136, 135, // elem 133
    1, 137, 136, // elem 134
    1, 138, 137, // elem 135
    1, 139, 138, // elem 136
    1, 140, 139, // elem 137
    1, 141, 140, // elem 138
    1, 142, 141, // elem 139
    1, 143, 142, // elem 140
    1, 144, 143, // elem 141
    1, 145, 144, // elem 142
    1, 146, 145, // elem 143
    1, 147, 146, // elem 144
    1, 148, 147, // elem 145
    1, 149, 148, // elem 146
    1, 150, 149, // elem 147
    1, 151, 150, // elem 148
    1, 152, 151, // elem 149
    1, 153, 152, // elem 150
    1, 154, 153, // elem 151
    1, 155, 154, // elem 152
    1, 156, 155, // elem 153
    1, 157, 156, // elem 154
    1, 158, 157, // elem 155
    1, 159, 158, // elem 156
    1, 160, 159, // elem 157
    1, 161, 160, // elem 158
    1, 162, 161, // elem 159
    1, 163, 162, // elem 160
    1, 164, 163, // elem 161
    1, 165, 164, // elem 162
    1, 166, 165, // elem 163
    1, 167, 166, // elem 164
    1, 168, 167, // elem 165
    1, 169, 168, // elem 166
    1, 170, 169, // elem 167
    1, 171, 170, // elem 168
    1, 172, 171, // elem 169
    1, 173, 172, // elem 170
    1, 174, 173, // elem 171
    1, 175, 174, // elem 172
    1, 176, 175, // elem 173
    1, 177, 176, // elem 174
    1, 178, 177, // elem 175
    1, 179, 178, // elem 176
    1, 180, 179, // elem 177
    1, 181, 180, // elem 178
    1, 182, 181, // elem 179
    1, 183, 182, // elem 180
    1, 184, 183, // elem 181
    1, 185, 184, // elem 182
    1, 186, 185, // elem 183
    1, 187, 186, // elem 184
    1, 188, 187, // elem 185
    1, 189, 188, // elem 186
    1, 190, 189, // elem 187
    1, 191, 190, // elem 188
    1, 192, 191, // elem 189
    1, 193, 192, // elem 190
    1, 194, 193, // elem 191
    1, 195, 194, // elem 192
    1, 196, 195, // elem 193
    1, 197, 196, // elem 194
    1, 198, 197, // elem 195
    1, 199, 198, // elem 196
    1, 200, 199, // elem 197
    1, 201, 200, // elem 198
    1, 202, 201, // elem 199
    1, 203, 202, // elem 200
    1, 204, 203, // elem 201
    1, 205, 204, // elem 202
    1, 206, 205, // elem 203
    1, 207, 206, // elem 204
    1, 208, 207, // elem 205
    1, 209, 208, // elem 206
    1, 210, 209, // elem 207
    1, 211, 210, // elem 208
    1, 212, 211, // elem 209
    1, 213, 212, // elem 210
    1, 214, 213, // elem 211
    1, 215, 214, // elem 212
    1, 216, 215, // elem 213
    1, 217, 216, // elem 214
    1, 218, 217, // elem 215
    1, 219, 218, // elem 216
    1, 220, 219, // elem 217
    1, 221, 220, // elem 218
    1, 222, 221, // elem 219
    1, 223, 222, // elem 220
    1, 224, 223, // elem 221
    1, 225, 224, // elem 222
    1, 226, 225, // elem 223
    1, 227, 226, // elem 224
    1, 228, 227, // elem 225
    1, 229, 228, // elem 226
    1, 230, 229, // elem 227
    1, 231, 230, // elem 228
    1, 232, 231, // elem 229
    1, 233, 232, // elem 230
    1, 234, 233, // elem 231
    1, 235, 234, // elem 232
    1, 236, 235, // elem 233
    1, 237, 236, // elem 234
    1, 238, 237, // elem 235
    1, 239, 238, // elem 236
    1, 240, 239, // elem 237
    1, 241, 240, // elem 238
    1, 122, 241, // elem 239
    2, 122, 123, // elem 240
    2, 123, 3, // elem 241
    3, 123, 124, // elem 242
    3, 124, 4, // elem 243
    4, 124, 125, // elem 244
    4, 125, 5, // elem 245
    5, 125, 126, // elem 246
    5, 126, 6, // elem 247
    6, 126, 127, // elem 248
    6, 127, 7, // elem 249
    7, 127, 128, // elem 250
    7, 128, 8, // elem 251
    8, 128, 129, // elem 252
    8, 129, 9, // elem 253
    9, 129, 130, // elem 254
    9, 130, 10, // elem 255
    10, 130, 131, // elem 256
    10, 131, 11, // elem 257
    11, 131, 132, // elem 258
    11, 132, 12, // elem 259
    12, 132, 133, // elem 260
    12, 133, 13, // elem 261
    13, 133, 134, // elem 262
    13, 134, 14, // elem 263
    14, 134, 135, // elem 264
    14, 135, 15, // elem 265
    15, 135, 136, // elem 266
    15, 136, 16, // elem 267
    16, 136, 137, // elem 268
    16, 137, 17, // elem 269
    17, 137, 138, // elem 270
    17, 138, 18, // elem 271
    18, 138, 139, // elem 272
    18, 139, 19, // elem 273
    19, 139, 140, // elem 274
    19, 140, 20, // elem 275
    20, 140, 141, // elem 276
    20, 141, 21, // elem 277
    21, 141, 142, // elem 278
    21, 142, 22, // elem 279
    22, 142, 143, // elem 280
    22, 143, 23, // elem 281
    23, 143, 144, // elem 282
    23, 144, 24, // elem 283
    24, 144, 145, // elem 284
    24, 145, 25, // elem 285
    25, 145, 146, // elem 286
    25, 146, 26, // elem 287
    26, 146, 147, // elem 288
    26, 147, 27, // elem 289
    27, 147, 148, // elem 290
    27, 148, 28, // elem 291
    28, 148, 149, // elem 292
    28, 149, 29, // elem 293
    29, 149, 150, // elem 294
    29, 150, 30, // elem 295
    30, 150, 151, // elem 296
    30, 151, 31, // elem 297
    31, 151, 152, // elem 298
    31, 152, 32, // elem 299
    32, 152, 153, // elem 300
    32, 153, 33, // elem 301
    33, 153, 154, // elem 302
    33, 154, 34, // elem 303
    34, 154, 155, // elem 304
    34, 155, 35, // elem 305
    35, 155, 156, // elem 306
    35, 156, 36, // elem 307
    36, 156, 157, // elem 308
    36, 157, 37, // elem 309
    37, 157, 158, // elem 310
    37, 158, 38, // elem 311
    38, 158, 159, // elem 312
    38, 159, 39, // elem 313
    39, 159, 160, // elem 314
    39, 160, 40, // elem 315
    40, 160, 161, // elem 316
    40, 161, 41, // elem 317
    41, 161, 162, // elem 318
    41, 162, 42, // elem 319
    42, 162, 163, // elem 320
    42, 163, 43, // elem 321
    43, 163, 164, // elem 322
    43, 164, 44, // elem 323
    44, 164, 165, // elem 324
    44, 165, 45, // elem 325
    45, 165, 166, // elem 326
    45, 166, 46, // elem 327
    46, 166, 167, // elem 328
    46, 167, 47, // elem 329
    47, 167, 168, // elem 330
    47, 168, 48, // elem 331
    48, 168, 169, // elem 332
    48, 169, 49, // elem 333
    49, 169, 170, // elem 334
    49, 170, 50, // elem 335
    50, 170, 171, // elem 336
    50, 171, 51, // elem 337
    51, 171, 172, // elem 338
    51, 172, 52, // elem 339
    52, 172, 173, // elem 340
    52, 173, 53, // elem 341
    53, 173, 174, // elem 342
    53, 174, 54, // elem 343
    54, 174, 175, // elem 344
    54, 175, 55, // elem 345
    55, 175, 176, // elem 346
    55, 176, 56, // elem 347
    56, 176, 177, // elem 348
    56, 177, 57, // elem 349
    57, 177, 178, // elem 350
    57, 178, 58, // elem 351
    58, 178, 179, // elem 352
    58, 179, 59, // elem 353
    59, 179, 180, // elem 354
    59, 180, 60, // elem 355
    60, 180, 181, // elem 356
    60, 181, 61, // elem 357
    61, 181, 182, // elem 358
    61, 182, 62, // elem 359
    62, 182, 183, // elem 360
    62, 183, 63, // elem 361
    63, 183, 184, // elem 362
    63, 184, 64, // elem 363
    64, 184, 185, // elem 364
    64, 185, 65, // elem 365
    65, 185, 186, // elem 366
    65, 186, 66, // elem 367
    66, 186, 187, // elem 368
    66, 187, 67, // elem 369
    67, 187, 188, // elem 370
    67, 188, 68, // elem 371
    68, 188, 189, // elem 372
    68, 189, 69, // elem 373
    69, 189, 190, // elem 374
    69, 190, 70, // elem 375
    70, 190, 191, // elem 376
    70, 191, 71, // elem 377
    71, 191, 192, // elem 378
    71, 192, 72, // elem 379
    72, 192, 193, // elem 380
    72, 193, 73, // elem 381
    73, 193, 194, // elem 382
    73, 194, 74, // elem 383
    74, 194, 195, // elem 384
    74, 195, 75, // elem 385
    75, 195, 196, // elem 386
    75, 196, 76, // elem 387
    76, 196, 197, // elem 388
    76, 197, 77, // elem 389
    77, 197, 198, // elem 390
    77, 198, 78, // elem 391
    78, 198, 199, // elem 392
    78, 199, 79, // elem 393
    79, 199, 200, // elem 394
    79, 200, 80, // elem 395
    80, 200, 201, // elem 396
    80, 201, 81, // elem 397
    81, 201, 202, // elem 398
    81, 202, 82, // elem 399
    82, 202, 203, // elem 400
    82, 203, 83, // elem 401
    83, 203, 204, // elem 402
    83, 204, 84, // elem 403
    84, 204, 205, // elem 404
    84, 205, 85, // elem 405
    85, 205, 206, // elem 406
    85, 206, 86, // elem 407
    86, 206, 207, // elem 408
    86, 207, 87, // elem 409
    87, 207, 208, // elem 410
    87, 208, 88, // elem 411
    88, 208, 209, // elem 412
    88, 209, 89, // elem 413
    89, 209, 210, // elem 414
    89, 210, 90, // elem 415
    90, 210, 211, // elem 416
    90, 211, 91, // elem 417
    91, 211, 212, // elem 418
    91, 212, 92, // elem 419
    92, 212, 213, // elem 420
    92, 213, 93, // elem 421
    93, 213, 214, // elem 422
    93, 214, 94, // elem 423
    94, 214, 215, // elem 424
    94, 215, 95, // elem 425
    95, 215, 216, // elem 426
    95, 216, 96, // elem 427
    96, 216, 217, // elem 428
    96, 217, 97, // elem 429
    97, 217, 218, // elem 430
    97, 218, 98, // elem 431
    98, 218, 219, // elem 432
    98, 219, 99, // elem 433
    99, 219, 220, // elem 434
    99, 220, 100, // elem 435
    100, 220, 221, // elem 436
    100, 221, 101, // elem 437
    101, 221, 222, // elem 438
    101, 222, 102, // elem 439
    102, 222, 223, // elem 440
    102, 223, 103, // elem 441
    103, 223, 224, // elem 442
    103, 224, 104, // elem 443
    104, 224, 225, // elem 444
    104, 225, 105, // elem 445
    105, 225, 226, // elem 446
    105, 226, 106, // elem 447
    106, 226, 227, // elem 448
    106, 227, 107, // elem 449
    107, 227, 228, // elem 450
    107, 228, 108, // elem 451
    108, 228, 229, // elem 452
    108, 229, 109, // elem 453
    109, 229, 230, // elem 454
    109, 230, 110, // elem 455
    110, 230, 231, // elem 456
    110, 231, 111, // elem 457
    111, 231, 232, // elem 458
    111, 232, 112, // elem 459
    112, 232, 233, // elem 460
    112, 233, 113, // elem 461
    113, 233, 234, // elem 462
    113, 234, 114, // elem 463
    114, 234, 235, // elem 464
    114, 235, 115, // elem 465
    115, 235, 236, // elem 466
    115, 236, 116, // elem 467
    116, 236, 237, // elem 468
    116, 237, 117, // elem 469
    117, 237, 238, // elem 470
    117, 238, 118, // elem 471
    118, 238, 239, // elem 472
    118, 239, 119, // elem 473
    119, 239, 240, // elem 474
    119, 240, 120, // elem 475
    120, 240, 241, // elem 476
    120, 241, 121, // elem 477
    121, 241, 122, // elem 478
    121, 122, 2 // elem 479
};

int QtRigidCylinderObject::init(
	double _x,
	double _y,
	double _z,
	double h,
	double r,
	QVector3D& c
	)
{
    color = c;

    gl.glGenVertexArrays(1, &vao);
    gl.glBindVertexArray(vao);

    const float* pcn = cylinder_nodes;
    NodeCoord &nd0 = node_coords[0];
    nd0.x = 0.0;
    nd0.y = 0.0;
    nd0.z = pcn[2] * h;
    pcn += 3;
    NodeCoord&nd1 = node_coords[1];
    nd1.x = 0.0;
    nd1.y = 0.0;
    nd1.z = pcn[2] * h;
    pcn += 3;
    for (size_t n_id = 2; n_id < cylinder_node_num; ++n_id)
    {
        NodeCoord &nc = node_coords[n_id];
        nc.x = pcn[0] * r;
        nc.y = pcn[1] * r;
        nc.z = pcn[2] * h;
        pcn += 3;
    }

    vbo_index_num = 3 * cylinder_elem_num;

    const unsigned int* pce = cylinder_elems;
    NodeData *pnd = node_datas;
    GLfloat v12_x, v12_y, v12_z, v13_x, v13_y, v13_z;
    GLfloat v_norm, _nx, _ny, _nz;
    for (size_t e_id = 0; e_id < cylinder_elem_num; ++e_id)
    {
        NodeCoord &n1 = node_coords[pce[0]];
        NodeCoord &n2 = node_coords[pce[1]];
        NodeCoord &n3 = node_coords[pce[2]];
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

int QtRigidCylinderObject::update(
	double x,
	double y,
	double z
	)
{
    form_model_mat(Point3D(x, y, z));
    return 0;
}

void QtRigidCylinderObject::draw(QOpenGLShaderProgram& shader)
{
    shader.bind();
    shader.setUniformValue("model_mat", model_mat);
    shader.setUniformValue("g_color", color);

    gl.glBindVertexArray(vao);

    //gl.glFrontFace(GL_CW);
    //gl.glCullFace(GL_BACK);
    gl.glDrawArrays(GL_TRIANGLES, 0, vbo_index_num);
}

inline void QtRigidCylinderObject::form_model_mat(const Point3D& cen)
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
