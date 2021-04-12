#ifndef __Mat_Model_Utils_h__
#define __Mat_Model_Utils_h__

#ifndef __DOUBLE_PRECISION_MAT_MODEL_
#define ffmat(num) (num##f)
#else
#define ffmat(num) (num)
#endif

#define PI (ffmat(3.14159265358979311599796346854))
#define sqrt2 (ffmat(1.41421356237309514547462185874))
#define sqrt3 (ffmat(1.73205080756887719317660412344))
#define sqrt6 (ffmat(2.44948974278317788133563226438))

#endif