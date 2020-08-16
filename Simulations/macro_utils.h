#ifndef __macro_utils_h__
#define __macro_utils_h__

#define INIT_BC_TEMPLATE(name, type)    \
	inline void init_ ## name ## s(size_t num) \
	{                                   \
		if (name ## s)                  \
		{                               \
			if (name ## _num < num)		\
				delete[] name ## s;		\
			else                        \
			{                           \
				name ## _num = num;	    \
				return;                 \
			}                           \
		}                               \
		name ## s = new type ## [num];  \
		name ## _num = num;             \
	}                                   \
	inline void clear_ ## name ## s()   \
	{                                   \
		if (name ## s)                  \
		{                               \
			delete[] name ## s;         \
			name ## s = nullptr;        \
		}                               \
		name ## _num = 0;               \
	}                                   \
	inline type *get_ ## name ## s()    \
	{ return name ## s; }               \
	inline size_t get_ ## name ## _num() \
	{ return name ## _num; }
	
#define N_tol (1.0e-10)

#define ATOM_INC(var) ((var).fetch_add(1, std::memory_order_relaxed))

#endif