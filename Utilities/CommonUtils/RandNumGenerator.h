#ifndef __Rand_Num_Generator__
#define __Rand_Num_Generator__

#include <cstdlib>

class RandNumGenerator
{
private:
	// singleton
	RandNumGenerator() { srand(1); }
	~RandNumGenerator() {}
	// no copy
	RandNumGenerator(RandNumGenerator& other) {}
	RandNumGenerator(const RandNumGenerator& other) {}
	RandNumGenerator& operator=(RandNumGenerator& other);
	RandNumGenerator& operator=(const RandNumGenerator& other);
	
	static RandNumGenerator inst;

	// generate random int between [0, RAND_MAX]
	inline int gen_int_rand() { return rand(); }
	// generate random int between [0, max]
	inline int gen_int_rand(int max) { return rand() % (max + 1); }
	// generate random int between [min, max]
	inline int gen_int_rand(int min, int max) { return min + rand() % (max-min+1); }
	// generate random double between [0.0, 1.0]
	inline double gen_double_rand()
	{
		return 1.0 / double(RAND_MAX) * double(rand());
	}
	// generate random double between [0.0, max]
	inline double gen_double_rand(double max)
	{
		return max / double(RAND_MAX) * double(rand());
	}
	// generate random double between [min, max]
	inline double gen_double_rand(double min, double max)
	{
		return min + (max - min) / double(RAND_MAX) * double(rand());
	}

public:
	inline static void set_seed(unsigned int seed) { srand(seed); }
	inline static int get_int() { return inst.gen_int_rand(); }
	inline static int get_int(int max) { return inst.gen_int_rand(max); }
	inline static int get_int(int min, int max) { return inst.gen_int_rand(min, max); }
	inline static double get_double() { return inst.gen_double_rand(); }
	inline static double get_double(double max) { return inst.gen_double_rand(max); }
	inline static double get_double(double min, double max) { return inst.gen_double_rand(min, max); }
};

#endif