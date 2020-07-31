#ifndef __Thread_Barrier_Fixed_Num_h__
#define __Thread_Barrier_Fixed_Num_h__

#include <atomic>
#include <thread>
#include <iostream>

class ThreadBarrierFixedNum
{
protected:
	unsigned int thread_num;
	std::atomic<size_t> generation;
	// number of threads left to wait for
	std::atomic<unsigned int> thread_left_num;

public:
	explicit ThreadBarrierFixedNum(unsigned int num = 1) :
		thread_num(num), thread_left_num(num), generation(0) {}
	~ThreadBarrierFixedNum() {}
	ThreadBarrierFixedNum(const ThreadBarrierFixedNum &other) = delete;
	ThreadBarrierFixedNum &operator=(const ThreadBarrierFixedNum &other) = delete;

	inline unsigned int get_thread_num() { return thread_num; }
	inline size_t get_generation() { return generation.load(); } // for debug
	inline unsigned int get_thread_left_num() { return thread_left_num.load(); } // for debug

	inline void set_thread_num(unsigned int num)
	{
		thread_num = num;
		thread_left_num = num;
	}

	inline void wait()
	{
		const size_t cur_generation = generation;
		thread_left_num.fetch_sub(1, std::memory_order_relaxed);
		while (generation.load(std::memory_order_relaxed) == cur_generation);
	}

	inline void wait_for_others()
	{
		while (thread_left_num.load(std::memory_order_relaxed) > 1);
	}

	inline void lift_barrier()
	{
		thread_left_num.store(thread_num, std::memory_order_relaxed);
		generation.fetch_add(1, std::memory_order_relaxed);
	}
};

#endif