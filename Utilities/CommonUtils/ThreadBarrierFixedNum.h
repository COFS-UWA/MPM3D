#ifndef __Thread_Barrier_Fixed_Num_h__
#define __Thread_Barrier_Fixed_Num_h__

#include <atomic>
#include <thread>

class ThreadBarrierFixedNum
{
protected:
	// total thread number
	unsigned int thread_num;
	// how many thread left to wait for
	std::atomic<unsigned int> thread_left_num;
	std::atomic<size_t> generation;

public:
	explicit ThreadBarrierFixedNum(unsigned int num = 1) :
		thread_num(num), thread_left_num(num), generation(0) {}
	~ThreadBarrierFixedNum() {}
	ThreadBarrierFixedNum(const ThreadBarrierFixedNum &other) = delete;
	ThreadBarrierFixedNum &operator=(const ThreadBarrierFixedNum &other) = delete;

	inline unsigned int get_thread_num() { return thread_num; }
	inline void set_thread_num(unsigned int num)
	{
		thread_num = num;
		thread_left_num = num;
	}

	inline void wait_at_barrier()
	{
		const unsigned long long cur_generation = generation.load();
		if (!--thread_left_num)
		{	// all threads have arrive
			thread_left_num = thread_num;
			++generation;
		}
		else
		{	// wait for other threads
			while (generation.load() == cur_generation);
		}
	}

	inline void wait()
	{
		const unsigned long long cur_generation = generation.load();
		--thread_left_num;
		while (generation.load() == cur_generation);
			//std::this_thread::yield();
	}

	inline void wait_for_others()
	{
		while (thread_left_num.load() > 1);
			//std::this_thread::yield();
	}

	inline void lift_barrier()
	{
		thread_left_num.store(thread_num);
		++generation;
	}
};

#endif