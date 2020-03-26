#ifndef __Memory_Utils_Memory_Buffer_h__
#define __Memory_Utils_Memory_Buffer_h__

#include <string>

namespace MemoryUtils
{
	/*===============================================
	 * Class MemoryBuffer
	 *-----------------------------------------------
	 * A simple buffer to hold a chunk of raw memory 
	 * data of randow size.
	 *===============================================*/
	class MemoryBuffer
	{
	protected:
		char *buffer;
		size_t size;
		size_t capacity;
	public:
		MemoryBuffer(size_t init_size = 0) :
			buffer(nullptr), size(0), capacity(0)
		{
			if (init_size)
			{
				buffer = new char[init_size];
				capacity = init_size;
			}
		}
		~MemoryBuffer() { clear(); }
		void clear(void)
		{
			if (buffer)
			{
				delete[] buffer;
				buffer = nullptr;
			}
			size = 0;
			capacity = 0;
		}
		void reserve(size_t reserve_size)
		{
			if (reserve_size > capacity)
			{
				char *tmp = new char[reserve_size];
				if (buffer)
				{
					if (size)
						memcpy(tmp, buffer, size);
					delete[] buffer;
				}
				buffer = tmp;
				capacity = reserve_size;
			}
		}
		inline void expand(size_t expand_size)
		{
			reserve(capacity + expand_size);
		}
		inline void *get_buffer(void) const { return buffer; }
		inline size_t get_size(void) const { return size; }
		inline size_t get_capacity(void) const { return capacity; }
		inline void reset(void)	{ size = 0; }
		void add_data(const void *data, size_t data_size)
		{
			size_t new_size = size + data_size;
			if (new_size > capacity)
				reserve(new_size);
			memcpy(buffer + size, data, data_size);
			size = new_size;
		}
	};
};

#endif