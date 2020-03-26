#ifndef __Memory_Utils_Item_Array_hpp__
#define __Memory_Utils_Item_Array_hpp__

#include <string>

// must be power of two
#define MEMORY_ALIGNMENT sizeof(void *)

// padding to round up the cloest power of two
#define MEMORY_ALIGNMENT_PADDING(address) \
	((MEMORY_ALIGNMENT - ((address) & (MEMORY_ALIGNMENT - 1))) & (MEMORY_ALIGNMENT - 1))

// round up the cloest power of two
#define MEMORY_ALIGNED_ADDRESS(address) \
	((address) + MEMORY_ALIGNMENT_PADDING(address))

namespace MemoryUtils
{
	/*=================================================
	Class ItemArray
	---------------------------------------------------
	Note:
	1. Continuous memory;
	2. Set reasonable large page_size for efficiency;
	3. The memory address will **change** dynamically
	   with reallocation;
	4. fold must > 0.
	==================================================*/
	template<typename Item, size_t fold = 2, size_t pre_alloc_size = 0>
	class ItemArray
	{
	protected:
		char *mem;
		Item *start, *cur, *mem_end;
		size_t base_page_size, page_size;
		char pre_alloc_mem[MEMORY_ALIGNMENT + pre_alloc_size * sizeof(Item)];

	public:
		ItemArray(size_t init_page_size = pre_alloc_size) :
			mem(pre_alloc_mem),
			start((Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem))),
			cur(start), mem_end(start + pre_alloc_size),
			base_page_size(init_page_size ? init_page_size : 1),
			page_size(base_page_size) {}
		~ItemArray() { clear(); }

		inline Item *get_mem(void) const { return start; }
		inline size_t get_num(void) const { return size_t(cur - start); }
		inline size_t get_capacity(void) const { return size_t(mem_end - start); }
		inline Item *get_nth_item(size_t id) { return start + id; }
		inline Item &operator[] (size_t id) { return start[id]; }
		inline void set_page_size(size_t init_page_size)
		{
			base_page_size = init_page_size ? init_page_size : 1;
			page_size = base_page_size;
		}
		inline size_t get_page_size(void) const { return page_size; }

		// Alloc memory for one item
		inline Item *alloc(void)
		{
			if (cur < mem_end)
				return cur++;
			expand_memory();
			return cur++;
		}
		// Add one item
		inline void add(Item &item)
		{
			if (cur < mem_end)
			{
				*(cur++) = item;
				return;
			}
			expand_memory();
			*(cur++) = item;
		}
		inline void add(Item *pitem)
		{
			if (cur < mem_end)
			{
				*(cur++) = *pitem;
				return;
			}
			expand_memory();
			*(cur++) = *pitem;
		}

		// Alloc multiple items at one time
		// Return the address of empty of memory
		inline Item *alloc(size_t num)
		{
			Item *res = cur;
			cur += num;
			if (cur < mem_end)
				return res;
			// not enough memory
			if (page_size < num)
				page_size = num;
			cur = res;
			expand_memory();
			res = cur;
			cur += num;
			return res;
		}

		// Ensure that memory has capacity not less than num
		inline void reserve(size_t num)
		{
			if (num > size_t(mem_end - start))
			{
				if (page_size < num)
					page_size = num;
				expand_memory();
			}
		}
		inline void reset(void)	{ cur = start; }
		void clear(void)
		{
			// reset page size
			page_size = base_page_size;
			// reset memory
			if (mem != pre_alloc_mem)
			{
				delete[] mem;
				mem = pre_alloc_mem;
				start = (Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem));
				mem_end = start + pre_alloc_size;
			}
			cur = start;
		}

	protected:
		void expand_memory(void)
		{
			size_t capacity_new = size_t(mem_end - start) + page_size;
			page_size *= fold;
			char *mem_new = new char[MEMORY_ALIGNMENT + capacity_new * sizeof(Item)];
			Item *start_new = (Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem_new));
			size_t content_len = size_t(cur - start);
			if (content_len)
				memcpy(start_new, start, content_len * sizeof(Item));
			if (mem != pre_alloc_mem)
				delete[] mem;
			mem = mem_new;
			start = start_new;
			cur = start + content_len;
			mem_end = start + capacity_new;
		}
	};
}

#undef MEMORY_ALIGNMENT
#undef MEMORY_ALIGNMENT_PADDING
#undef MEMORY_ALIGNED_ADDRESS

#endif