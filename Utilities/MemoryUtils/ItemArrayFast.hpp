#ifndef __Meomory_Utils_Item_Array_Fast_hpp__
#define __Meomory_Utils_Item_Array_Fast_hpp__

namespace MemoryUtils
{
	/*===============================================
	* Class ItemArrayFast
	*-----------------------------------------------
	* Only hold a simple chunk of memory.
	* No copy operation when resize.
	*===============================================*/
	template <typename Item, size_t pre_alloc_size = 0>
	class ItemArrayFast
	{
	protected:
		Item *mem_buf;
		size_t mem_size;
		Item pre_alloc_buf[pre_alloc_size];

	public:
		ItemArrayFast() : mem_buf(pre_alloc_buf), mem_size(pre_alloc_size) {}
		~ItemArrayFast()
		{
			if (mem_buf != pre_alloc_buf)
				delete[] mem_buf;
			mem_buf = pre_alloc_buf;
			mem_size = pre_alloc_size;
		}
		inline size_t resize(size_t _size)
		{
			if (_size > mem_size)
			{
				if (mem_buf != pre_alloc_buf) delete[] mem_buf;
				mem_buf = new Item[_size];
				mem_size = _size;
			}
		}
		inline Item *get_mem() noexcept { return mem_buf; }
		inline size_t get_size() const noexcept { return mem_size; }
		inline Item &operator[] (size_t id) { return mem_buf[id]; }
	};

	template <typename Item>
	class ItemArrayFast<Item, 0>
	{
	protected:
		Item *mem_buf;
		size_t mem_size;

	public:
		ItemArrayFast() : mem_buf(nullptr), mem_size(0) {}
		~ItemArrayFast() { clear(); }
		inline void clear()
		{
			if (mem_buf) delete[] mem_buf;
			mem_buf = nullptr;
			mem_size = 0;
		}
		inline Item *resize(size_t _size)
		{
			if (_size > mem_size)
			{
				if (mem_buf) delete[] mem_buf;
				mem_buf = new Item[_size];
				mem_size = _size;
			}
			return mem_buf;
		}
		inline Item *get_mem() noexcept { return mem_buf; }
		inline size_t get_size() const noexcept { return mem_size; }
		inline Item &operator[] (size_t id) { return mem_buf[id]; }
	};
};

#endif