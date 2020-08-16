#ifndef __Fixed_Size_Mem_Allocator_h__
#define __Fixed_Size_Mem_Allocator_h__

#include <mutex>

template <typename Item>
class FixedSizeMemAllocator
{
protected:
	struct MemPageHeader { MemPageHeader *next; };
	struct MemHeader { MemHeader* next; };

	size_t mem_size;
	MemHeader* empty_mems;
	MemHeader* used_mems, * last_used_mems;
	
	size_t mem_num_per_page;
	size_t mem_page_size;
	MemPageHeader *mem_pages;

	std::mutex alloc_lock; // for thread-saft allocation
	
	inline void add_empty_mem(MemHeader *mem)
	{
		mem->next = empty_mems;
		empty_mems = mem;
	}

	inline MemHeader* del_empty_mem()
	{
		MemHeader* res = empty_mems;
		empty_mems = empty_mems->next;
		return res;
	}

	inline void add_used_mem(MemHeader *mem)
	{
		if (!used_mems)
			last_used_mems = mem;
		mem->next = used_mems;
		used_mems = mem;
	}

	// only called when there is no empty mems
	inline void alloc_mem_page()
	{
		char* new_mem = new char[mem_page_size];
		
		MemPageHeader* new_page = reinterpret_cast<MemPageHeader *>(new_mem);
		new_page->next = mem_pages;
		mem_pages = new_page;

		char* cur_mem = new_mem + sizeof(MemPageHeader);
		MemHeader *mem_hd;
		for (size_t m_id = 0; m_id < mem_num_per_page; ++m_id)
		{
			mem_hd = reinterpret_cast<MemHeader*>(cur_mem);
			add_empty_mem(mem_hd);
			cur_mem += mem_size;
		}
	}

	inline void del_all_mem_pages()
	{
		MemPageHeader *tmp_page;
		MemPageHeader *cur_page = mem_pages;
		while (cur_page)
		{
			tmp_page = cur_page;
			cur_page = cur_page->next;
			delete[] (char *)tmp_page;
		}
	}

public:
	FixedSizeMemAllocator():
		mem_page_size(0), mem_pages(nullptr),
		mem_size(0), mem_num_per_page(1),
		empty_mems(nullptr), used_mems(nullptr) {}

	~FixedSizeMemAllocator() { del_all_mem_pages(); }

	inline void init_mem(size_t _mem_item_num, size_t _mem_num_per_page)
	{
		mem_size = sizeof(MemHeader) + _mem_item_num * sizeof(Item);
		mem_num_per_page = _mem_num_per_page;
		mem_page_size = sizeof(MemPageHeader) + mem_size * mem_num_per_page;
		alloc_mem_page();
	}

	// needs to be thread safe
	inline Item *alloc_mem()
	{
		std::lock_guard<std::mutex> alloc_guard(alloc_lock);
		if (!empty_mems)
			alloc_mem_page();
		MemHeader *res = del_empty_mem();
		add_used_mem(res);
		return reinterpret_cast<Item*>((char *)res + sizeof(MemHeader));
	}

	inline void reset() noexcept
	{
		if (used_mems)
		{
			last_used_mems->next = empty_mems;
			empty_mems = used_mems;
			used_mems = nullptr;
		}
	}

	inline Item* get_first_mem() noexcept { return last_used_mems; }
};

#endif