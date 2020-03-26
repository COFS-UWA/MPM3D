#ifndef __Memory_Utils_Item_Stack_hpp__
#define __Memory_Utils_Item_Stack_hpp__

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
	/*================================================
	Class ItemStack
	--------------------------------------------------
	1. Use linked memory pool;
	2. Set pre_alloc_size != 0 to get in-stack memory;
	3. set reasonable large page_size for efficiency.
	=================================================*/
	template <typename Item, size_t fold = 1, size_t pre_alloc_size = 0>
	class ItemStack
	{
	protected:
		struct MemPageHeader
		{
			MemPageHeader *prev, *next;
			Item *start, *end;
		};
		// first in-stack memory page
		union
		{
			char first_page_mem[sizeof(MemPageHeader)
								+ MEMORY_ALIGNMENT
								+ sizeof(Item)*(pre_alloc_size+1)];
			MemPageHeader first_page;
			struct
			{
				MemPageHeader *first_page_prev;
				MemPageHeader *first_page_next;
				Item *first_page_start;
				Item *first_page_end;
			};
		};
		size_t base_page_size, page_size;
		bool need_optimized;
		MemPageHeader *cur_page;
		Item *cur, *start, *end;

	public:
		ItemStack(size_t init_page_size = pre_alloc_size) :
			base_page_size(init_page_size ? init_page_size : 1), page_size(base_page_size),
			need_optimized(false),
			first_page_prev(&first_page), first_page_next(&first_page),
			first_page_start((Item *)MEMORY_ALIGNED_ADDRESS(size_t(first_page_mem)+sizeof(MemPageHeader))),
			first_page_end(first_page_start + pre_alloc_size + 1),
			cur_page(&first_page), cur(first_page_start),
			start(first_page_start), end(first_page_end) {}
		~ItemStack() { clear(); }
		inline void set_page_size(size_t init_page_size) noexcept
		{
			base_page_size = init_page_size ? init_page_size : 1;
			page_size = base_page_size;
		}
		inline size_t get_page_size(void) const noexcept { return page_size; }
		inline void push(Item *pitem)
		{
			if (cur < end)
			{
				*cur = *pitem;
				++cur;
				return;
			}
			if (cur_page->next == &first_page)
				alloc_new_memory_page();
			cur_page = cur_page->next;
			start = cur_page->start;
			cur = start;
			end = cur_page->end;
			*cur = *pitem;
			++cur;
			return;
		}
		inline void push(Item &item)
		{
			if (cur < end)
			{
				*cur = item;
				++cur;
				return;
			}
			if (cur_page->next == &first_page)
				alloc_new_memory_page();
			cur_page = cur_page->next;
			start = cur_page->start;
			cur = start;
			end = cur_page->end;
			*cur = item;
			++cur;
			return;
		}
		inline Item *pop(void)
		{
			--cur;
			if (cur < start)
			{
				if (cur_page == &first_page)
				{
					cur = start;
					return nullptr;
				}
				cur_page = cur_page->prev;
				start = cur_page->start;
				end = cur_page->end;
				cur = end - 1;
			}
			return cur;
		}
		inline bool is_empty(void) const
		{
			return cur_page == &first_page && cur = start;
		}
		inline void reset(void)
		{
			cur_page = &first_page;
			start = first_page_start;
			cur = start;
			end = first_page_end;
		}
		// contract allocated memory pages into one
		void reset_optimize(size_t additional_stack_size = 0)
		{
			if (need_optimized)
			{
				// cal total allocated size
				size_t total_size = 0;
				for (MemPageHeader *pg_iter = first_page.next;
					 pg_iter != &first_page; pg_iter = pg_iter->next)
					total_size += pg_iter->end - pg_iter->start;
				total_size += additional_stack_size;
				// clear list
				clear();
				// alloc new page
				union { char *mem; MemPageHeader *mem_page; };
				// take memory of the first page as "preallocated space"
				total_size += pre_alloc_size + 1;
				size_t char_num = sizeof(MemPageHeader) + MEMORY_ALIGNMENT + sizeof(Item)*(total_size + 1);
				mem = new char[char_num];
				first_page_start = (Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem) + sizeof(MemPageHeader));
				first_page_end = first_page_start + total_size;
				mem_page->start = first_page_end;
				mem_page->end = mem_page->start + 1;
				// add page to list
				mem_page->next = &first_page;
				mem_page->prev = &first_page;
				first_page.next = mem_page;
				first_page.prev = mem_page;
				need_optimized = false;
			}
			// reset state
			cur_page = &first_page;
			start = first_page_start;
			end = first_page_end;
			cur = start;
		}
		void clear(void)
		{
			// reset page size
			page_size = base_page_size;
			// clear list
			MemPageHeader *&top_page = first_page.next;
			MemPageHeader *tmp_page = top_page;
			while (tmp_page != &first_page)
			{
				top_page = top_page->next;
				delete[] (char*)tmp_page;
				tmp_page = top_page;
			}
			// reused in-stack memory
			first_page_start = ((Item *)MEMORY_ALIGNED_ADDRESS(size_t(first_page_mem) + sizeof(MemPageHeader)));
			first_page_end = first_page_start + pre_alloc_size + 1;
			first_page.prev = &first_page;
			first_page.next = &first_page;
			need_optimized = false;
			// reset state
			cur_page = &first_page;
			start = first_page_start;
			cur = start;
			end = first_page_end;
		}
		// contract into one memory page
		inline void optimize(size_t additional_stack_size = 0)
		{
			MemPageHeader *pg_iter;
			if (need_optimized)
			{
				// cal total allocated size
				size_t total_size = 0;
				for (pg_iter = first_page.next; pg_iter != &first_page; pg_iter = pg_iter->next)
					total_size += pg_iter->end - pg_iter->start;
				total_size += additional_stack_size;
				// alloc new page (if necessary)
				union { char *mem; MemPageHeader *mem_page; };
				total_size += pre_alloc_size;
				size_t char_num = sizeof(MemPageHeader) + MEMORY_ALIGNMENT + sizeof(Item)*(total_size+1);
				mem = new char[char_num];
				// copy content from old pages to new page
				mem_page->start = (Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem) + sizeof(MemPageHeader));
				char *cur_page_content = (char *)mem_page->start;
				size_t cur_offset = 0, cur_page_size;
				for (pg_iter = &first_page; pg_iter != cur_page; pg_iter = pg_iter->next)
				{
					cur_offset += pg_iter->end - pg_iter->start;
					cur_page_size = (pg_iter->end - pg_iter->start) * sizeof(Item);
					if (cur_page_size)
						memcpy(cur_page_content, pg_iter->start, cur_page_size);
					cur_page_content += cur_page_size;
				}
				cur_offset += cur - start;
				cur_page_size = (cur - start) * sizeof(Item);
				if (cur_page_size)
					memcpy(cur_page_content, start, cur_page_size);
				// clear original page
				clear();
				// set page memory range
				first_page.start = mem_page->start;
				first_page.end = first_page.start + total_size;
				mem_page->start = first_page.end;
				mem_page->end = mem_page->start + 1;
				// add page to list
				mem_page->next = &first_page;
				mem_page->prev = &first_page;
				first_page.next = mem_page;
				first_page.prev = mem_page;
				// modify end, cur, start
				cur_page = &first_page;
				start = first_page.start;
				cur = start + cur_offset;
				end = first_page.end;
				// reset state
				need_optimized = false;
			}
		}

	protected:
		void alloc_new_memory_page(void)
		{
			// alloc new page
			union { char *mem; MemPageHeader *mem_page; };
			size_t char_num = sizeof(MemPageHeader) + MEMORY_ALIGNMENT + sizeof(Item)*page_size;
			mem = new char[char_num];
			// init new page
			mem_page->start = (Item *)MEMORY_ALIGNED_ADDRESS(size_t(mem) + sizeof(MemPageHeader));
			mem_page->end = mem_page->start + page_size;
			mem_page->next = &first_page;
			MemPageHeader *&ori_last_page = first_page.prev;
			mem_page->prev = ori_last_page;
			ori_last_page->next = mem_page;
			ori_last_page = mem_page;
			need_optimized = true;
			// update page_size
			page_size *= fold;
		}
	};
};

#undef MEMORY_ALIGNMENT
#undef MEMORY_ALIGNMENT_PADDING
#undef MEMORY_ALIGNED_ADDRESS

#endif