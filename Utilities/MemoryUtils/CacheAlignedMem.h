#ifndef __Cache_Aligned_Memory_h__
#define __Cache_Aligned_Memory_h__

#define Cache_Alignment 0x40 // 64

#define Cache_Mask (0x40 - 1)

#define Cache_Alignment_Padding(address) \
	((Cache_Alignment - ((address) & Cache_Mask)) & Cache_Mask)

#define Cache_Offset_Padding(address, offset) \
	((Cache_Alignment + (offset) - ((address) & Cache_Mask)) & Cache_Mask)

#define Cache_Aligned_Address(address) \
	((void *)(size_t(address) + Cache_Alignment_Padding(size_t(address))))

#define Cache_Offset_Address(address, offset) \
	((void *)(size_t(address) + Cache_Offset_Padding(size_t(address), offset)))

template <typename Item>
inline Item *cache_aligned(Item *address)
{
	return (Item *)Cache_Aligned_Address(address);
}

template <typename Item>
inline Item* cache_aligned_offset(Item* address, size_t offset)
{
	return Cache_Offset_Address(address, offset);
}

class CacheAlignedMem
{
protected:
	char *mem_head;

public:
	inline explicit CacheAlignedMem() : mem_head(nullptr) {}
	inline ~CacheAlignedMem() noexcept
	{
		if (mem_head)
		{
			delete[] mem_head;
			mem_head = nullptr;
		}
	}

	CacheAlignedMem(const CacheAlignedMem &other) = delete;
	CacheAlignedMem(CacheAlignedMem& other) noexcept
	{
		free();
		mem_head = other.mem_head;
		other.mem_head = nullptr;
	}

	inline void *raw_address() const noexcept { return mem_head; }
	inline void *aligned_address() const noexcept
	{
		return Cache_Aligned_Address(mem_head);
	}
	inline void *aligned_address(size_t offset) const noexcept
	{
		return Cache_Offset_Address(mem_head, offset & Cache_Mask);
	}

	inline void free() noexcept
	{
		if (mem_head)
		{
			delete[] mem_head;
			mem_head = nullptr;
		}
	}

	inline void *alloc(size_t mem_size)
	{
		free();
		mem_head = new char[mem_size + Cache_Alignment];
		return Cache_Aligned_Address(mem_head);
	}

	inline void* alloc(size_t mem_size, size_t offset)
	{
		free();
		mem_head = new char[mem_size + Cache_Alignment];
		return Cache_Offset_Address(mem_head, offset & Cache_Mask);
	}

	template <typename Item>
	inline Item* alloc(size_t num)
	{
		free();
		mem_head = new char[sizeof(Item) * num + Cache_Alignment];
		return (Item *)Cache_Aligned_Address(mem_head);
	}

	template <typename Item>
	inline Item *alloc(size_t num, size_t offset)
	{
		free();
		mem_head = new char[sizeof(Item) * num + Cache_Alignment];
		return Cache_Offset_Address(mem_head, offset & Cache_Mask);
	}
};

#undef Cache_Mask
#undef Cache_Alignment_Padding
#undef Cache_Offset_Padding
#undef Cache_Aligned_Address
#undef Cache_Offset_Address

#endif