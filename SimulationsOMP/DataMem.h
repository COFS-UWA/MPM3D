#pragma once

namespace Util
{
	class DataMem
	{
	public:
		enum DataSize
		{
			Byte1 = 1,
			Byte2 = 2,
			Byte4 = 4,
			Byte8 = 8,
		};

		DataMem();
		DataMem(size_t num, DataSize dsize);
		DataMem(size_t num);
		DataMem(DataMem& other);
		DataMem(const DataMem& other) = delete;
		~DataMem();

		void* alloc(size_t num);
		void* alloc(size_t num, DataSize dsize);
		template <typename Type>
		Type* alloc(size_t num)
		{
			free();
			if (num)
			{
				_size = num * sizeof(Type);
				_mem = new char[_size];
				return (Type*)_mem;
			}
			return nullptr;
		}

		void free();
		bool copy_from(const DataMem& other);
		void move_from(DataMem& other);

		inline void* mem() noexcept { return _mem; }
		template <typename Type>
		inline Type* mem() noexcept { return (Type*)_mem; }
		inline size_t size() const noexcept { return _size; }

	protected:
		void* _mem;
		size_t _size;
	};
}