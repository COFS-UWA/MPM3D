#include "SimulationsOMP_pcp.h"

#include <cstring>

#include "DataMem.h"

namespace Util
{
	DataMem::DataMem() : _mem(nullptr), _size(0) {}

	DataMem::DataMem(size_t num, DataSize dsize) : _mem(nullptr), _size(0) { alloc(num, dsize); }

	DataMem::DataMem(size_t num) : _mem(nullptr), _size(0) { alloc(num); }

	DataMem::DataMem(DataMem& other) { move_from(other); }

	DataMem::~DataMem() { free(); }

	void* DataMem::alloc(size_t num)
	{
		free();
		if (num)
		{
			_size = num;
			_mem = new char[_size];
			return _mem;
		}
		return nullptr;
	}

	void* DataMem::alloc(size_t num, DataSize dsize)
	{
		free();
		if (num)
		{
			_size = num * size_t(dsize);
			_mem = new char[_size];
			return _mem;
		}
		return nullptr;
	}

	void DataMem::free()
	{
		if (_mem)
		{
			delete[] _mem;
			_mem = nullptr;
			_size = 0;
		}
	}

	bool DataMem::copy_from(const DataMem& other)
	{
		if (_size < other._size)
			return false;
		memcpy(_mem, other._mem, other._size);
		return true;
	}

	void DataMem::move_from(DataMem& other)
	{
		_mem = other._mem;
		_size = other._size;
		other._mem = nullptr;
		other._size = 0;
	}

}