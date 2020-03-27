#ifndef __Memory_Utils_Pre_Alloc_String_Buffer_hpp__
#define __Memory_Utils_Pre_Alloc_String_Buffer_hpp__

#include <string>

namespace MemoryUtils
{
	/*===============================================
	 * Class PreAllocPreAllocStringBuffer
	 *-----------------------------------------------
	 *   A simple buffer for C style string. This
	 *   buffer has initial static memory chunk with
	 *   size of "init_buffer_size" to avoid memory
	 *   allocation
	 * Note:
	 *   1. The "init_buffer_size" must not be zero
	 *   2. The "length" variables doesn't include 
	 *      '\0' at the end, but "capacity" variables
	 *      take it into consideration.
	 *===============================================*/
	template<size_t init_buffer_size>
	class PreAllocStringBuffer
	{
	protected:
		char *str;
		size_t capacity;
		size_t length;
		char ori_buffer[init_buffer_size];
	public:
		PreAllocStringBuffer() :
			str(ori_buffer), capacity(init_buffer_size), length(0) { str[0] = '\0'; }
		~PreAllocStringBuffer() { if (str != ori_buffer) delete[] str; }

		inline const char *get_string() const { return str; }
		inline const char *c_str() const { return str; }
		inline size_t get_length() const { return length; }
		inline size_t get_capacity() const { return capacity; }

		inline void reserve(size_t reserve_string_length)
		{
			if (reserve_string_length && (++reserve_string_length) > capacity)
				reserve_raw(reserve_string_length);
		}
		inline void expand(size_t expand_string_length)
		{
			if (expand_string_length)
				reserve_raw(capacity + expand_string_length);
		}
		inline void reset()
		{
			length = 0;
			str[0] = '\0';
		}
		void clear()
		{
			if (str != ori_buffer)
			{
				delete[] str;
				str = ori_buffer;
				capacity = init_buffer_size;
			}
			str[0] = '\0';
			length = 0;
		}

		void add_string(const char *another_str)
		{
			size_t new_capacity = length + strlen(another_str) + 1;
			if (new_capacity > capacity)
				reserve_raw(new_capacity);
			strcpy(str + length, another_str);
			length = new_capacity - 1;
		}
		void add_char(char another_char)
		{
			size_t new_capacity = length + 2;
			if (new_capacity > capacity)
				reserve_raw(new_capacity);
			str[length] = another_char;
			++length;
			str[length] = '\0';
		}

	protected:
		void reserve_raw(size_t reserve_size)
		{
			char *tmp = new char[reserve_size];
			if (length)
				strcpy(tmp, str);
			if (str != ori_buffer)
				delete[] str;

			str = tmp;
			capacity = reserve_size;
		}

	public:
		inline PreAllocStringBuffer &operator<< (PreAllocStringBuffer &another_buf)
		{
			add_string(another_buf.get_string());
			return *this;
		}
		inline PreAllocStringBuffer &operator<< (const char *another_str)
		{
			add_string(another_str);
			return *this;
		}
		inline PreAllocStringBuffer &operator<< (char another_char)
		{
			add_char(another_char);
			return *this;
		}
		inline PreAllocStringBuffer &operator<< (std::string &another_str)
		{
			add_string(another_str.c_str());
			return *this;
		}
		inline PreAllocStringBuffer &operator<< (size_t num)
		{
			char str_tmp[25];
			snprintf(str_tmp, 25, "%llu", num);
			add_string(str_tmp);
			return *this;
		}
		inline PreAllocStringBuffer<init_buffer_size> &operator= (const char *another_str)
		{
			reset();
			add_string(another_str);
			return *this;
		}
		inline bool operator== (PreAllocStringBuffer &another_buf)
		{
			return strcmp(str, another_buf.str) == 0;
		}
		inline bool operator!= (PreAllocStringBuffer &another_buf)
		{
			return strcmp(str, another_buf.str) != 0;
		}
		inline bool operator== (const char *another_str)
		{
			return strcmp(str, another_str) == 0;
		}
		inline bool operator!= (const char *another_str)
		{
			return strcmp(str, another_str) != 0;
		}
		inline bool operator== (std::string &another_str)
		{
			return strcmp(str, another_str.c_str()) == 0;
		}
		inline bool operator!= (std::string &another_str)
		{
			return strcmp(str, another_str.c_str()) != 0;
		}
	};

	template<size_t init_buffer_size>
	std::ostream &operator<< (std::ostream& out, 
		PreAllocStringBuffer<init_buffer_size> &str_buf)
	{
		out << str_buf.get_string();
		return out;
	}
};

#endif