#ifndef __Para_Util_Serial_Sort_h__
#define __Para_Util_Serial_Sort_h__

namespace ParaUtil
{
namespace Internal
{
    constexpr size_t max_data_num_for_inseration_sort = 32;
    
    void insertion_sort(size_t keys[], size_t values[], size_t data_num);
    
    void in_place_lsd_sort(size_t keys[], size_t values[], size_t data_num,
        size_t key_bit_offset, size_t radix_bit_num, size_t sum_bin[], size_t end_bin[]);
    
    void serial_sort(size_t keys[], size_t values[], size_t data_num,
        size_t key_bit_num, size_t tmp_keys[], size_t tmp_values[]);
}
}

#endif