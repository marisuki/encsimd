#ifndef SCALAR_H
#define SCALAR_H


#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>


template<typename dtaIn, typename rleIn, typename decOut> 
class scalar {
private:
    dtaIn* dta; rleIn* rle; int minimum; int len; int minimumRLE;
    //int attr;
    int dec_sz;
    decOut ps; //tmp[i] -> attr^i: filtered
    decOut* dec; 
    void decode_ps() {
        int pos = 0, loc = 0, curr_dta = dta[pos] + minimum;
        int rle_rest = rle[0] + minimumRLE;
        while(rle_rest) {
            dec[loc+1] = dec[loc] + curr_dta;
            loc++; rle_rest--;
            if(rle_rest == 0) {
                pos++;
                if(pos < len) {
                    rle_rest = rle[pos] + minimumRLE; 
                    curr_dta = dta[pos] + minimum; 
                }
            }
        }
        dec_sz = loc;
    }
    void apply(int An) {
        for(int i=0;i<len;i++){
            dec[i] += An;
        }
    }
public:
    scalar(dtaIn* delta, rleIn* rle, decOut* decBuff, 
            decOut start, int mDta, int mRLE, int len) {
        this->dta = delta;
        this->rle = rle;
        this->dec = decBuff;
        this->ps = start;
        this->minimumRLE = mRLE;
        this->minimum = mDta;
        this->len = len;
    }
    void decode() {decode_ps();}
    void finalize(int An) {apply(An);}
};

// struct DTARLE { // implementation
//     int* dta; int* rle; int start; int minimum_dta; int minimum_rle;
//     int cps_num; int total_len;
//     int* decoded_data; Schema schema; int attr;
//     DTARLE(){}
//     DTARLE(Schema schema, int attr) {
//         this->attr = attr;
//         cps_num = schema.total_compressed_len_(attr);
//         total_len = schema.total_decoded_len_(attr);
//         decoded_data = (int*) malloc((total_len + 16)*sizeof(int));
//         dta = (int*) malloc((cps_num+ 16)*sizeof(int));
//         rle = (int*) malloc((cps_num+ 16)*sizeof(int));
//         schema.readFiles(attr, cps_num, dta, rle);
//         start = schema.start_(attr);
//         minimum_dta = schema.minimum_delta_(attr);
//         minimum_rle = schema.minimum_rle_(attr);
//     }

//     std::vector<Task> decode_vec(int pc, int b, bool expt=false) {
//         if(pc*b > cps_num) {b = cps_num/pc; }
//         int dec_loc[100];
//         std::vector<Task> exp;
//         for(int i=0;i<pc;i++) {for(int j=0;j<b;j++) dec_loc[i+1] += rle[i*b+j];}
//         #pragma omp parallel for
//         for(int i=0;i<pc;i++) {
//             if(!expt) this->decode_vec16(i*b, dec_loc[i]);
//             else exp.insert(exp.begin(), Task(i, attr, start, i*b, dec_loc[i]));
//         }
//         return exp;
//     }

//     void decode_scalar16() {
//         this->decode_scalar16(0, 0);
//     }
    
//     void decode_scalar16(int first_pair_pos, int first_decoded_pos) {
//         int16_t* dta0 = (int16_t*) dta;
//         int16_t* rle0 = (int16_t*) rle;
//         int16_t* decode = (int16_t*) decoded_data;
//         int loc = first_pair_pos; int rle_rest = rle0[loc] + minimum_rle;
//         int write_bias = first_decoded_pos + 1; 
//         if(first_decoded_pos == 0) decode[first_decoded_pos] = start; 
//         else {
//             decode[first_decoded_pos] = (dta0[loc] + minimum_dta);
//             rle_rest --;
//             if(rle_rest == 0) {
//                 loc ++; rle_rest = rle0[loc] + minimum_rle;
//             }
//         }
//         while(loc < total_len) {
//             decode[write_bias] = decode[write_bias-1] + (dta0[loc] + minimum_dta);
//             rle_rest --;
//             if(rle_rest == 0) {
//                 loc ++; rle_rest = rle0[loc] + minimum_rle;
//             }
//         }
//     }

//     void decode_vec16() {
//         this->decode_vec16(0, 0);
//     }

//     void decode_vec16(int first_pair_pos, int first_decoded_pos) {
//         int16_t* dta0 = (int16_t*) dta;
//         int16_t* rle0 = (int16_t*) rle;
//         int16_t* decode = (int16_t*) decoded_data;
//         int loc = first_pair_pos; int rle_rest = rle0[loc] + minimum_rle;
//         int write_bias = first_decoded_pos + 1; 
//         if(first_decoded_pos == 0) decode[first_decoded_pos] = start; 
//         else {
//             decode[first_decoded_pos] = (dta0[loc] + minimum_dta);
//             rle_rest --;
//             if(rle_rest == 0) {
//                 loc ++; rle_rest = rle0[loc] + minimum_rle;
//             }
//         }
//         int16_t buffer[16]; int buffer_pos = 0;
//         while(loc < total_len) {
//             buffer[buffer_pos++] = dta0[loc];
//             rle_rest --;
//             if(buffer_pos == 16) {
//                 __m256i* buf = (__m256i*) buffer;
//                 __m256i vecDta = _mm256_loadu_si256(buf);
//                 vecDta = _mm256_add_epi16(vecDta, _mm256_set1_epi16(minimum_dta));
//                 __m256i ps = partialsum16(vecDta);
//                 ps = _mm256_add_epi16(ps, _mm256_set1_epi16(start));
//                 _mm256_storeu_si256((__m256i*)(decode + write_bias), ps);
//                 write_bias += 16;
//                 start = decode[write_bias - 1];
//                 memset(buffer, 0, 16*sizeof(int16_t));
//             }
//             if(rle_rest == 0) {
//                 loc ++; rle_rest = rle0[loc] + minimum_rle;
//             }
//         }
//         if(buffer_pos!=0) {
//             __m256i* buf = (__m256i*) buffer;
//             __m256i vecDta = _mm256_loadu_si256(buf);
//             vecDta = _mm256_add_epi16(vecDta, _mm256_set1_epi16(minimum_dta));
//             __m256i ps = partialsum16(vecDta);
//             ps = _mm256_add_epi16(ps, _mm256_set1_epi16(start));
//             _mm256_storeu_si256((__m256i*)(decode + write_bias), ps);
//             write_bias += buffer_pos;
//             start = decode[write_bias - 1];
//         }
//     }
//     bool* filter(int c1, int c2, int start_pos, int end_pos) {
//         if(start_pos == -1) {start_pos = 0;}
//         if(end_pos == -1) {end_pos=total_len;}
//         bool* res = (bool*) malloc((end_pos - start_pos + 16)*sizeof(bool));
//         for(int i=start_pos;i<end_pos;i++) {
//             if(decoded_data[i] > c1 && decoded_data[i] < c2) {
//                 res[i] = true;
//             }
//         }
//         return res;
//     }
// };

#pragma once
#endif