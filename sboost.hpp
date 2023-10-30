#ifndef SBOOST_H
#define SBOOST_H


#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>

//__m256i mask_pos_epi32_avx2[8];
//__m256i mask_pos_epi16_avx2[16];

__m256i mask_pos_epi32_avx2(int pos) {
    switch(pos) {
        case 1: return  _mm256_set_epi32(-1, 0, 0, 0, 0, 0, 0, 0);
        case 2: return  _mm256_set_epi32(-1, -1, 0, 0, 0, 0, 0, 0);
        case 3: return  _mm256_set_epi32(-1, -1, -1, 0, 0, 0, 0, 0);
        case 4: return  _mm256_set_epi32(-1, -1, -1, -1, 0, 0, 0, 0);
        case 5: return  _mm256_set_epi32(-1, -1, -1, -1, -1, 0, 0, 0);
        case 6: return  _mm256_set_epi32(-1, -1, -1, -1, -1, -1, 0, 0);
        case 7: return  _mm256_set_epi32(-1, -1, -1, -1, -1, -1, -1, 0);
        default: return _mm256_set_epi32(-1, -1, -1, -1, -1, -1, -1, -1);
    }
}

__m256i mask_pos_epi16_avx2(int pos) {
    switch(pos) {
        case 1: return  _mm256_set_epi16(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        case 2: return  _mm256_set_epi16(-1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        case 3: return  _mm256_set_epi16(-1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        case 4: return  _mm256_set_epi16(-1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        case 5: return  _mm256_set_epi16(-1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        case 6: return  _mm256_set_epi16(-1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        case 7: return  _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        case 8: return  _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0);
        case 9: return  _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0);
        case 10: return _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0);
        case 11: return _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0);
        case 12: return _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0);
        case 13: return _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0);
        case 14: return _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0);
        case 15: return _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0);
        default: return _mm256_set_epi16(-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    }
}

__m256i partialsum16(__m256i b) {
    __m256i bp = _mm256_bslli_epi128(b, 2);
    __m256i s1 = _mm256_hadd_epi16(b, bp);
    __m256i s2 = _mm256_sllv_epi64(s1, _mm256_set1_epi64x(16));
    __m256i s3 = _mm256_hadd_epi16(s1, s2);
    __m256i s4 = _mm256_and_si256(s3, _mm256_set1_epi32(0xffff));
    __m256i result = _mm256_hadd_epi16(s3, s4);
    __m256i inv = _mm256_setr_epi16(0xF0E, 0xD0C, 0xB0A,
0x908, 0x706, 0x504, 0x302, 0x100, 0xF0E, 0xD0C, 0xB0A, 0x908,
0x706, 0x504, 0x302, 0x100);
    return _mm256_shuffle_epi8(result, inv);
}

__m256i partialsum32(__m256i b) {
    __m256i bp = permute1(b);
    __m256i s1 = _mm256_add_epi32(b, bp);
    __m256i s2 = permute2(s1);
    __m256i s3 = _mm256_add_epi32(s1, s2);
    __m256i s4 = permute4(s3);
    __m256i result = _mm256_add_epi32(s3, s4);
    return result;
}


template<typename dtaIn, typename rleIn, typename decOut> 
class sboost {
private:
    dtaIn* dta; rleIn* rle; int minimum; int len; int minimumRLE;
    int p_inst=-1;
    //int attr;
    decOut ps; //tmp[i] -> attr^i: filtered
    decOut* dec; decOut An; __m256i Anv;
    __m256i_u* dec_buffer;
    void decode_sboost() {
        int pos = 0, loc = 0;//, curr_dta = dta[pos] + minimum;
        int rle_rest = rle[pos] + minimumRLE;
        decOut* buffer = (decOut*) malloc(8*sizeof(int32_t));
        memset(buffer, 0, sizeof(buffer));
        int buffer_loc = 0;
        if(p_inst == -1) p_inst == 8;
        decOut curr = (decOut) dta[pos] + minimum;
        __m256i x; decOut an = 0;
        __m256i_u* dec_buffer = (__m256i_u*) dec;
        while(true) {
            if(buffer_loc == 0 && rle_rest >= p_inst) {
                if(p_inst == 16) {
                    x = _mm256_set1_epi16(curr); x = partialsum16(x);
                    x = _mm256_add_epi16(x, _mm256_set1_epi16(an));
                    _mm256_storeu_si256(dec_buffer, x);
                    an = dec[loc + 15];
                    dec_buffer ++; loc += 16;
                }
                else {
                    x = _mm256_set1_epi32(curr); x = partialsum32(x);
                    x = _mm256_add_epi32(x, _mm256_set1_epi32(an));
                    _mm256_storeu_si256(dec_buffer, x);
                    an = dec[loc + 7];
                    dec_buffer ++; loc += 8;
                }
                rle_rest -= p_inst;
            } else if (buffer_loc + rle_rest < p_inst) {
                while(rle_rest--) {
                    buffer[buffer_loc++] = curr;
                }
            } else {
                x = _mm256_loadu_si256(buffer);
                int tmp = 0;
                if(p_inst==16) {
                    tmp = 16 - buffer_loc;
                    x = _mm256_and_si256(x, 
                        _mm256_and_si256(mask_pos_epi16_avx2(tmp), 
                                            _mm256_set1_epi16(curr)));
                    x = partialsum16(x);
                    x = _mm256_add_epi16(x, _mm256_set1_epi16(an));
                    _mm256_storeu_si256(dec_buffer, x);
                    an = dec[loc + 15];
                    dec_buffer ++; loc += 16; buffer_loc = 0;
                    memset(buffer, 0, sizeof(buffer));
                } else {
                    tmp = 32 - buffer_loc;
                    x = _mm256_and_si256(x, 
                        _mm256_and_si256(mask_pos_epi32_avx2(tmp),
                                            _mm256_set1_epi32(curr)));
                    x = partialsum32(x);
                    x = _mm256_add_epi32(x, _mm256_set1_epi32(an));
                    _mm256_storeu_si256(dec_buffer, x);
                    an = dec[loc + 7];
                    dec_buffer ++; loc += 8; buffer_loc = 0;
                    memset(buffer, 0, sizeof(buffer));
                }
                rle_rest -= tmp;
            }
            //dec[loc+1] = dec[loc] + curr_dta;
            //loc++; rle_rest--;
            if(rle_rest == 0) {
                pos++;
                if(pos < len) {
                    rle_rest = rle[pos] + minimumRLE; 
                    curr = dta[pos] + minimum; 
                } else {
                    break;
                }
            }
        }
        if(buffer_loc != 0) {
            for(int i=0;i<buffer_loc;i++) {
                this->dec[loc+1] = this->dec[loc] + buffer[i]; loc ++;
            }
        }
    }
    void apply(int An) {
        __m256i_u* dec_buffer = (__m256i_u*) dec;
        __m256i x;
        for(int i=0;i<int(len/p_inst);i++){
            if(p_inst==16) {
                x = _mm256_loadu_epi16(dec_buffer);
                x = _mm256_add_epi16(x, _mm256_set1_epi16(An));
                _mm256_storeu_epi16(dec_buffer, x);
                dec_buffer++;
            }
        }
        if(this->len%p_inst) {
            for(int i=this->len-this->len%p_inst;i<this->len;i++) {
                dec[i] += An;
            }
        }
    }
public:
    sboost(dtaIn* delta, rleIn* rle, decOut* decBuff, 
            decOut start, int mDta, int mRLE, int len) {
        this->dta = delta;
        this->rle = rle;
        this->dec = decBuff;
        this->ps = start;
        this->minimumRLE = mRLE;
        this->minimum = mDta;
        this->len = len;
    }
    void set_pinst(int p_inst) {this->p_inst = p_inst; }
    void decode() {decode_sboost();}
    void finalize(int An) {apply(An);}
    void set_streaming(decOut An) {
        this->dec_buffer = (__m256i_u*) dec;
        this->An = An; this->Anv = _mm256_set1_epi16(An);
    }
    __m256i next() {
        __m256i x = _mm256_loadu_epi16(dec_buffer);
        return _mm256_add_epi16(x, Anv);
    }
};

class sboost_naive {

};

#pragma once
#endif