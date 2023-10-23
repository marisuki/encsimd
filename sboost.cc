#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>


template<typename dtaIn, typename rleIn, typename decOut> 
class sboost {
private:
    dtaIn* dta; rleIn* rle; int minimum; int len; int minimumRLE;
    //int attr;
    decOut ps; //tmp[i] -> attr^i: filtered
    decOut* dec;
    void decode_sboost() {
        int pos = 0, loc = 0, curr_dta = dta[pos] + minimum;
        int rle_rest = rle[0] + minimumRLE;
        while(true) {
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
    }
    void apply(int An) {
        for(int i=0;i<loc;i++){
            dec[i] += An;
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
    void decode() {decode_sboost();}
    void finalize(int An) {apply(An);}
};