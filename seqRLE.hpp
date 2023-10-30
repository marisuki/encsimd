#ifndef SEQRLE_H
#define SEQRLE_H


#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>


template<typename dtaIn, typename rleIn> class sequentialRLE {
private:
    dtaIn* dta; rleIn* rle; int loc; int minimum; int len; int minimumRLE;
    long rle_rest; int attr;
    long ps; long tmp[10]; //tmp[i] -> attr^i: filtered
    long res[10];
    std::set<long> agg_i;
    
    void set_start(int start) {ps = start;}
    void aggregate_func(int valid_n) {
        for(int x: agg_i) {
            if (x == 0) {
                tmp[0] += valid_n;
            } else if(x == 1) {
                tmp[1] += (valid_n*(valid_n+1+2*ps)*dta[loc]) / 2;
            } else if(x == 2) {
                long added = ps + valid_n * dta[loc];
                long res1 = (added*(added+1)*(2*added+1))/6;
                long res2 = (ps*(ps+1)*(2*ps+1))/6;
                tmp[2] += res1 - res2;
            } else{
                // TODO.
            }
        }
        //this->skip(valid_n);
    }
public:
    sequentialRLE(dtaIn* dta, rleIn* rle, std::set<int> agg, int mDta, int mRLE, int len);
    sequentialRLE(dtaIn* dta, rleIn* rle, std::set<int> agg, int mDta, int mRLE, int len, int start);
    void skip(int skip);
    void aggregate(int valid_n);
    void finalize(int An);
    long extract(int agg_x);
};

template<typename dtaIn, typename rleIn> 
sequentialRLE<dtaIn, rleIn>::sequentialRLE(dtaIn* dta, rleIn* rle, std::set<int> agg, int mDta, int mRLE, int len) {
    this->dta = dta; this->rle = rle; this-> agg_i = agg;
    memset(tmp, 0, 10*sizeof(long));
    memset(res, 0, 10*sizeof(long));
    ps = 0; loc = 0; rle_rest = (long) rle[loc];
    this->agg_i.insert(0);
    if(this->agg_i.find(2) != this->agg_i.end()) this->agg_i.insert(1);
    this->len = len;
    this->minimumRLE = mRLE;
    this->minimum = mDta;
}

template<typename dtaIn, typename rleIn> 
sequentialRLE<dtaIn, rleIn>::sequentialRLE(dtaIn* dta, rleIn* rle, std::set<int> agg, int mDta, int mRLE, int len, int start) {
    this->dta = dta; this->rle = rle; this-> agg_i = agg;
    memset(tmp, 0, 10*sizeof(long));
    memset(res, 0, 10*sizeof(long));
    loc = 0; rle_rest = (long) rle[loc];
    set_start(start);
    this->agg_i.insert(0);
    if(this->agg_i.find(2) != this->agg_i.end()) this->agg_i.insert(1);
    this->len = len;
    this->minimumRLE = mRLE;
    this->minimum = mDta;
}

template<typename dtaIn, typename rleIn> 
void sequentialRLE<dtaIn, rleIn>::skip(int skip) {
    if(skip > rle_rest) {
        ps += (dta[loc] + minimum)*rle_rest;
        int next = skip-rle_rest;
        loc++; 
        if(loc == len) return;
        rle_rest = (int) rle[loc] + minimumRLE; 
        this->skip(next);
    } else {
        ps += (dta[loc] + minimum)*skip;
        rle_rest -= skip;
    }
}

template<typename dtaIn, typename rleIn> 
void sequentialRLE<dtaIn, rleIn>::aggregate(int valid_n) {
    if(valid_n >= rle_rest) {
        aggregate_func(rle_rest);
        ps += (dta[loc] + minimum)*rle_rest;
        int next = valid_n-rle_rest;
        loc++; 
        if(loc == len) return;
        rle_rest = (int) rle[loc] + minimumRLE; 
        this->aggregate(next);
    } else {
        aggregate_func(valid_n);
        ps += (dta[loc] + minimum)*valid_n;
        rle_rest -= valid_n;
    }
}

template<typename dtaIn, typename rleIn> 
void sequentialRLE<dtaIn, rleIn>::finalize(int An) {
    for(int x: agg_i) {
        if (x == 0) { this->res[0] = this->tmp[0]; }
        if (x == 1) {
            this->res[x] += this->tmp[0]*An;
        } else if(x == 2) {
            this->res[x] += this->tmp[0]*An*An + 2*An*tmp[1];
        } else {
            // TODO
        }
    }
}

template<typename dtaIn, typename rleIn> 
long sequentialRLE<dtaIn, rleIn>::extract(int agg_x) {
    return this->res[agg_x];
}

// int main() {
//     int16_t delta[100];
//     int16_t rle[100];
//     for(int i=0;i<100;i++) {
//         delta[i] = (int16_t) random_i32();
//         rle[i] = std::max((int16_t) 1, (int16_t) random_i32());
//     }
//     std::set<int> aggfun;
//     aggfun.insert(0); aggfun.insert(1); aggfun.insert(2);  
//     sequentialRLE<int16_t, int16_t> seqsc(delta, rle, aggfun, 20, 1, 100);
//     for(int i=0;i<200;i++) {
//         int pos = random_i32()%16;
//         int valid = random_i32()%16;

//     }
//     return 0;
// }

#pragma once
#endif