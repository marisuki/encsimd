#ifndef AGGREG_H
#define AGGREG_H

#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>

inline double average16(int16_t* data, int len, bool* filter) {
    int res = 0; int cnt = 0;
    for(int i=0;i<len;i++) {
        if(filter[i]) {res += data[i]; cnt += 1;}
    }
    return ((double)res)/((double) cnt);
}

inline double average(int* data, int len, bool* filter) {
    int res = 0; int cnt = 0;
    for(int i=0;i<len;i++) {
        if(filter[i]) {res += data[i]; cnt += 1;}
    }
    return ((double)res)/((double) cnt);
}

inline double cross_average16(int16_t* data1,int16_t* data2, int len, bool* filter) {
    int res = 0; int cnt = 0;
    for(int i=0;i<len;i++) {
        if(filter[i]) {res += (data1[i] + data2[i])/2; cnt += 1;}
    }
    return -1;
}

inline double cross_average(int* data1,int* data2, int len, bool* filter) {
    int res = 0; int cnt = 0;
    for(int i=0;i<len;i++) {
        if(filter[i]) {res += (data1[i] + data2[i])/2; cnt += 1;}
    }
    return -1;
}

inline double cov16(int16_t* data1, int16_t* data2, int len, bool* filter) {
    int avg1 = 0, avg2 = 0, sq1 = 0, sq2 = 0, ab = 0; int cnt = 0;
    for(int i=0;i<len;i++) {
        if(filter[i]) {
            avg1 += data1[i]; avg2 += data2[i];
            sq1 += data1[i]*data1[i]; sq2 += data2[i]*data2[i];
            ab += data1[i]*data2[i];
            cnt += 1;
        }
    }
    return ((double)(cnt*ab - avg1*avg2))/
            (std::sqrt(cnt*sq1-avg1*avg1)*std::sqrt(cnt*sq2-avg2*avg2));
}

inline double cov(int* data1, int* data2, int len, bool* filter) {
    int avg1 = 0, avg2 = 0, sq1 = 0, sq2 = 0, ab = 0; int cnt = 0;
    for(int i=0;i<len;i++) {
        if(filter[i]) {
            avg1 += data1[i]; avg2 += data2[i];
            sq1 += data1[i]*data1[i]; sq2 += data2[i]*data2[i];
            ab += data1[i]*data2[i];
            cnt += 1;
        }
    }
    return ((double)(cnt*ab - avg1*avg2))/
            (std::sqrt(cnt*sq1-avg1*avg1)*std::sqrt(cnt*sq2-avg2*avg2));
}

#pragma once
#endif