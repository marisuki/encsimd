#ifndef COST_H
#define COST_H

#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>
#include <vector>
#include <string>
#include "schema.hpp"
#include "target.hpp"


enum Cx {
    C_mul = 5, C_add = 1, C_div = 11, C_ext = 1, C_ld = 1, C_shuffle = 1, C_and = 1
};

// used to decide seqsc or vagg.
struct Cost {
    Schema schema; Filter tf; Filter range;
    int p_inst; int p_core; int N;
    int RLE0; int l_rle;
    Cost(){}
    Cost(int size, Schema schema, Filter tf, Filter range): N(size), schema(schema), tf(tf), range(range) {}
    double cost_seqsc(int attr) {
        RLE0 = schema.minimum_rle_(attr);
        l_rle = schema.bitlen_RLE(attr);
        double cost = 0, selectivity = (tf.conj_c2 - tf.conj_c1)*1.0/(range.conj_c2 - range.conj_c1);
        double expect_RLE = RLE0 + (1 << (l_rle - 1));
        return (N*selectivity*(3*Cx(C_mul) + 4*Cx(C_add) + Cx(C_div)))/expect_RLE + (N*(1.0-selectivity)*(Cx(C_mul) + 2*Cx(C_add)))/expect_RLE; 
    }
    double cost_vagg(int attr) {
        double ps = N*1.0*(3*Cx(C_add) + 4*Cx(C_shuffle) + Cx(C_and) + Cx(C_ld) + Cx(C_ext))/16.0;
        double dec = N*1.0*(Cx(C_ld) + Cx(C_ext))/p_inst;
        double agg_additional = N*1.0*(Cx(C_ext) + Cx(C_add))/p_inst;
        return ps + dec + agg_additional;
    }
};


#pragma once
#endif