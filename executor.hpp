#ifndef EXEC_H
#define EXEC_H

#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>
#include <vector>
#include <string>
#include "schema.hpp"
#include "target.hpp"

class Manager {
private:
    __m256i active[16];
    int16_t buffer[20][16];
public:
    Manager(){}
    
};

class Executor {
private:
    Plan plan;
    Schema schema;
    void ps(std::set<int> cols);
public:
    Executor(){}
    Executor(Schema schema) {}
    void set_plan(Plan plan) {this->plan = plan; }
    void decoding_cols(std::set<int> cols);
    void execute();
    void execute1();
};

void Executor::execute() {
    // evaluate timestamp filters.
    Filter tf = plan.time_filter();
    std::set<int> decoded;

    Workset work;
    while(plan.hasNext()) {
        work = plan.next();
        std::set<int> cols = work.required_decCols();
        for(int col: cols) {
            if(decoded.find(col) == decoded.end()) {
                // TODO vector decode
            } else {
                // TODO load vector
            }
            decoded.merge(cols);
        }
        for(Label l: work.ILP_works()) {

        }
    }
}

#pragma once
#endif