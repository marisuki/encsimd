#ifndef TARGET_H
#define TARGET_H

#include <set>
#include <vector>

/**
 * @brief the aggregation target to solve, \sum x^op or \sum x*y
 * op = 0, 1, 2, 3, ...
 * timestamp: 0, attr1-n: 1...n, attr_num = n;
 */
class Label {
private:
    int attr1; int attr2; int op; int id;
public:
    Label() {}
    Label(int id, int attr, int op): id(id), attr1(attr), attr2(-1), op(op) {}
    Label(int id, int attr1, int attr2, int op): id(id), attr1(attr1), attr2(attr2), op(op) {}
    inline bool operator== (const Label &l) const {
        return this->attr1 == l.attr1 && this->attr2 == l.attr2 && op == l.op;
    }
    inline bool operator< (const Label &l) const {
        return this->attr1 < l.attr1 && this->attr2 <= l.attr2 && op <= l.op;
    }
    int attr() {return this->attr1;}
    int oper() {return this->op;}
    int attr_comp() {return this->attr2;}
    int identify() {return this->id; }
};

struct Filter {
    int filter_dim;
    int conj_c1; // > c1
    int conj_c2; // < c2
    Filter() {}
    Filter(int dim, int c1, int c2): filter_dim(dim), conj_c1(c1), conj_c2(c2) {}
};

class Workset {
private:
    std::set<Label> ilp_works; 
    std::set<int> required_cols;
    int agg_method = 0; // 0 decode -> sync -> add -> agg, 1 seqScan -> sync -> add
public:
    Workset() {}
    //Workset(int regNumber) { this->num_regs = regNumber; }
    std::set<Label> ILP_works() { return ilp_works; }
    std::set<int> required_decCols() { return required_cols; }
    void set_cols(std::set<int> required_cols) { this->required_cols = required_cols; }
    void set_agg(bool seqScan) { if(seqScan) agg_method = 1; }
    void insert(Label l) { ilp_works.insert(l); }
    void add_labels(std::set<Label> ls) {ilp_works.merge(ls);}
    bool isSeqSc() {return agg_method == 1;}
    //inline bool check_full() { return ilp_works.size() < num_regs; }
};

class Plan {
private:
    Filter tf;
    std::vector<Workset> works; int pos = 0;
    int iter = 0;
public:
    Plan(){}
    void set_filter(Filter tf) {
        this->tf = tf;
        //works = (Workset*) malloc(int(aggs*1.0/(tot_regs-1.0) + 0.5)*sizeof(Workset));
    }
    void add_workset(Workset work) { works.insert(works.end(), work); }
    Workset next() { return works[iter++]; }
    bool hasNext() { return iter < works.size(); }
    Filter time_filter() { return tf; }
};


#pragma once
#endif