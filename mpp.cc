#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include <fstream>
#include <map>
#include <set>
#include <list>
#include <queue>
#include "mpp.hpp"
#include "sboost.hpp"
#include "scalar.hpp"
#include "schema.hpp"
#include "seqRLE.hpp"
#include "aggregation.hpp"
#include "cost.hpp"
#include "target.hpp"

class PlanGenerator {
private:
    std::vector<Label> labels;
    int tot_reg;
    Schema schema;
    Plan plan;
    Filter time_range, tf;
public:
    PlanGenerator() {}
    PlanGenerator(int total_registers, Schema schema, std::vector<Label> aggs): 
        tot_reg(total_registers), schema(schema), labels(aggs) {}
    void set_range(Filter tf, Filter time_range_meta) {
        this->tf = tf; this->time_range = time_range_meta;
    }
    Plan generate();
};

struct node {
    int attr; int cnt; int deg;
    node(){}
    node(int attrx, int cntx, int deg): attr(attrx), cnt(cntx), deg(deg) {}
    void minus(int val) {this->cnt -= val;}
    inline bool operator<(const node &n) const {
        return deg < n.deg && cnt > n.cnt; // heuristic rule.
    }
};


Plan PlanGenerator::generate() {
    Plan plan;
    int attr_num = schema.get_attr_num();
    std::map<int, std::vector<int> > self; 
    std::map<int, std::vector<int> > outer;
    // if(this->labels.size() < this->tot_reg - 1) {
    //     Workset ws = Workset(); //ws.add_labels(std::set<Label>(this->labels));
    //     plan.add_workset(ws);
    // }
    for(int i=0;i<this->labels.size();i++) {
        Label lab = this->labels[i];
        if(lab.attr_comp() == -1) {
            if(self.find(lab.attr())!= self.end()) {
                std::map<int, std::vector<int> >::iterator p = self.find(lab.attr());
                p->second.insert((p->second).end(), i);
            } else {
                std::vector<int> v; v.insert(v.end(), i);
                self.insert(self.end(), std::pair<int, std::vector<int> >(lab.attr(), v));
            }
        } else {
            if(outer.find(lab.attr())!= outer.end()) {
                std::map<int, std::vector<int> >::iterator p = outer.find(lab.attr());
                p->second.insert((p->second).end(), i);
            } else {
                std::vector<int> v; v.insert(v.end(), i);
                outer.insert(outer.end(), std::pair<int, std::vector<int> >(lab.attr(), v));
            }

            if(outer.find(lab.attr_comp())!= outer.end()) {
                std::map<int, std::vector<int> >::iterator p = outer.find(lab.attr_comp());
                p->second.insert((p->second).end(), i);
            } else {
                std::vector<int> v; v.insert(v.end(), i);
                outer.insert(outer.end(), std::pair<int, std::vector<int> >(lab.attr_comp(), v));
            }
        }
    }
    for(int i=1;i<=attr_num;i++){
        if(outer.find(i) == outer.end()) {
            
        }
    }
}

// Plan PlanGenerator::generate() {
//     int attr_num = schema.get_attr_num();
//     int stat[attr_num + 4][attr_num + 4]; 
//     int sum[attr_num+4][2];
//     memset(stat, 0, (attr_num + 4)*(attr_num + 4)*sizeof(int));
//     memset(sum, 0, sizeof(sum));
//     for(auto lab: this->labels) {
//         if(lab.attr_comp() != -1) {
//             stat[lab.attr()][lab.attr_comp()] |= lab.oper();
//             stat[lab.attr_comp()][lab.attr()] |= lab.oper(); // sum x^2, sum xy
//             sum[lab.attr()][0] += 1;
//             sum[lab.attr_comp()][0] += 1;sum[lab.attr_comp()][1] += 1;
//         } else {
//             //if(stat[lab.attr()][lab.attr()]==0) sum[lab.attr()][1] += 1;
//             stat[lab.attr()][lab.attr()] |= (1 << lab.oper()); // sum x
//             sum[lab.attr()][0] += 1;
//         }
//     }
//     std::set<int> decoded;
//     std::list<node> rest;
//     //std::queue<std::pair<int,int>> rest;
//     const int available_regs = tot_reg - 1;
//     for(int i=1;i<=attr_num;i++) rest.insert(rest.end(), node(i, sum[i][0],sum[i][1]));
//     int label_id = 0;
//     while(!rest.empty()) {
//         rest.sort();
//         int rest_reg = available_regs;
//         std::list<node>::iterator it = rest.begin();
//         Workset work(available_regs);
//         while(it!=rest.end()) {
//             int curr = it->attr; decoded.insert(curr);
//             int tmp = 0;
//             while(stat[curr][curr]) {
//                 if(stat[curr][curr]&1) {
//                     work.insert(Label(label_id++, curr, tmp));
//                 }
//                 tmp++; stat[curr][curr] >>= 1;
//                 rest_reg --;
//             }
//             if(rest_reg > it->cnt) {
//                 rest_reg -= it->cnt; sum[curr] = 0;
//                 for(int i=1;i<=attr_num;i++) {
//                     if(stat[curr][i]) {
//                         if(i != curr) {
//                             work.insert(Label(label_id++, curr, i));
//                         }
//                         stat[curr][i] = 0; stat[i][curr] = 0;
//                         decoded.insert(i);
//                     }
//                 }
//             } else {
//                 for(int i=1;i<=attr_num && rest_reg;i++) {
//                     if(decoded.find(i)!=decoded.end()) {
//                         work.insert(Label(label_id++, curr, i));
//                         rest_reg--; stat[curr][i]
//                     }
//                 }
//             }
//         }
//     }
// }


int main() {
    std::list<node> lab;
    // lab.insert(Label(1, 1, 0));
    // lab.insert(Label(2, 2, 0));
    // //lab.insert(BinaryLabel(1, 2, 0));
    // lab.insert(Label(3, 2, 0));
    lab.insert(lab.end(), node(1, 5, 1));
    lab.insert(lab.end(), node(3, 4, 2));
    lab.sort();
    std::list<node>::iterator it = lab.begin();

    std::map<int, std::vector<int> > self;
    std::vector<int> v; v.insert(v.end(), 1);
    self.insert(std::pair<int, std::vector<int> >(1, v));
    std::map<int, std::vector<int> >::iterator p = self.find(1);
    (*p).second.insert((*p).second.end(), 2);
    printf("%d\n",  self[1].size());
    // it++;
    // printf("%d\n", it->cnt);
    return 0;
}