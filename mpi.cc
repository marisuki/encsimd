#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <assert.h>
#include "mpp.hpp"
#include <fstream>
#include <map>
#include <set>
#include <list>
#include<queue>


void read_from_encoded_file(std::string path, long bias, long length, int32_t* buffer) {
    int32_t tmp; 
    long count = 0;
    std::ifstream file(path, std::ifstream::in | std::ifstream::binary | std::ifstream::beg);
    while(file.is_open() && !file.eof()) {
        file.read((char*)&tmp, sizeof tmp);
        printf("%d ", tmp);
        if(bias != 0) {
            bias --; continue;
        }
        buffer[count++] = tmp;
        if(count == length) break;
    }
    file.close();
}

// inline int log2(int x) {
//     int cnt = 0; 
//     while(x) {
//         x -= (x & (1 << cnt));
//         cnt ++;
//     }
//     return cnt;
// }

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

struct AttrList {
    int* attr;
    AttrList() {}
    AttrList(int length) {
        attr = (int*) malloc((length + 2)*sizeof(int));
    }
};

struct Label {
    int* relatedAttrs;
    int card;
    std::string labelType;
};

struct BinRegisterLabel {
    int attr1; int attr2;
    // opcode: 0: mul
    int g_op; int result;
    BinRegisterLabel(){}
    BinRegisterLabel(int a1, int a2) {
        if(a1 < a2) {attr1 = a1; attr2 = a2;}
        else {attr2 = a1; attr1 = a2;}
    }
};

struct UnaryRegisterLabel {
    int attr; int id; int result;
    UnaryRegisterLabel() {}
    UnaryRegisterLabel(int attr): attr(attr) {}
};

struct Filter {
    int filter_dim;
    int conj_c1; // > c1
    int conj_c2; // < c2
    Filter() {}
    Filter(int dim, int c1, int c2): filter_dim(dim), conj_c1(c1), conj_c2(c2) {}
};

struct Aggfun {
    int aggdim1; int aggdim2;
    // 0: ret tpval of aggdim[0], 1: opcode g=0 ^2, 2: [0]*[1], 3: ret 1.(cnt)
    int opcode_g; 
    // 0: sum.
    //int opcode_f; 
    Aggfun() {}
    Aggfun(int dim1, int dim2, int op_g): aggdim1(dim1), aggdim2(dim2), opcode_g(op_g) {}
};

struct Header {
    int start; int total_len; 
    int minimum; int bitlen; int cps_len;
    Header(){}
    Header(int st, int datasetsz, int minimum, int blen, int compress_len): 
        start(st), total_len(datasetsz), minimum(minimum), bitlen(blen), cps_len(compress_len) {}
};

struct LabelTask {
    int labels[2]; int tid;
    bool isSeqSc;
    bool isSeqAdd; bool isSeqPS;
    LabelTask() {}
    LabelTask(int tid, int attr1, int attr2, bool seqsc, bool seqadd, bool seqps){
        if(attr1 < attr2)
        {labels[0] = attr1; labels[1] = attr2;}
        else {
            labels[0] = attr2; labels[1] = attr1;
        }
        this->tid = tid;
        isSeqAdd = seqadd;
        isSeqPS = seqps;
        isSeqSc = seqsc;
    }
};

struct Task {
    int attr; int start; int st_pos;
    int st_dec_pos; int tid;
    Task(){}
    Task(int tid, int attr, int start, int st_pos, int st_dec_pos):
        tid(tid), attr(attr), start(start), st_pos(st_pos), st_dec_pos(st_dec_pos) {}
};

struct Schema {
    std::string timestampFiles; // multiple: cpu/t/1, 2, 3...
    std::vector<std::string> attrFiles; // each files common: cpu/1/...
    Header timeDeltaHeader; Header timeRLEHeader;
    Header *attrDeltaheaders; Header *attrRLEheaders;
    int fread_time_bias = 0; 
    int time_counter = 0;
    int* fread_bias;
    int* file_counter;
    int attr_num;
    Schema() {}
    Schema(std::string timestampFiles, std::vector<std::string> attrFiles, int attr_num) {
        this->timestampFiles = timestampFiles;
        this->attrFiles = attrFiles;
        this->attr_num = attr_num;
        this->attrDeltaheaders = (Header *) malloc((attr_num + 4)*(sizeof(Header)));
        this->attrRLEheaders = (Header *) malloc((attr_num + 4)*(sizeof(Header)));
        this->fread_bias = (int*) malloc((attr_num + 4)*(sizeof(int)));
        this->file_counter = (int*) malloc((attr_num + 4)*(sizeof(int)));
        init();
    }
    void init() {
        int32_t buffer[10]; memset(buffer, 0, 10*sizeof(int32_t));

        read_from_encoded_file(current_file(-1, true), fread_time_bias, 5, buffer);
        this->timeDeltaHeader = Header(buffer[0], buffer[1], buffer[2], buffer[3], buffer[4]);
        //printf("Dt %d %d %d %d\n", buffer[0], buffer[1], buffer[2], buffer[3], buffer[4]);
        memset(buffer, 0, 10*sizeof(int32_t));

        //printf("%s\n", current_file(-1, false).c_str());
        read_from_encoded_file(current_file(-1, false), fread_time_bias, 5, buffer);
        this->timeRLEHeader = Header(buffer[0], buffer[1], buffer[2], buffer[3], buffer[4]);
        //printf("Rt %d %d %d %d\n", buffer[0], buffer[1], buffer[2], buffer[3], buffer[4]);
        fread_time_bias += 5;
        for(int i=0;i<this->attr_num;i++) {
            read_from_encoded_file(current_file(i, true), fread_bias[i], 5, buffer);
            this->attrDeltaheaders[i] = Header(buffer[0], buffer[1], buffer[2], buffer[3], buffer[4]);
            read_from_encoded_file(current_file(i, false), fread_bias[i], 5, buffer);
            this->attrRLEheaders[i] = Header(buffer[0], buffer[1], buffer[2], buffer[3], buffer[4]);
            //printf("R %d %d %d %d\n", buffer[0], buffer[1], buffer[2], buffer[3]);
            fread_bias[i] += 5;
        }
    }
    void readFiles(int attr, int &length, int32_t* bufferDelta, int32_t* bufferRLE) {
        int bias = read_upd_bias(attr, length);
        read_from_encoded_file(current_file(attr, true), bias, length, bufferDelta);
        read_from_encoded_file(current_file(attr, false), bias, length, bufferRLE);
    }
    int read_upd_bias(int attr, int &inc) {
        int res = 0;
        if(attr == -1) {
            int enclen = timeDeltaHeader.cps_len - (fread_time_bias - 4);
            if(inc == -1) inc = enclen;
            else inc = std::min(inc, enclen);
            res = fread_time_bias; fread_time_bias += inc;
        } else if (attr > 0 && attr < attr_num) {
            int enclen = attrDeltaheaders[attr].cps_len - (fread_bias[attr] - 4);
            if(inc == -1) inc = enclen;
            else inc = std::min(inc, enclen);
            res = fread_bias[attr]; fread_bias[attr] += inc; 
        }
        return res;
    }
    int curr_rest_len(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.cps_len - (fread_time_bias - 4);
        } else {
            return attrDeltaheaders[attr].cps_len - (fread_time_bias - 4);
        }
    }
    int encoding_len(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.cps_len - 1;
        } else {
            return attrDeltaheaders[attr].cps_len - 1;
        }
    }
    int minimum_delta_(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.minimum;
        } else {
            return attrDeltaheaders[attr].minimum;
        }
    }
    int minimum_rle_(int attr) {
        if(attr == -1) {
            return timeRLEHeader.minimum;
        } else {
            return attrRLEheaders[attr].minimum;
        }
    }
    int bitlen_(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.bitlen;
        } else {
            return attrDeltaheaders[attr].bitlen;
        }
    }
    int start_(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.start;
        } else {
            return attrDeltaheaders[attr].start;
        }
    }
    int total_decoded_len_(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.total_len;
        } else {
            return attrDeltaheaders[attr].total_len;
        }
    }
    int total_compressed_len_(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.cps_len;
        } else {
            return attrDeltaheaders[attr].cps_len;
        }
    }
    std::string current_file(int attr, bool delta) {
        std::string res;
        char bytes[10];
        if(attr == -1) {
            //itoa(time_counter, bytes, 10);
            std::string tmp = timestampFiles;
            if(delta) {//.append(bytes)
                res = tmp.append("t").append("-rle-d.bp0");
            }
            else {
                //res = timestampFiles.append("t").append(bytes).append("-rle-r.b");
                res = tmp.append("t").append("-rle-r.bp0");
            }
        } else {
            itoa(attr + 1, bytes, 10);
            std::string tmp = attrFiles.at(attr);
            if(delta) {
                //res = attrFiles[attr].append("t").append(bytes).append("-rle-d.b");
                res = tmp.append(bytes).append("-rle-d.bp0");
            }
            else {
                //res = attrFiles[attr].append("t").append(bytes).append("-rle-r.b");
                res = tmp.append(bytes).append("-rle-r.bp0");
            }
        }
        return res;
    }
    int cost() {

    }
};

struct Job {
    bool vector_decoding[16];
    bool scalar_decoding[16];
    bool sequential_sc[16];
    UnaryRegisterLabel ul[16]; 
    BinRegisterLabel bl[16];
    Job(){}
};

int* rle_aggr_position_decision(void* data, int bitlen, int elementsz) {
    int* pos = (int32_t*) malloc((elementsz+16) * sizeof(int32_t));
    memset(pos, 0, (elementsz+16) * sizeof(int32_t));
    if(bitlen == 32) {
        int32_t *tmp = (int32_t*) data;
        pos[0] = tmp[0];
        for(int i=1;i<elementsz;i++) {
            pos[i] = pos[i-1] + tmp[i];
        }
    } else if (bitlen == 16) {
        int16_t *tmp = (int16_t*) data;
        pos[0] = tmp[0];
        for(int i=1;i<elementsz;i++) {
            pos[i] = pos[i-1] + tmp[i];
        }
    } else if (bitlen == 8) {
        int8_t *tmp = (int8_t*) data;
        pos[0] = tmp[0];
        for(int i=1;i<elementsz;i++) {
            pos[i] = pos[i-1] + tmp[i];
        }
    } else if (bitlen == 4) {
        char *tmp = (char*) data;
        pos[0] = int(tmp[0]);
        for(int i=1;i<elementsz;i++) {
            pos[i] = pos[i-1] + tmp[i];
        }
    }
    return pos;
}

class TSAcc {
public:
    //int* data; int* delta; int* rle; 
    int* latest; int* result; 
    int datasz; int slicesz;
    Filter filter; 
    int num_decode_dim; int schema_attrs; // exclude timestamp
    int cores; int aggs;
    Schema schema;
    std::set<LabelTask> tasks;
    int num_registers; Job jobs[10]; int num_jobs;
    TSAcc(){}
    TSAcc(int slice_sz, int core, 
            Schema schema, 
            std::set<LabelTask> tasks, 
            const int num_register) {
        this->slicesz = slice_sz;
        this->cores = core;
        this->num_registers = num_register;
        this->tasks = tasks;
        this->schema = schema;
        //delta = (int*) malloc(num_decode_dim*(slice_sz*core+100)*sizeof(int));
        //data = (int*) malloc(num_decode_dim*(slice_sz*core+100)*sizeof(int));
        //rle = (int*) malloc(num_decode_dim*(slice_sz*core+100)*sizeof(int));
        latest = (int*) malloc((cores + 10)*sizeof(int));
        result = (int*) malloc((500)*sizeof(int));
    }
    //void loadDelta(std::string path, long bias, long length);
    //void loadRLE(std::string path, long bias, long length);

    //void sboost_decode();
    //void scalar_decode();

    void schedule(bool naivePlan);
    Job nextJob(int id);
    void monitor();

    std::map<int, std::set<LabelTask>> planning(bool naive);
    //std::set<int> related_attrs();

    void load(int attr, int &b, int* buffdta, int*buffrle);

    void partialsum(int opcode);
    void partialsumvector16(int start, int* read_dta, int*read_rle, 
            int16_t* write, int total, int bitlen); //op: 1
    void partialsumscalar16(int start, int* read_dta, int*read_rle, 
            int16_t *write, int total, int bitlen); //op: 2
    void partialsumscalar32(int start, int* read_dta, int*read_rle, int *write, int total, int bitlen);

    void sequentialaggRLE(); //op: 3

    void run(int pc, int b);
    void run2(int pc, int b);

    void add(int opcode);
    void addv();//op: 4
    void adds();//op: 5

    void predicate();
};

struct sequentialRLE {
    int16_t* dta; int16_t* rle; int loc; int minimum;
    int rle_rest;  int attr;
    int ps; int tmp[10]; //tmp[i] -> attr^i: filtered
    std::set<int> agg_i;
    sequentialRLE(int16_t* dta, int16_t* rle, std::set<int> agg, int minimum) {
        this->dta = dta; this->rle = rle; this-> agg_i = agg;
        memset(tmp, 0, 10*sizeof(int));
        ps = 0; loc = 0; rle_rest = rle[loc];
    }
    void set_start(int start) {ps = start;}
    void skip(int skip) {
        if(skip > rle_rest) {
            ps += (dta[loc] + minimum)*rle_rest;
            int next = skip-rle_rest;
            rle_rest = 0; loc++;
            this->skip(next);
        } else {
            ps += (dta[loc] + minimum)*skip;
            rle_rest -= skip;
        }
    }
    void aggregate(int valid_n) {
        for(int x: agg_i) {
            if (x == 0) {
                tmp[0] += valid_n;
            } else if(x == 1) {
                tmp[0] += (valid_n*(valid_n+1+2*ps)*dta[loc]) / 2;
            } else if(x == 2) {
                tmp[0] += (valid_n*(valid_n+1)*(2*valid_n+1)*(dta[loc]*dta[loc])/6);
            } else{
                // TODO.
            }
        }
        this->skip(valid_n);
    }
};



void TSAcc::load(int attr, int &b, int* buffdta, int*buffrle) {
    schema.readFiles(attr, b, buffdta, buffrle);
}

std::set<int> related_attrs(std::set<LabelTask> task) {
    std::set<int> res;
    for(LabelTask t: task) {
        res.insert(t.labels[0]);
        res.insert(t.labels[1]);
    }
    return res;
}

std::map<int, std::set<LabelTask>> TSAcc::planning(bool naive) {
    std::set<int> attrs = related_attrs(this->tasks);
    std::map<int, std::set<LabelTask>> res;
    std::set<int> decoded;
    //res.insert(std::pair<int, std::set<int>>(1, related_attrs));
    if(tasks.size() + 2 + attrs.size() < this->num_registers) {
        num_jobs = 1;
        res.insert(std::pair<int, std::set<LabelTask>>(1, this->tasks));
        return res;
    } else {
        // 1 turn: filter -> others -> cache filter, cache reusable attrs
        // other turns: if all in cache -> reload, else decode.
        // concentrate ps in one turn.
        std::queue<LabelTask> que; for(LabelTask t: this->tasks) que.push(t);
        std::set<int> considered;
        while(!que.empty()) {
            LabelTask t = que.front(); int tid = t.tid;
            que.pop();
            if(considered.find(tid)!=considered.end()) continue;
            considered.insert(tid);
            // TODO
        }
    }
}

void partialsumscalar16(int start, int* read_dta, int*read_rle, int16_t *write, int total, int bitlen) {
    //int slicebias = (attr+1)*this->slicesz;
    int pos = 1; int last = start;
    for(int i=1;i<total;i++) {
        int delta = read_dta[i];
        int rle = read_rle[i];
        while(rle--) {
            write[pos] = last + delta;
            last = write[pos];
            pos++;
        }
    }
}

void partialsumscalar32(int start, int* read_dta, int*read_rle, int *write, int total, int bitlen) {
    //int slicebias = (attr+1)*this->slicesz;
    int pos = 1; int last = start;
    for(int i=1;i<total;i++) {
        int delta = read_dta[i];
        int rle = read_rle[i];
        while(rle--) {
            write[pos] = last + delta;
            last = write[pos];
            pos++;
        }
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

void partialsumvector16(int start, int* read_dta, int*read_rle, 
            int16_t* write, int total, int bitlen) {
    if(bitlen == 16) {
        int16_t* dta = (int16_t*) read_dta;
        int16_t* rle = (int16_t*) read_rle;
        int loc = 0; int rle_rest = rle[0];
        int write_bias = 0;
        int16_t buffer[16]; int buffer_pos = 0;
        while(loc < total) {
            buffer[buffer_pos++] = dta[loc];
            rle_rest --;
            if(buffer_pos == 16) break;
            if(rle_rest == 0) {
                loc ++; rle_rest = rle[loc];
            }
        }
        __m256i* buf = (__m256i*) buffer;
        __m256i vecDta = _mm256_loadu_si256(buf);
        __m256i ps = partialsum16(vecDta);
        ps = _mm256_add_epi16(ps, _mm256_set1_epi16(start));
        _mm256_storeu_si256((__m256i*)(write + write_bias), ps);
        write_bias += 16;
        start = write[write_bias - 1];
        memset(buffer, 0, 16*sizeof(int16_t));
    }
}

// void partialsumscalar16(int start, int* read_dta, int*read_rle, 
//             int16_t* write, int total, int bitlen) {
//     if(bitlen == 16) {
//         int16_t* dta = (int16_t*) read_dta;
//         int16_t* rle = (int16_t*) read_rle;
//         int loc = 0; int rle_rest = rle[0];
//         int write_bias = 1; write[0] = start + dta[0];
//         while(loc < total) {
//             write[write_bias] = write[write_bias-1] + dta[loc];
//             rle_rest --;
//             if(rle_rest == 0) {
//                 loc ++; rle_rest = rle[loc];
//             }
//         }
//     }
// }



struct DTARLE { // implementation
    int* dta; int* rle; int start; int minimum_dta; int minimum_rle;
    int cps_num; int total_len;
    int* decoded_data; Schema schema; int attr;
    DTARLE(){}
    DTARLE(Schema schema, int attr) {
        this->attr = attr;
        cps_num = schema.total_compressed_len_(attr);
        total_len = schema.total_decoded_len_(attr);
        decoded_data = (int*) malloc((total_len + 16)*sizeof(int));
        dta = (int*) malloc((cps_num+ 16)*sizeof(int));
        rle = (int*) malloc((cps_num+ 16)*sizeof(int));
        schema.readFiles(attr, cps_num, dta, rle);
        start = schema.start_(attr);
        minimum_dta = schema.minimum_delta_(attr);
        minimum_rle = schema.minimum_rle_(attr);
    }

    std::vector<Task> decode_vec(int pc, int b, bool expt=false) {
        if(pc*b > cps_num) {b = cps_num/pc; }
        int dec_loc[100];
        std::vector<Task> exp;
        for(int i=0;i<pc;i++) {for(int j=0;j<b;j++) dec_loc[i+1] += rle[i*b+j];}
        #pragma omp parallel for
        for(int i=0;i<pc;i++) {
            if(!expt) this->decode_vec16(i*b, dec_loc[i]);
            else exp.insert(exp.begin(), Task(i, attr, start, i*b, dec_loc[i]));
        }
        return exp;
    }

    void decode_scalar16() {
        this->decode_scalar16(0, 0);
    }
    
    void decode_scalar16(int first_pair_pos, int first_decoded_pos) {
        int16_t* dta0 = (int16_t*) dta;
        int16_t* rle0 = (int16_t*) rle;
        int16_t* decode = (int16_t*) decoded_data;
        int loc = first_pair_pos; int rle_rest = rle0[loc] + minimum_rle;
        int write_bias = first_decoded_pos + 1; 
        if(first_decoded_pos == 0) decode[first_decoded_pos] = start; 
        else {
            decode[first_decoded_pos] = (dta0[loc] + minimum_dta);
            rle_rest --;
            if(rle_rest == 0) {
                loc ++; rle_rest = rle0[loc] + minimum_rle;
            }
        }
        while(loc < total_len) {
            decode[write_bias] = decode[write_bias-1] + (dta0[loc] + minimum_dta);
            rle_rest --;
            if(rle_rest == 0) {
                loc ++; rle_rest = rle0[loc] + minimum_rle;
            }
        }
    }

    void decode_vec16() {
        this->decode_vec16(0, 0);
    }

    void decode_vec16(int first_pair_pos, int first_decoded_pos) {
        int16_t* dta0 = (int16_t*) dta;
        int16_t* rle0 = (int16_t*) rle;
        int16_t* decode = (int16_t*) decoded_data;
        int loc = first_pair_pos; int rle_rest = rle0[loc] + minimum_rle;
        int write_bias = first_decoded_pos + 1; 
        if(first_decoded_pos == 0) decode[first_decoded_pos] = start; 
        else {
            decode[first_decoded_pos] = (dta0[loc] + minimum_dta);
            rle_rest --;
            if(rle_rest == 0) {
                loc ++; rle_rest = rle0[loc] + minimum_rle;
            }
        }
        int16_t buffer[16]; int buffer_pos = 0;
        while(loc < total_len) {
            buffer[buffer_pos++] = dta0[loc];
            rle_rest --;
            if(buffer_pos == 16) {
                __m256i* buf = (__m256i*) buffer;
                __m256i vecDta = _mm256_loadu_si256(buf);
                vecDta = _mm256_add_epi16(vecDta, _mm256_set1_epi16(minimum_dta));
                __m256i ps = partialsum16(vecDta);
                ps = _mm256_add_epi16(ps, _mm256_set1_epi16(start));
                _mm256_storeu_si256((__m256i*)(decode + write_bias), ps);
                write_bias += 16;
                start = decode[write_bias - 1];
                memset(buffer, 0, 16*sizeof(int16_t));
            }
            if(rle_rest == 0) {
                loc ++; rle_rest = rle0[loc] + minimum_rle;
            }
        }
        if(buffer_pos!=0) {
            __m256i* buf = (__m256i*) buffer;
            __m256i vecDta = _mm256_loadu_si256(buf);
            vecDta = _mm256_add_epi16(vecDta, _mm256_set1_epi16(minimum_dta));
            __m256i ps = partialsum16(vecDta);
            ps = _mm256_add_epi16(ps, _mm256_set1_epi16(start));
            _mm256_storeu_si256((__m256i*)(decode + write_bias), ps);
            write_bias += buffer_pos;
            start = decode[write_bias - 1];
        }
    }
    bool* filter(int c1, int c2, int start_pos, int end_pos) {
        if(start_pos == -1) {start_pos = 0;}
        if(end_pos == -1) {end_pos=total_len;}
        bool* res = (bool*) malloc((end_pos - start_pos + 16)*sizeof(bool));
        for(int i=start_pos;i<end_pos;i++) {
            if(decoded_data[i] > c1 && decoded_data[i] < c2) {
                res[i] = true;
            }
        }
        return res;
    }
};

struct Consume {
    int pos[100]; int rle_pos[100];
    int curr_start[100];
    int total_read[100];
    Consume(){
        memset(pos, 0, sizeof(int)*100);
        memset(rle_pos, 0, sizeof(int)*100);
        memset(curr_start, 0, sizeof(int)*100);
        memset(total_read, 0, sizeof(int)*100);
    }
 };



void TSAcc::run(int pc, int b) {
    std::map<int, std::set<LabelTask>> res = planning(false);
    std::set<int> attrs = related_attrs(res.at(1));
    int timestamp_delta, timestamp_rle, sz=1;
    int timestamp_st = schema.start_(-1);
    this->load(-1, sz, &timestamp_delta, &timestamp_rle);
    //int16_t*data = (int16_t*) malloc(schema.attr_num*(schema.total_decoded_len_(-1)+100)*sizeof(int16_t));
    std::map<int, int*> data; 
    std::map<int, int*> dta; std::map<int, int*> rle;
    for(int attr: attrs) {

    }
}

class etsqp {
    Schema schema; std::map<int, int> label;
    std::map<int, DTARLE> data; int reg_num;
    etsqp() {}
    etsqp(Schema schema, std::map<int, int> label): schema(schema), label(label) {}
    std::vector<Task> plan() {
        
    }
};

// void run(int pc, int b) {
//     std::map<int, std::set<int>> res = planning(false);
//     std::map<int, int> result;
//     std::set<int> related_attrs = res.at(1);
//     int64_t start_time = timeSinceEpochMillisec();
//     int16_t*data = (int16_t*) malloc(schema.attr_num*(schema.total_decoded_len_(-1)+100)*sizeof(int16_t));

//     int timestamp_delta, timestamp_rle, sz=1;
//     int timestamp_st = schema.start_(-1);
//     this->load(-1, sz, &timestamp_delta, &timestamp_rle);
//     int threads = omp_get_num_threads();
//     int each_share = int(schema.total_decoded_len_(-1)/threads);
//     if(filter.conj_c1 > timestamp_st + timestamp_delta*timestamp_rle ||
//         filter.conj_c2 < timestamp_st) return;

//     #pragma omp parallel for // num_threads(pc)
//     for(int attr: related_attrs) { // 4, 8, 16 attrs
//         int*delta = (int*) malloc((schema.total_compressed_len_(attr)+100)*sizeof(int));
//         int*rle = (int*) malloc((schema.total_compressed_len_(attr)+100)*sizeof(int));
//         int readlen = -1;
//         this->load(attr, readlen, delta, rle);
//         if(schema.bitlen_(attr) < 16) {
//             this->partialsumscalar16(schema.start_(attr), delta, rle, data + attr*schema.total_decoded_len_(-1), readlen);
//         }
//     }

//     #pragma vector
//     #pragma omp parallel for
//     for(int slice_i = 0; slice_i < threads; slice_i++) {
//         __m256i r[30];
//         int st_ = each_share*slice_i; int ed_ = std::min(each_share * slice_i, schema.total_decoded_len_(-1));
//         int time_st = timestamp_st + st_*timestamp_delta;
//         int curr_rle = each_share; int skip = 0;
//         while(time_st < filter.conj_c1) {
//             time_st += timestamp_delta; curr_rle --;
//             skip++;
//         }
//         if(time_st >= filter.conj_c1) {
//             for(int ui=0;ui<unary_num;ui++) {
//                 r[ui] = ZERO; result.insert(std::pair<int,int>(ui, 0));
//             }
//             for(int ui=unary_num;ui<unary_num+binary_num;ui++) {
//                 r[ui] = ZERO; result.insert(std::pair<int,int>(ui, 0));
//             } int pos = binary_num + unary_num;
//             while(time_st < filter.conj_c2 && curr_rle > 0) {
//                 for(int attr:related_attrs) {
//                     int16_t* mov = data+attr*(schema.total_decoded_len_(-1)+100);
//                     r[pos+attr] = _mm256_loadu_epi16(mov+skip);
//                 }
//                 time_st += timestamp_delta*16;
//                 curr_rle -= 16; skip += 16;
//                 for(int ui=0;ui<unary_num;ui++) {
//                     r[ui] = _mm256_add_epi16(r[ui], r[pos+unary[ui].attr]);
//                 }
//                 for(int ui=unary_num;ui<unary_num+binary_num;ui++) {
//                     r[ui] = _mm256_add_epi16(r[ui], 
//                         _mm256_mullo_epi16(r[pos+binary[ui-unary_num].attr1], 
//                                     r[pos+binary[ui-unary_num].attr2])
//                     );
//                 }
//             }
//             int16_t buff[16];
//             for(int ui=0;ui<unary_num;ui++) {
//                 _mm256_storeu_epi16(buff,r[ui]); 
//                 int sum = 0;
//                 for(int i=0;i<16;i++) sum += buff[i];
//                 result.insert(std::pair<int,int>(ui, sum));
//             }
//             for(int ui=unary_num;ui<unary_num+binary_num;ui++) {
//                 _mm256_storeu_epi16(buff,r[ui]); 
//                 int sum = 0;
//                 for(int i=0;i<16;i++) sum += buff[i];
//                 result.insert(std::pair<int,int>(ui, sum));
//             }
//         }
//     }
// }

// void TSAcc::run2(int pc, int b) {
//     std::map<int, std::set<int>> res = planning(false);
//     int*data = (int*) malloc(num_decode_dim*(pc*this->slicesz+100)*sizeof(int));
//     int*delta = (int*) malloc(num_decode_dim*(this->slicesz+100)*sizeof(int));
//     int*rle = (int*) malloc(num_decode_dim*(this->slicesz+100)*sizeof(int));
//     Consume monitor = Consume();
//     if(res.size() == 1) {
//         std::set<int> related_attrs = res.at(1);
//         related_attrs.insert(-1);
//         for(int attr: related_attrs) {
//             if(attr == -1) {
//                 //this->load(-1, schema.curr_rest_len(-1));
//                 monitor.curr_start[0] = schema.start_(-1);
//             } else {
//                 int runs = schema.total_decoded_len_(attr)/(pc*b);
//                 monitor.curr_start[attr+1] = schema.start_(attr);
//                 int totrd = std::min(schema.curr_rest_len(attr), this->slicesz);
//                 this->load(attr, totrd);
//                 monitor.total_read[attr+1] = totrd;
//             }
//         }
        
//         //this->load(std::set<int>({-1}), b);
//     int64_t start_time = timeSinceEpochMillisec();
//     #pragma omp parallel for num_threads(pc)
//         for(int attr: related_attrs) {
//             if(attr == -1) { // timestamp, high RLE
                
//             } else { 
//                 // series attributes
//                 if(omp_get_thread_num() == 0) {
//                     this->partialsumscalar(
//                         monitor.curr_start[attr+1], 
//                         delta + omp_get_thread_num()*monitor.total_read[attr+1]/pc,
//                         rle + omp_get_thread_num()*monitor.total_read[attr+1]/pc,
//                         data + omp_get_thread_num()*monitor.)
//                 }
//                 else {

//                 }
//             }
//         }

//     } 
// }




// void TSAcc::schedule(bool naivePlan) {
//     //int num_pack = 0;
//     std::map<int, std::set<int>> pack = TSAcc::planning(naivePlan);
//     if(num_jobs == 1) {
//         // vectorized all temporary
//         // 1 ps of all related dim -> solve all agg
//         jobs[0] = Job();
//         memset(jobs[0].vector_decoding, 1, 16);
//         memset(jobs[0].scalar_decoding, 0, 16);
//         memset(jobs[0].sequential_sc, 0, 16);
//         for(int i=0;i<binary_num;i++) jobs[0].bl[i] = binary[i];
//         for(int i=0;i<unary_num;i++) jobs[0].ul[i] = unary[i];
//     }
    
//     std::set<int> cached;
//     int32_t* cache = (int32_t*) malloc((this->schema.attr_num*64)*(sizeof int32_t));
//     for(int pid = 0; pid<num_jobs; pid++) {
//         std::set<int> labels = pack.at(pid);

//     }
// }

// class MultipleTasks {
//     TSAcc tsacc[8]; int count;
//     MultipleTasks(){}
//     MultipleTasks(TSAcc singleFileTask) {
//         count = 0;
//         addTask(singleFileTask);
//     }
//     void addTask(TSAcc singleFileTask) {
//         tsacc[count++] = singleFileTask;
//     }
//     void run() {
//         for(int i=0;i<count;i++) tsacc[i].schedule(false);
//         int tmp = count; int pos = 0; int lb;
//         while(tmp) {
//             lb = tmp - (tmp & (1 << pos));
//             exe(lb, tmp);
//             tmp = lb;
//             pos ++;
//         }
//     }
//     void exe(int lb, int ub) {
//         for(int jid = lb; jid<ub; jid++) {

//         }
//     }
// }



void query(int qid, int cols, DTARLE time, DTARLE* dim, 
            Schema schema, double select) {
    bool* filter = time.filter(
            schema.start_(-1), 
            time.decoded_data[int(select*schema.total_decoded_len_(-1))],
            0, schema.total_decoded_len_(-1));
    if(qid == 3) { // aggregation average
        for(int i=0;i<cols;i++) {
            average(dim[i].decoded_data, dim[i].total_len, filter);
        }
    } else if(qid == 4) { // downsample
        bool* another = time.filter(
            time.decoded_data[int((1-select)*schema.total_decoded_len_(-1))],
            time.decoded_data[int(schema.total_decoded_len_(-1))],
            0, schema.total_decoded_len_(-1)
        );
        for(int i=0;i<cols;i++) {
            average(dim[i].decoded_data, dim[i].total_len, filter);
            average(dim[i].decoded_data, dim[i].total_len, another);
        }
    } else if(qid == 6) {
        for(int i=0;i<cols;i++) {
            for(int j=0;j<cols;j++) {
                cross_average(
                    dim[i].decoded_data, dim[j].decoded_data,
                    dim[i].total_len, filter
                );
            }
        }
    } else if(qid == 7) {
        for(int i=0;i<cols;i++) {
            for(int j=0;j<cols;j++) {
                cov(
                    dim[i].decoded_data, dim[j].decoded_data,
                    dim[i].total_len, filter
                );
            }
        }
    }
}

void Ours(Schema schema, int qid, int cols, double selectivity, int pc, int b) {

}

void Scalar(Schema schema, int qid, int cols, double selectivity) {
    DTARLE *data = (DTARLE*) malloc((schema.attr_num+2)*sizeof(DTARLE));
    DTARLE time = DTARLE(schema, -1);
    time.decode_scalar16();
    for(int i=0;i<schema.attr_num;i++) {
        data[i] = DTARLE(schema, i);
        data[i].decode_scalar16();
    }
    query(qid, cols, time, data, schema, selectivity);
}


void SBoost(Schema schema, int qid, int cols, double selectivity, int pc) {
    DTARLE *data = (DTARLE*) malloc((schema.attr_num+2)*sizeof(DTARLE));
    DTARLE time = DTARLE(schema, -1);
    time.decode_vec16();
    for(int i=0;i<schema.attr_num;i++) {
        data[i] = DTARLE(schema, i);
        data[i].decode_vec16();
    }
}

int main() {
    std::vector<std::string> arr;
    std::string tmp = "./data/climate/";
    for(int i=1;i<12;i++) {
        std::string t1 = tmp;
        char buff[123];
        t1.append(itoa(i, buff, 10));
        arr.insert(arr.begin()+i-1, t1);
    }
    Schema schema = Schema("./data/climate/", arr, 12);
    for(int i=1;i<10;i++) {
        int64_t time = timeSinceEpochMillisec();
        Scalar(schema, 1, i, 0.5);
        printf("%ld\n", timeSinceEpochMillisec()-time);
    }
    
    // TSAcc acc = TSAcc();
    // std::string attrf[4];
    // for(int i=0;i<4;i++) attrf[i] = "./data/climate/";
    // Schema schema = Schema("./data/climate/", attrf, 4);
    // int32_t buffd[100]; int32_t buffr[100];
    // int readline=20;
    // schema.readFiles(1, readline, buffd, buffr);
    // for(int i=0;i<20;i++) {
    //     printf("%d %d, ", buffd[i], buffr[i]);
    // }
    // printf("\n%d %d\n", schema.timeDeltaHeader.start, schema.timeRLEHeader.bitlen);
    // printf("\n%d %d\n", schema.attrDeltaheaders[1].start, schema.attrDeltaheaders[1].bitlen);
    // //acc.loadDelta("./data/cpu/test-t-rle-r.b", 0, 1000);
    //acc.loadDelta("./data/climate/test-t-rle-d.b", 0, 1000);
}