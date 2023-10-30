#ifndef SCHEMA_H
#define SCHEMA_H

#include "mpp.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <set>
#include <vector>
#include <string>


void read_from_encoded_file(std::string path, long bias, long length, int32_t* buffer) {
    int32_t tmp; 
    long count = 0;
    std::ifstream file(path, std::ifstream::in | std::ifstream::binary | std::ifstream::beg);
    while(file.is_open() && !file.eof()) {
        file.read((char*)&tmp, sizeof tmp);
        //printf("%d ", tmp);
        if(bias != 0) {
            bias --; continue;
        }
        buffer[count++] = tmp;
        if(count == length) break;
    }
    file.close();
}

struct Header {
    int start; int total_len; 
    int minimum; int bitlen; int cps_len;
    Header(){}
    Header(int st, int datasetsz, int minimum, int blen, int compress_len): 
        start(st), total_len(datasetsz), minimum(minimum), bitlen(blen), cps_len(compress_len) {}
};

class Schema {
private:
    std::string timestampFiles; // multiple: cpu/t/1, 2, 3...
    std::vector<std::string> attrFiles; // each files common: cpu/1/...
    Header timeDeltaHeader; Header timeRLEHeader;
    Header *attrDeltaheaders; Header *attrRLEheaders;
    int fread_time_bias = 0; 
    int time_counter = 0;
    int* fread_bias;
    int* file_counter;
    int attr_num;
public:
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
    int get_attr_num() { return this->attr_num; }
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
        //printf("xx\n");
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
    int bitlen_RLE(int attr) {
        if(attr == -1) {
            return timeDeltaHeader.bitlen;
        } else {
            return attrRLEheaders[attr].bitlen;
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
        //char bytes[10];
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
            std::string tmpx = std::to_string(attr + 1);
            //std::itoa(attr + 1, bytes, 10);
            std::string tmp = attrFiles.at(attr);
            if(delta) {
                //res = attrFiles[attr].append("t").append(bytes).append("-rle-d.b");
                res = tmp.append(tmpx).append("-rle-d.bp0");
            }
            else {
                //res = attrFiles[attr].append("t").append(bytes).append("-rle-r.b");
                res = tmp.append(tmpx).append("-rle-r.bp0");
            }
        }
        return res;
    }
};


#pragma once
#endif
