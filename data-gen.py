import random
import math
import struct


BP_LEN = {4, 8, 16, 24, 32}

def bitsz(data):
    return 32*len(data)

def ftoi(data, precision=10):
    return [int(data[i]*precision) for i in range(len(data))]

def delta(data, lth):
    output = [data[0]]
    for i in range(1, lth):
        output.append(data[i]-data[i-1])
    return output

def rle(data, lth):
    output = [data[0]]
    rle = [1]
    for i in range(1, lth):
        if data[i] == output[-1]:
            rle[-1] += 1
        else:
            output.append(data[i])
            rle.append(1)
    return [output, rle]

class writer:
    cur = 0
    buf = 0
    def write_bit(self):
        self.buf << 1
        self.buf |= 1
        self.cur += 1
    
    def skip_bit(self):
        self.buf << 1
        self.cur += 1
    
    def check(self):
        if self.cur == 31: return True
        else: 
            return False
    
    def export_and_restart(self):
        self.cur = 0
        self.buf = 0
        self.buffer = buf
        return self.buffer

class reader:
    cur = 0
    buf = -1
    lth = 0
    def init_rd(self, bp_len, data):
        self.lth = bp_len
        self.buf = data
    
    def read_next(self):
        if self.buf & (1 << self.cur): 
            self.cur += 1
            return True
        else: 
            self.cur += 1
            return False
    
    def check(self):
        if self.cur == self.lth: return True
        else: return False

def aligned_bitpack(data, lth):
    minimum = min(data[1:])
    maximum = max(data[1:])
    if maximum-minimum == 0:
        bitlen = 1
    else: bitlen = math.ceil(math.log2(maximum-minimum))
    buffer = [data[0], lth, minimum, bitlen]
    data = [x-minimum for x in data]
    bitlen = min([x for x in BP_LEN if x > bitlen])
    encsz = int(32/bitlen)
    for i in range(int((lth-1)/encsz)): 
        # 1...8, 9...16
        tmp = 0
        for j in range(encsz):
            tmp |= (data[i*encsz+j+1]<<((encsz-j-1)*bitlen))
        buffer.append(tmp)
    return [len(buffer)+1] + buffer

def bitpack(data, lth):
    minimum = min(data[1:])
    maximum = max(data[1:])
    if maximum-minimum == 0:
        bitlen = 1
    else: bitlen = math.ceil(math.log2(maximum-minimum))
    bitlen = min([x for x in BP_LEN if x > bitlen])
    buffer1 = [data[0], lth, minimum, bitlen, len(data)]
    #wt = writer()
    buffer = []
    pack = int(32/bitlen)
    MAX_VAL = (1 << bitlen) - 1
    for st in range(1, lth):
        pk = 0
        for j in range(0, pack):
            if st + j == lth: break
            pk = pk << bitlen
            pk |= ((data[st + j] - minimum) & (MAX_VAL))
        buffer.append(pk)
    # rd = reader()
    # for i in range(1, lth):
    #     rd.init_rd(bitlen, data[i] - minimum)
    #     while rd.check():
    #         if rd.read_next():
    #             wt.write_bit()
    #         else:
    #             wt.skip_bit()
    #         if wt.check():
    #             buffer.append(wt.export_and_restart())
    return [buffer1, buffer]

def write_file(header, data, path, increase, boost=10):
    fop = open(path, "wb")
    header[0] += increase
    header[2] += increase
    for i in range(len(header)):
        fop.write(header[i].to_bytes(4, byteorder='little', signed=True))
    for x in range(boost):
        for i in range(len(data)):
            fop.write(data[i].to_bytes(4, byteorder='little', signed=True))
    fop.close()

def read_file(path):
    data = []
    with open(path, "rb") as fop:
        byte = fop.read()
        for i in range(1, int(len(byte)/32)+1):
            nxt_i32 = int.from_bytes(byte[(i-1)*32:i*32], byteorder='little', signed=True)
            data.append(nxt_i32)
    return data

def csv(path, skipheader=True, replace1Col=True):
    data = []
    lid = -1
    fop = open(path, "r")
    for line in fop:
        if skipheader and lid == -1:
            lid = 0
        else:
            line = line.split(",")
            if replace1Col: cols = [lid]
            else: cols = [int(line[0])]
            for i in range(1, len(line)):
                cols.append(float(line[i]))
            data.append(cols)
            lid += 1
    #print(data)
    return data

def datagen(pathr, path, extend_col=3, boost=1, extend_partition=2, factor_int=100):
    origin = csv(pathr)
    columnar = []
    for ci in range(len(origin[0])):
        tmp = [x[ci] for x in origin]
        columnar.append(tmp)
        if ci != 0:
            for i in range(1, extend_col):
                pos = int(len(tmp)/(1<<i))
                columnar.append(tmp[pos:] + tmp[:pos])
    print(len(columnar), len(origin[0]))
    for ci in range(len(columnar)):
        exp = columnar[ci]
        increase = int((exp[-1] - exp[0])*factor_int)
        if ci == 0:
            dt = delta(exp, len(exp))
            runrle = rle(dt, len(dt))
            #print(runrle)
            bpd = [bitpack(runrle[0], len(runrle[0])), 
                   bitpack(runrle[1], len(runrle[1]))]
            #print(bpd)
            write_file(bpd[0][0], bpd[0][1], path+"t-rle-d.bp0", 0, boost)
            write_file(bpd[1][0], bpd[1][1], path+"t-rle-r.bp0", 0, boost)
            for partition in range(1, extend_partition):
                write_file(bpd[0][0], bpd[0][1], path+"t-rle-d.bp"+str(partition), increase*(partition), boost)
                write_file(bpd[1][0], bpd[1][1], path+"t-rle-r.bp"+str(partition), 0, boost)
            #bpd = bitpack(dt, len(dt))
            #print(bpd)
            #write_file(bpd[0], bpd[1], path+"t-donly-d.b", boost)
            #break
        else:
            exp = ftoi(exp, factor_int)
            dt = delta(exp, len(exp))
            runrle = rle(dt, len(dt))
            bpd = [bitpack(runrle[0], len(runrle[0])), 
                   bitpack(runrle[1], len(runrle[1]))]
            write_file(bpd[0][0], bpd[0][1], path + str(ci) +"-rle-d.bp0", 0, boost)
            write_file(bpd[1][0], bpd[1][1], path + str(ci) +"-rle-r.bp0", 0, boost)
            for partition in range(1, extend_partition):
                write_file(bpd[0][0], bpd[0][1], path + str(ci) +"-rle-d.bp"+str(partition), increase*(partition), boost)
                write_file(bpd[1][0], bpd[1][1], path + str(ci) +"-rle-r.bp"+str(partition), 0, boost)
    return

if __name__ == "__main__":
    #datagen("./data/cpu.txt", "./data/cpu1p/", extend_col=3, boost=1,extend_partition=1, factor_int=100)
    #datagen("./data/cpu.txt", "./data/cpu1M/", extend_col=3, boost=250, extend_partition=1, factor_int=100)
    datagen("./data/iot.climate.csv", "./data/climate1p/", extend_col=3, boost=1, extend_partition=1, factor_int=100)
    #print(read_file("./data/gen50kd.bn"))
        