import random
import math


BP_LEN = {4, 8, 16, 24, 32, 64}

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
    def write_bit():
        buf |= 1
        cur += 1
    
    def skip_bit():
        cur += 1
    
    def check():
        if cur == 31: return True
        else: return False
    
    def export_and_restart():
        cur = 0
        buf = 0
        buffer = buf
        return buffer

class reader:
    cur = 0
    buf = -1
    lth = 0
    def init_rd(bp_len, data):
        lth = bp_len
        buf = data
    
    def read_next():
        if buf & (1 << cur): 
            cur += 1
            return True
        else: 
            cur += 1
            return False
    
    def check():
        if cur == lth: return True
        else: return False

def bitpack(data, lth, aligned=True):
    minimum = min(data[1:])
    maximum = max(data[1:])
    bitlen = math.ceil(math.log2(maximum-minimum))
    if aligned: bitlen = min([x for x in BP_LEN if x > bitlen])
    buffer = [data[0], lth, minimum, bitlen]
    wt = writer()
    rd = reader()
    for i in range(1, lth):
        rd.init_rd(bitlen, data[i])
        while rd.check():
            if rd.read_next():
                wt.write_bit()
            else:
                wt.skip_bit()
            if wt.check():
                buffer.append(wt.export_and_restart())
    return buffer

def write_file(data, dlen, path):
    fop = open(path, "wb")
    for i in range(len(data)):
        fop.write(data[i])
    fop.close()



