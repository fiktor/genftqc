import numpy as np
import random
import time
import json


class symplectic_pauli_string:
    def __init__(self, x_comp, z_comp):
        self.x_comp = x_comp
        self.z_comp = z_comp
    

    
def omega(v, w, n):
    v1 = int("".join(str(x) for x in v[:n]), 2)
    v2 = int("".join(str(x) for x in v[n:]), 2)
    
    w1 = int("".join(str(x) for x in w[:n]), 2)
    w2 = int("".join(str(x) for x in w[n:]), 2)
    
    
    return bin((v1 & w2) ^ (w1 & v2)).count('1') & 1



def round_down_pow2(v, n):
    k = 1
    while k < n:
        v |= v >> k
        k *= 2
    v += 1
    return v // 2
    
    
    
    
    
    
    
def generate_random_bitstring(n):
    
    bitstring = []
    
    for i in range(2*n):
        rand = random.random()
        if rand >= 0.5:
            bitstring.append(1)
        else:
            bitstring.append(0)
            
    bitstring = np.array(bitstring)
    
    return bitstring

def bitstring_is_zero(v):
    bitstring_is_zero = True
    
    for bit in v:
        if bit == 1:
            bitstring_is_zero = False
            break
            
    return bitstring_is_zero


def add_bitstrings(bitstring1, bitstring2):
    bitstring_sum = []
    
    for i in range(len(bitstring1)):
        bitstring_sum.append(int((bitstring1[i] + bitstring2[i]) % 2)) 
    
    return np.array(bitstring_sum)


def flip_bit(bitstring, index):
    bitstring[index] = (bitstring[index] + 1) % 2
    return bitstring

def validate_v_w(reference, target, delta):
    if omega(reference, target, int(len(reference) / 2)) == delta:
        return reference, target
    else:
        for i in range(len(reference)):
            if reference[i] == 1:
                target = flip_bit(target, i)
                break
            else:
                continue
    
    return reference, target


                
def validate_v_v(reference, target):
    if omega(reference, target, int(len(reference) / 2)) == 0:
        return reference, target
    else:
        for i in range(len(reference)):
            if reference[i] == 1:
                target = flip_bit(target, i)
                break
            else:
                continue
                
    return reference, target

def validate_w_w(reference, target):
    if omega(reference, target, int(len(reference) / 2)) == 0:
        return reference, target
    else:
        for i in range(len(reference)):
            if reference[i] == 1:
                target = flip_bit(target, i)
                break
            else:
                continue
                
    return reference, target
    

def generate_v(V, W, n):
    v = generate_random_bitstring(n)
    while bitstring_is_zero(v):
        v = generate_random_bitstring(n)
    
    # validate v with previous strings in V and W
    # omega(vi, wj) = delta_ij
    # omega(vi, vj) = 0
    # omega(wi, wj) = 0
    for j in range(len(V)):
        V[j], v = validate_v_v(V[j], v)
        W[j], v = validate_v_w(W[j], v, 0)
        
    # Check that v is not the zero string; if it is, restart the process
    while bitstring_is_zero(v):
        v = generate_v(V, W, n)
        
    return v

def generate_v2(V, W, n):
    v = generate_random_bitstring(n)
    while bitstring_is_zero(v):
        v = generate_random_bitstring(n)
    
    # validate v with previous strings in V and W
    # omega(vi, wj) = delta_ij
    # omega(vi, vj) = 0
    # omega(wi, wj) = 0
    for j in range(len(V)):
        ## Need to work out proper syntax
        # v += omega(v, V[j]) * W[j]
        if omega(v, V[j], n):
            v = add_bitstrings(v, W[j])
        
        # v += omega(v, W[j]) * V[j]
        if omega(v, W[j], n):
            v = add_bitstrings(v, V[j])
        
    # Check that v is not the zero string; if it is, restart the process
    if bitstring_is_zero(v):
        v = generate_v2(V, W, n)
        
    return v


def generate_w(V, W, n):
    w = generate_random_bitstring(n)
    while bitstring_is_zero(w):
        w = generate_random_bitstring(n)
    
    # validate v with previous strings in V and W
    # omega(vi, wj) = delta_ij
    # omega(vi, vj) = 0
    # omega(wi, wj) = 0
    for j in range(len(W)):
        V[j], w = validate_v_w(V[j], w, 0)
        W[j], w = validate_w_w(W[j], w)
        
    last_index = len(W)
    V[last_index], w = validate_v_w(V[last_index], w, 1)
        
    
        
    # Check that w is not the zero string; if it is, restart the process
    while bitstring_is_zero(w):
        w = generate_w(V, W, n)
        
    return w

def fix_w(w, v, n):
    
    v1 = int("".join(str(x) for x in v[:n]), 2)
    v2 = int("".join(str(x) for x in v[n:]), 2)
    
    w1 = int("".join(str(x) for x in w[:n]), 2)
    w2 = int("".join(str(x) for x in w[n:]), 2)
    
    if v1 != 0:
        w2 ^= round_down_pow2(v1, n)
    else:
        assert v2 != 0
        w1 ^= round_down_pow2(v2, n)
            
    
    w1 = f"{w1:0{n}b}"
    w2 = f"{w2:0{n}b}"
    w = w1 + w2
#     print(w1)
#     print(w2)
#     print(w)
    
    new_w = np.array([int(char) for char in w])
     
    assert omega(v, new_w, n)
    
    return new_w


def generate_w2(V, W, n):
    w = generate_random_bitstring(n)
    while bitstring_is_zero(w):
        w = generate_random_bitstring(n)
        
    ## Update this comment # w += (1-omega(w, V[last_index])) * V[last_index]
    last_index = len(V) - 1
    if omega(w, V[last_index], n) == 0:        
        w = fix_w(w, V[last_index], n)
    
    # validate v with previous strings in V and W
    # omega(vi, wj) = delta_ij
    # omega(vi, vj) = 0
    # omega(wi, wj) = 0
    for j in range(len(W)):
        ## Wokring here
        # w += omega(w, V[j]) * W[j]
        if omega(w, V[j], n):
            w = add_bitstrings(w, W[j])
        
        # w += omega(w, W[j]) * V[j]
        if omega(w, W[j], n):
            w = add_bitstrings(w, V[j])
            
    
        
    
        
    # Check that w is not the zero string; if it is, restart the process
    if bitstring_is_zero(w):
        w = generate_w2(V, W, n)
        
    return w



    

def generate_bitstrings(n):
    
    V = []
    W = []
    
    for i in range(n):
        
        v = generate_v2(V, W, n)
        V.append(v)
        
        w = generate_w2(V, W, n)
        W.append(w)
        
    return V, W
        
     
def generate_pauli_string(bitstring, n):
    
    z_comp = bitstring[:n]
    x_comp = bitstring[n:]
    
    pauli_string = ""
    
    for i in range(n):
        if z_comp[i] == 0 and x_comp[i] == 0:
            pauli_string += 'I'
            
        elif z_comp[i] == 0 and x_comp[i] == 1:
            pauli_string += 'X'
        elif z_comp[i] == 1 and x_comp[i] == 1:
            pauli_string += 'Y'
        else:
            pauli_string += 'Z'
            
    return pauli_string
    
    
    
    
    
        
        
def generate_pauli_strings(V, W, n):
    
    V_paulis = []
    W_paulis = []
    
    for bitstring in V:
        pauli_string = generate_pauli_string(bitstring, n)
        V_paulis.append(pauli_string)
        
    for bitstring in W:
        pauli_string = generate_pauli_string(bitstring, n)
        W_paulis.append(pauli_string)
        
    return V_paulis, W_paulis
        
    
    
def generate_code(n, k):
    V, W = generate_bitstrings(n)
    
#     print(V)
#     print()
#     print(W)
    
    # need to create pauli strings out of the bitstring representation
    V_paulis, W_paulis = generate_pauli_strings(V, W, n)
    
#     print(V_paulis)
#     print()
#     print(W_paulis)
    
    # need to use n and k to determine which strings correspond to stabilizers and logical operators
    
    stabilizer_list = []
    logical_ops = []
    
    # adding logical X ops to the list
    for i in range(k):
        logical_ops.append(V_paulis[i])
        
    # adding logical Z ops to the list
    for i in range(k):
        logical_ops.append(W_paulis[i])
        
    for i in range(k, len(V_paulis)):
        stabilizer_list.append(V_paulis[i])
        
        
#     print(f"New [[{n}, {k}, d]] Code:")
#     print()
#     print(f"stab_gen = {stabilizer_list}")
#     print(f"logical_ops = {logical_ops}")
    
    
    # need to return code information in format compatible with genFTQC
    
    code = {
        "stab_gen" : stabilizer_list,
        "logical_ops" : logical_ops
    }
    
    return code


def generate_codes(n, ks, num_codes):
    
    start_time = time.time()
    
    codes = []
    for i in range(num_codes):
        k = random.choice(ks)
        codes.append(generate_code(n, k))
    
    
    end_time = time.time()

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time} seconds")
    
    return codes
        
    
def write_to_json():
    ns = [4, 5, 6, 7, 8, 9, 10, 11, 12]
    ks = [[1, 2], [1, 2], [1, 2], [1, 2, 3], [1, 2, 3], [1, 2, 3], [1, 2, 3, 4], [1, 2, 3, 4], [1, 2, 3, 4]]
    num_codes = [30, 100, 400, 1600, 6000, 20000, 100000, 400000, 471870]
    #num_codes = [3, 3, 3, 3, 3, 3, 3, 3, 3]
    codes = []
    
    for i in range(len(ns)):
        codes += generate_codes(ns[i], ks[i], num_codes[i])
        
        
    with open("codes.json", "w") as f:
        json.dump(codes, f, indent=4)  # indent=4 for pretty formatting
    
    
    
write_to_json()

#generate_codes(9, 1, 10000)
    
#generate_code(5,2)