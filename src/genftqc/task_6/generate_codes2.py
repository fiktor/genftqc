import numpy as np
import random

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

    

def generate_pauli_strings(n):
    
    V = []
    W = []
    
    for i in range(n):
        
        v = generate_v(V, W, n)
        V.append(v)
        
        w = generate_w(V, W, n)
        W.append(w)
        
    return V, W
        
        
            
        
    
    
def generate_code(n, k):
    
    V, W = generate_pauli_strings(n)
    
    # need to use n and k to determine which strings correspond to stabilizers and logical operators
    
    # need to create pauli strings out of the bitstring representation
    
    # need to return code information in format compatible with genFTQC

