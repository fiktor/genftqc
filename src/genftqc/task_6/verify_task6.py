import numpy as np
import time


## Define Pauli Logic
def create_symplectic_representation(string):
    
    symplectic_string_x = []
    symplectic_string_z = []
    
    for char in string:
        if char == 'X':
            symplectic_string_x.append(1)
            symplectic_string_z.append(0)
        elif char == 'Y':
            symplectic_string_x.append(1)
            symplectic_string_z.append(1)
        elif char == 'Z':
            symplectic_string_x.append(0)
            symplectic_string_z.append(1)
        else:
            symplectic_string_x.append(0)
            symplectic_string_z.append(0)
            
    symplectic_string = symplectic_string_z + symplectic_string_x
    
    return np.array(symplectic_string)

def generate_pauli_string(bitstring):
    n = int(len(bitstring) / 2)
    
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

# To multiply Pauli Strings we simply create symplectic representations and add the strings using modulo 2 addition
def add_bitstrings(bitstring1, bitstring2):
    bitstring_sum = []
    
    for i in range(len(bitstring1)):
        bitstring_sum.append(int((bitstring1[i] + bitstring2[i]) % 2)) 
    
    return np.array(bitstring_sum)  


def add_multiple_bitstrings(bitstrings, n):
    
    temp = np.zeros(2*n)
    
    for bitstring in bitstrings:
        temp = add_bitstrings(temp, bitstring)
        
    return temp


def select_generators(generator_list, number):
    
    bits = []
    for i, c in enumerate(bin(number)[:1:-1], 1):
        if c == '1':
            bits.append(i-1)
    #print(bits)
    generators = []
    for index in bits:
        generators.append(generator_list[index])
        
    return generators

# generator_list = ['XXII', 'IIXX', 'ZZII', 'IIZZ']
# number = 15

# print(select_generators(generator_list, number))


def gen_logical_op_set(logical_op, stabilizer_generators):
    n = len(logical_op)
    
    symplectic_alternative_reps = []
    
    syplectic_logical_op = create_symplectic_representation(logical_op)
    
    symplectic_stabilizer_generators = []
    
    for stabilizer in stabilizer_generators:
        symplectic_stabilizer_generators.append(create_symplectic_representation(stabilizer))
        
    num_stabilizers = 2**len(stabilizer_generators)
    
    ## generate entire stabilizer set
    symplectic_stabilizer_set = []
    for i in range(num_stabilizers):
        # add stabilizer generators corresponding to index with one in the binary form to a list such that they can be multiplied together to generate a stabilizer
        generators = select_generators(symplectic_stabilizer_generators, i)
        symplectic_stabilizer_set.append(add_multiple_bitstrings(generators, n))
        
    for stabilizer in symplectic_stabilizer_set:
        symplectic_alternative_reps.append(add_bitstrings(stabilizer, syplectic_logical_op))
        
        
    alternative_reps = []
    
    for alternative_rep in symplectic_alternative_reps:
        alternative_reps.append(generate_pauli_string(alternative_rep))
        
    alternative_reps = set(alternative_reps)
    
    return alternative_reps


def verify_logical_ops(given_set, new_set, stabilizers):
    ## go through list one by one to determine if logical operator is valid
    valid_logical_ops = True
    
    # create complete set of logical operators according to Stabilizers and verify original logical operator is in list
    
    alternative_reps = []
    
    for i in range(len(given_set)):
        alternative_reps = gen_logical_op_set(given_set[i], stabilizers)
        
        if new_set[i] in alternative_reps:
            continue
        else:
            valid_logical_ops = False
            break
            
    return valid_logical_ops


# stab_gen = ["ZIZIZXZXIZ", "XZIXIXZZZI", "XYZXYYXXZI", "XYYZYXYXIZ", "YXXXZIYIYX", "ZZZYZXZIIZ", "YYIZZIZXIY", "YXXXXXXIYI"]
# logical_ops = ["ZXXXYXIIII", "IZZYXIYIII", "ZYXIZZXYYZ", "IZYXXZZYII"]
# new_set = ["IXYXXIZXIZ", "ZZIYYXXXIZ", "IYYIIYYZYI", "ZZXXYYIZIZ"]



# stab_gen = ["XYZI", "YXIY", "IYIZ"]
# logical_ops = ["YXZY", "ZIXZ"]

# new_set = ["YZZX", "ZYXI"]

# start_time = time.time()

# for i in range(1000000):
#     temp_bool = verify_logical_ops(logical_ops, new_set, stab_gen)

# end_time = time.time()
# elapsed_time = end_time - start_time

# print(temp_bool)
        
# print(f"Elapsed time: {elapsed_time} seconds")
        
    
    

        
        