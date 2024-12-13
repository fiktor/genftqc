import numpy as np
import sympy

def single_commutation(A,B):
    paulis = {"I", "X", "Y", "Z"}
    if A not in paulis:
        raise Exception("Pauli String A consists of non-Paulis")
    if B not in paulis:
        raise Exception("Pauli String B consists of non-Paulis")
    if A == 'I' or B == 'I':
        return 0 
    elif A == B:
        return 0
    else:
        return 1


def check_commutivity(A, B):
    if len(A) != len(B):
        raise Exception("Pauli Strings are unequal dimension")
                    
    commute = 0
                        
    for i in range(len(A)):
        commute ^ single_commutation(A[i],B[i])
                                                
    return commute


def multiply_pauli(A,B):
    paulis = {"I", "X", "Y", "Z"}
    if A not in paulis:
        raise Exception("Pauli String A consists of non-Paulis")
    if B not in paulis:
        raise Exception("Pauli String B consists of non-Paulis")
        
    if A == "I":
        return B
    elif B == "I":
        return A
    elif A == B:
        return "I"
    elif (A == "X" and B == "Y") or (A == "Y" and B == "X"):
        return "Z"
    elif (A == "X" and B == "Z") or (A == "Z" and B == "X"):
        return "Y"
    elif (A == "Z" and B == "Y") or (A == "Y" and B == "Z"):
        return "X"


def multiply_pauli_strings(A,B):
    if len(A) != len(B):
        raise Exception("Pauli Strings are unequal dimension")
        
    prod = ""
    for i in range(len(A)):
        prod += multiply_pauli(A[i],B[i])
        
    return prod


def multiply_multiple_pauli_strings(pauli_strings, n):
    # for i in range(len(pauli_strings)):
    #     if len(i) != len(pauli_strings[0]):
    #         raise Exception("Pauli Strings are unequal dimension")
    
    final_pauli_string = ''
    for i in range(n):
        final_pauli_string += 'I'  
        
    for pauli_string in pauli_strings:
        final_pauli_string = multiply_pauli_strings(final_pauli_string, pauli_string)
        
    return final_pauli_string

def gen_pauli_strings(n, generators):
    pauli_strings_to_be_mult = []
    for i in range(4**n):
        pauli_string_gen_group = []
        for j in range(2*n):
            ### NOTE THAT THE INTEGER IN THE BITSTRING EXPRESSION BELOW MUST BE THE NUMBER OF GENERATORS IN THE GENERATOR LIST ###
            if f'{i:0{2*n}b}'[j] == '1':
                pauli_string_gen_group.append(generators[j])
            
        pauli_strings_to_be_mult.append(pauli_string_gen_group)

#     print(len(pauli_strings_to_be_mult))
#     print(pauli_strings_to_be_mult)

    full_pauli_string_grp = []

    for pauli_strings in pauli_strings_to_be_mult:
        full_pauli_string_grp.append(multiply_multiple_pauli_strings(pauli_strings, n))
#     print()
#     print(len(full_pauli_string_grp))
#     print(full_pauli_string_grp)

#     full_pauli_string_grp1 = list(set(full_pauli_string_grp))
#     print()
#     print(len(full_pauli_string_grp1))
#     print(full_pauli_string_grp1)

#     num_matching_entries = 0
#     for i in range(len(full_pauli_string_grp)):
#         for j in range(len(full_pauli_string_grp1)):
#             if full_pauli_string_grp[i] == full_pauli_string_grp1[j]:
#                 num_matching_entries += 1
            
#     print()
#     if num_matching_entries == int(np.sqrt(len(full_pauli_string_grp) * len(full_pauli_string_grp1))):
#         print("Indentical Lists")

    
    return full_pauli_string_grp

def gen_pauli_string_generators(n):
    generators = []
    num_generators = 2*n
    
    for i in range(n):
        generator = ''
        for j in range(n):
            if i == j:
                generator += 'X'
            else:
                generator += 'I'
                
        generators.append(generator)
        
    for i in range(n):
        generator = ''
        for j in range(n):
            if i == j:
                generator += 'Z'
            else:
                generator += 'I'
                
        generators.append(generator)
        
    return generators


# Generates all n qubit Pauli Strings, except the identity string
def generate_n_qubit_paulis(n):
    generators = gen_pauli_string_generators(n)
    #print(generators)
    pauli_strings = gen_pauli_strings(n, generators)

    return pauli_strings

#print(generate_n_qubit_paulis(4))


class QECC:
    
    def __init__(self, stabilizers : list, logical_x : list, logical_z : list):
        #stabilizer generators
        self.stabilizers = stablilizers
        #logical X operator generators
        self.logical_x = logical_x
        #logical Z operator generators
        self.logical_z = logical_z
        
    
    
    
def vectorize_pauli_string(A):
    n = len(A)
    vector = np.zeros(2*n, dtype = int)
    
    for i in range(len(A)):
        if A[i] == 'X':
            vector[i] = 1
        elif A[i] == 'Y':
            vector[i] = 1
            vector[n + i] = 1
        elif A[i] == 'Z':
            vector[n + i] = 1
        else:
            continue
            
    return vector

def linear_indep_subset(stabilizers):
    vectors = []
    for stabilizer in stabilizers:
        vector = vectorize_pauli_string(stabilizer)
        vectors.append(vector)
        
    matrix = np.array(vectors)
    #print(vectors)
    #M = sympy.Matrix([[1, 0, 1], [1, 1, 0], [0, 1, 1]])
    #M = sympy.Matrix(vectors)
    #print(M)
    #print(M.nullspace(iszerofunc=lambda x: x % 2 == 0))
    #print(matrix)
    print()
    
    _, inds = sympy.Matrix(matrix).T.rref()
    #rref, inds = sympy.Matrix(np.array([[1, 1, 0], [0, 1, 1], [1, 0, 1]])).T.rref()
    #rref, inds = sympy.Matrix(np.array([[0.9, -0.1, -0.2, 0], [-0.8, 0.9, -0.4, 0], [-0.1, -0.8, 0.6, 0]])).T.rref(iszerofunc=lambda x:abs(x)<1e-9)
    
#     print(rref.T)
#     print(inds)
    
    independent_stabilizers = []
    for i in inds:
        independent_stabilizers.append(stabilizers[i])
        
    return independent_stabilizers

### TESTING ###
test_stabilizers1 = ['IIIXXXX', 'IXXIIXX', 'XIXIXIX', 'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ']

test_stabilizers2 = ['IIIXXXX', 'IXXIIXX', 'XIXIXIX', 'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ', 'XXIIXXI', 'IIIYYYY', 'ZIZXYXY']

reduced1 = linear_indep_subset(test_stabilizers1)
reduced2 = linear_indep_subset(test_stabilizers2)

print(reduced1)
print(reduced2)
### END TESTING ###


def choose_stabilizers(n, k, num_stab_gen):
    
    pauli_strings = generate_n_qubit_paulis(n)
    
    # randomly choose n-k pauli strings from the list of pauli strings
    # check if they commute
        # if one doesn't commute, get rid of it and choose another random stabilizer
    # check if they are linearly independent
        # if they aren't, find the minimal set which is linearly independent and add more which are linearly independent s.t. we have n-k of them
    # if all conditions are satisfied, return them
    
    
    return


def get_logical_ops(stabilizers):
    return


def generate_qec_code(n, k):
    
    num_stab_gen = n - k
    
    stabilizers = choose_stabilizers(n, k, num_stab_gen)
    logical_x, logical_z = get_logical_ops(stabilizers)
    
    code = QECC(stabilizers, logical_x, logical_z)
    
    return code


    
def generate_codes(num_codes, min_n, max_n):
    physical_qubit_counts = list(range(min_n, max_n +1))
    
    # Need to decide how many [[n,k,d]] codes to create for each value of n
    # What distribution is good?
    num_codes_per_n = []
    
    codes = []
    
    for i in range(len(num_codes_per_n)):
        for j in range(num_codes_per_n[i]):
            codes.append(generate_qec_code(physical_qubit_counts[i], 1))
        
    
    
    
    