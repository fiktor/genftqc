import numpy as np

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

