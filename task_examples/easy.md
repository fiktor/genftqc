# Task #1
Goals:
* Introduce the notion of Q&A format.
* Introduce the standard gates:
  - Standard gates are those listed in
    [stdgates.inc](https://github.com/openqasm/openqasm/blob/main/examples/stdgates.inc)
    excluding backwards compatibility gates
    `+` popular gates implemented in hardware (iswap?) + `id` for identity (`gphase(0)`).
  - TODO: list extra gates.
* Introduce `ignore_global_phase=True` option.

`multiply_gates` task description:
* A list of gates is presented as `gates` argument;
* Number of qubits is implicit: 1..3;
* A product of these gates is guaranteed to be a single named gate
  (up to global phase if `ignore_global_phase=True`);
* That gate should be returned as an answer.

Sometimes multiple descriptions are valid. E.g. `crx(pi) q[0] q[1];` and `cx q[0] q[1];`.
In this case the shortest description gets reward $1.0$ and $0.1$ is subtracted
for each additional parameter or qubit used (in the example above 1.0 is given for the first
solution and 0.9 for the second). Wrong or syntactically incorrect solution gets 0 points.

## Example 1
```
Task: multiply_gates(gates=["x", "y", "z"])
Answer: gphase(pi / 2);
```

Answer `gphase(pi / 2) q[0];` is also accepted.

## Example 2
```
Task: multiply_gates(
  gates=["x q[0]", "swap q[0] q[1]", "cx q[0] q[1]", "x q[1]"],
  ignore_global_phase=True)
Answer: cx q[1] q[0];
```

# Task #5
```
Task: fill_blank(
    gate_set=standard_clifford_gates,
    script="""
        pragma forall qstate[2] qv;
        pragma precondition q == qv;
        qubit q[2];
        swap q[0], q[1];
        x q[1];
        swap q[0], q[1];
        x1? q[x2?];
        pragma postcondition q == qv;
        """)
Answer:
line("x1?") = "x q[0];"
```

# Task #6
Goals:
* Introduce Stabilizer desription of QEC code
* Connect stabilizer description of QEC code to code words
* Here, any circuit that generates a codeword is a valid answer

```
Task: gen_L_state(
    stab_gen = ['XXXX','ZZZZ'] 
    )
Answer:
qreg qubits[4];
U(pi/2, 0, pi) qubits[0];
cx qubits[0] qubits[1];
cx qubits[1] qubits[2];
cx qubits[2] qubits[3];
```

# Task #7
Goals:
* Generate encoding circuits given a description of a QEC code through stabilizers and logical operators
* Here, only circuits that generate the logical |0> state are valid answers
```
Task: gen_enc_circ_0(
    stab_gen = ['XXXX','ZZZZ'],
    qL1 = ['XIXI', 'ZZII'],
    qL2 = ['XXII', 'ZIZI']
    )
Answer:
qreg qubits[4];
U(pi/2, 0, pi) qubits[0];
cx qubits[0] qubits[1];
cx qubits[1] qubits[2];
cx qubits[2] qubits[3];
```

# Task #8
Goals:
* Generate a general encoding circuit that takes an arbitrary K qubit physical state on Q0-QK (where K is the number of logical qubits), and generates the logical equivalent of that state given a set of stabilizers and logical operators
* This task aims to teach the agent about stabilizer algebra

```
Task: gen_enc_circ(
    stab_gen = ['XXXX','ZZZZ'],
    qL1 = ['XIXI', 'ZZII'],
    qL2 = ['XXII', 'ZIZI']
    )
Answer:
qreg qubits[4];
U(pi/2, 0, pi) qubits[2];
cx qubits[1] qubits[2];
cx qubits[1] qubits[3];
cx qubits[2] qubits[3];
cx qubits[2] qubits[0];
cx qubits[0] qubits[1];
```

# Task #9
Goals:
* Introducing ideas of Fault Tolerance (i.e. the generated circuits should be able to sustain at least one single weight t error without the error spreading to an uncorrectable error)
* We will only look at X type errors in this task

```
Task: gen_ft_X_enc_circ_0(
    stab_gen = ['XXXX','ZZZZ'],
    qL1 = ['XIXI', 'ZZII'],
    qL2 = ['XXII', 'ZIZI']
)
Answer:
qreg qubits[5];
U(pi/2, 0, pi) qubits[1];
cx qubits[1] qubits[2];
cx qubits[2] qubits[3];
cx qubits[1] qubits[0];
cx qubits[3] qubits[4];
cx qubits[0] qubits[4];
flag_bit = measure qubits[4];
```

# Task #10
Goals:
* Introducing ideas of Fault Tolerance (i.e. the generated circuits should be able to sustain at least one single weight t error without the error spreading to an uncorrectable error)
* We will only look at Z type errors in this task

```
Task: gen_ft_Z_enc_circ_0(
    stab_gen = ['XXXX','ZZZZ'],
    qL1 = ['XIXI', 'ZZII'],
    qL2 = ['XXII', 'ZIZI']
)
Answer:
qreg qubits[5];
U(pi/2, 0, pi) qubits[1];
cx qubits[1] qubits[2];
cx qubits[2] qubits[3];
cx qubits[1] qubits[0];
cx qubits[3] qubits[4];
cx qubits[0] qubits[4];
flag_bit = measure qubits[4];
```

# Task #11
Goals:
* Introducing ideas of Fault Tolerance (i.e. the generated circuits should be able to sustain at least one single weight t error without the error spreading to an uncorrectable error)
* We will look at both X and Z type errors in this task

```
Task: gen_ft_enc_circ_0(
    stab_gen = ['XXXX','ZZZZ'],
    qL1 = ['XIXI', 'ZZII'],
    qL2 = ['XXII', 'ZIZI']
)
Answer:
qreg qubits[5];
U(pi/2, 0, pi) qubits[1];
cx qubits[1] qubits[2];
cx qubits[2] qubits[3];
cx qubits[1] qubits[0];
cx qubits[3] qubits[4];
cx qubits[0] qubits[4];
flag_bit = measure qubits[4];
```