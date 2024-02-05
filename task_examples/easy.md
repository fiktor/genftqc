# Task #1
```
Task: multiply_gates(gates=["x", "y", "z"])
Answer: gphase(pi / 2);
```

```
Task: fill_blank(
    gate_set=standard_clifford_gates,
    script="""
        pragma forall qubit qv
        pragma precondition q == qv
        qubit q[2];
        swap q[0], q[1];
        x q[1];
        swap q[0], q[1];
        x1? q[x2?];
        pragma postcondition q == qv
        """)
Answer:
line("x1?") = "x q[0];"
```
