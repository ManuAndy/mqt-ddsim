// Feynman -- quantum circuit toolkit
// grover_5.qc
//   Qubits: 7
//   H: 64
//   T: 224
//   X: 57
//   cnot: 224
//   Depth: 313
//   T depth: 96

// Permutation
// 0 1 2 3 4 5 6
// 4 1 2 3 5 0 6
OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[7];
// [0] [1] [2] [3] [4] [5]
x qubits[6];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[0];
x qubits[3];
x qubits[2];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
// [0] [1 2] [3] [4] [5]
cx qubits[3],qubits[4];
// [0] [1 2] [3 4] [5]
cx qubits[4],qubits[2];
// [0] [1 2 3 4] [5]
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
h qubits[1];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
// [0] [1 2 3 4 5]
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[1];
t qubits[1];
t qubits[0];
t qubits[6];
cx qubits[1],qubits[0];
// [0 1 2 3 4 5]
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
tdg qubits[1];
tdg qubits[0];
t qubits[6];
cx qubits[0],qubits[1];
tdg qubits[1];
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
cx qubits[1],qubits[0];
h qubits[1];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
x qubits[0];
x qubits[3];
x qubits[2];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
t qubits[4];
t qubits[5];
t qubits[0];
cx qubits[4],qubits[5];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[0];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
cx qubits[4],qubits[5];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[0];
x qubits[3];
x qubits[2];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[1];
t qubits[1];
t qubits[0];
t qubits[6];
cx qubits[1],qubits[0];
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
tdg qubits[1];
tdg qubits[0];
t qubits[6];
cx qubits[0],qubits[1];
tdg qubits[1];
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
cx qubits[1],qubits[0];
h qubits[1];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
x qubits[0];
x qubits[3];
x qubits[2];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
t qubits[4];
t qubits[5];
t qubits[0];
cx qubits[4],qubits[5];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[0];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
cx qubits[4],qubits[5];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[0];
x qubits[3];
x qubits[2];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[1];
t qubits[1];
t qubits[0];
t qubits[6];
cx qubits[1],qubits[0];
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
tdg qubits[1];
tdg qubits[0];
t qubits[6];
cx qubits[0],qubits[1];
tdg qubits[1];
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
cx qubits[1],qubits[0];
h qubits[1];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
x qubits[0];
x qubits[3];
x qubits[2];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
t qubits[4];
t qubits[5];
t qubits[0];
cx qubits[4],qubits[5];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[0];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
cx qubits[4],qubits[5];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[0];
x qubits[3];
x qubits[2];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[1];
t qubits[1];
t qubits[0];
t qubits[6];
cx qubits[1],qubits[0];
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
tdg qubits[1];
tdg qubits[0];
t qubits[6];
cx qubits[0],qubits[1];
tdg qubits[1];
cx qubits[0],qubits[6];
cx qubits[6],qubits[1];
cx qubits[1],qubits[0];
h qubits[1];
t qubits[4];
t qubits[5];
t qubits[1];
cx qubits[4],qubits[5];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[1];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[1];
cx qubits[1],qubits[4];
cx qubits[4],qubits[5];
h qubits[1];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
x qubits[0];
x qubits[3];
x qubits[2];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
t qubits[4];
t qubits[5];
t qubits[0];
cx qubits[4],qubits[5];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
tdg qubits[4];
tdg qubits[5];
t qubits[0];
cx qubits[5],qubits[4];
tdg qubits[4];
cx qubits[5],qubits[0];
cx qubits[0],qubits[4];
cx qubits[4],qubits[5];
h qubits[4];
t qubits[2];
t qubits[3];
t qubits[4];
cx qubits[2],qubits[3];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
tdg qubits[2];
tdg qubits[3];
t qubits[4];
cx qubits[3],qubits[2];
tdg qubits[2];
cx qubits[3],qubits[4];
cx qubits[4],qubits[2];
cx qubits[2],qubits[3];
h qubits[4];
x qubits[2];
x qubits[3];
x qubits[5];
x qubits[0];
h qubits[2];
h qubits[3];
h qubits[5];
h qubits[0];
