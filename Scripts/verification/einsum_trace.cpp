#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorSlice.h"
#include <iostream>
using namespace SCESS;
enum { I, P, Q };
int main() {
    Eigen::Tensor<double, 4> t(2, 2, 3, 3);
    for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int p = 0; p < 3; p++) for (int q = 0; q < 3; q++)
        t(a, b, p, q) = (a == b ? 1.0 : 100.0);           // diagonal (i,i) = 1, off-diagonal (i,i') = 100
    Eigen::Tensor<double, 2> r = Einsum<Indices<I, I, P, Q>, Indices<P, Q>>(t);   // as used for "iipq->pq"
    std::cout << "Einsum iipq->pq gives " << r(0, 0) << "  (trace would be 2, double sum 202)\n";
}
