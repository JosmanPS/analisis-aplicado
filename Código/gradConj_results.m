function gradConj_results

[A, g, x] = P1_matriz(10, 10);

tol = 1e-5;

[x,iter] = grad_conj(A, g, x, tol);

x
iter