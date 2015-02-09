function resultados_max_desc

n = 5;
ncond = 1.0d1;

[H, g, x] = matriz(n, ncond);

max_iter = 1000;
tol = 1e-8;

max_mad(H, g, x, max_iter, tol);



