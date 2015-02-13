function x = GC(n)

% =========================================================
%
% Probamos la programacion del metodo de gradiente conjugado
%
% 12 Febrero 2015
%
% Jose Manuel Proudinat Silva
% 130056
%
% Entrada:
% 	- n : indica el tamano de la matriz cuadrada n x n
%
% Output:
% 	- x  : vector solucion del problema de optimizacion
%
% =========================================================

% Calculamos los valores iniciales
[A, g] = P1_matriz(n, 10);
x = zeros(n,1);
tol = 1e-5;

n = size(A);
n = n(1);
r0 = g - A * x;
d = r0;
iter = 1;

% imprimimos
fprintf(1, '	%3i			%1.5e  \n', iter, norm(r0));

while(iter < n && norm(r0) > tol)

	% Actualizamos
	alpha = (r0' * r0) / (d' * A * d);
	x = x + alpha * d;
	r1 = r0 - alpha * A * d;
	Beta = (r1' * r1) / (r0' * r0);
	d = r1 + Beta * d;
	r0 = r1;
	iter = iter + 1;

	% imprimimos
	fprintf(1, '	%3i			%1.5e  \n', iter, norm(r0));

end
