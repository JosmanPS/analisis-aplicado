function [x, iter] = grad_conjug(n)

% =========================================================
%
% Este programa genera una matriz aleatoria con el numero
% de condicion indicada y el gradiente necesario para que
% la solucion optima de la cuadratica generada por este 
% matriz y vector tenga como solucion el vector de 1s
%
% 12 Febrero 2015
%
% Jose Manuel Proudinat Silva
% 130056
%
% Entrada:
% 	- n : indica el tamano de la matriz cuadrada n x n
% 	- ncond : numero de condicion de la matriz
%
% Output:
% 	- A  : matriz aleatoria generada con condicion 'ncond'
% 	- g  : gradiente generado para que la solucion minima de 
%          la cuadratica generada por A y g sea el vector 
%          de unos
% 	- x0 : vector aleatorio de tamano 'n' 
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

while(iter < n && norm(r0) > tol)

	% Actualizamos
	alpha = (r0' * r0) / (d' * A * d);
	x = x + alpha * d;
	r1 = r0 - alpha * A * d;
	Beta = (r1' * r1) / (r0' * r0);
	d = r1 + Beta * d;
	r0 = r1;
	iter = iter + 1;

end
