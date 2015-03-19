function x = grad_conj( H, g, tol1, tol2, maxit)

% =========================================================
%
% Este programa encuentra la solucion (x) del problema:
%				Hx = g
% a traves del metodo de gradiente conjugado. Notar que este
% problema es equivalente al de encontrar el minimos sobre
% una funcion cuadratica
%
% 8 Marzo 2015
%
% Jose Manuel Proudinat Silva
% 130056
%
% Entrada:
% 	- H : matriz del problema (Hessiana)
% 	- g : vector de respuestas (gradiente)
% 	- tol1 : tolerancia para la norma del residuo
% 	- tol2 : tolerancia para la norma cuadrada de 'd' 
%			 respecto a H
% 	- maxit : maximo numero de iteraciones
%
% Salida:
% 	- x : el vector solucion al problema 
%
% =========================================================

% Calculamos los valores iniciales
n = length(g);
x = zeros(n,1);
r = g - H * x;
d = r;
iter = 0;

% Guardamos los resultados para reducir operaciones
rTr = r' * r;
Hd = H * d;
dHd = d' * Hd;
norm_r = sqrt(rTr);

while(iter < maxit && norm_r > tol1 && dHd > tol2)

	% Actualizamos los valores
	alpha = rTr / dHd;
	x = x + alpha * d;
	r = r - alpha * Hd;
	Beta = (r' * r) / rTr;
	d = r + Beta * d;

	% Guardamos los resultados para reducir operaciones
	rTr = r' * r;
	Hd = H * d;
	dHd = d' * Hd;
	norm_r = sqrt(rTr);

	% Siguiente iteracion
	iter = iter + 1;

end
