function [A, g, x0] = P1_matriz(n, ncond)

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


% Tomamos la raiz de la condicion porque al crear B = A * A'
% la condicion de B es el cuadrado de A
ncond = sqrt(ncond); 

D = zeros(n);
randn('state', 1)

% Generamos una matriz aleatoria de tamano n y hacemos descomposcion en
% valores en singulares
R = randn(n);
[U, ~, V] = svd(R);        

% El primer valor de nuestra matriz diagonal es 1
D(1,1) = 1.0d0;

for i=2:n
    % Generamos los valores de la diagonal segun el numero de condicion
    D(i,i) = ncond^(-1/(n-i+1));
end
    
% Creamos la nueva A sdp
A = U*D*V';
A = A'*A;

% Generamos el vector g
g = -A * ones(n,1);

% Generamos el punto inicial
x0 = randn(n,1);

end