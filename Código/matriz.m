function [A, g, x0] = matriz(n, ncond)

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
g = randn(n,1);

% Generamos el punto inicial
x0 = randn(n,1);

end