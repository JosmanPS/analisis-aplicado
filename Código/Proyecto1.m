
% Generamos una funcion en R^10 con condicion [1e1, 1e2, 1e3]
n = 10;
nconds = [1e1 1e2 1e3];
[A, g, x0] = matriz(10, nconds(1));

% Predecir numero de pasos para alcanzar la condicion:
% ||x_{k+1} - x^*|| <= 10e-5