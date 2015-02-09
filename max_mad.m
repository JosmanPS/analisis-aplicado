function [x, norma, f, iter, info] = max_mad(H, g, x, max_iter, tol)

% Calculamos los valores iniciales
iter = 0;
grad = H*x + g;
norma = norm(grad);
f = g'*x + (x'*H*x)/2;
info = 'Ha comenzado el programa';
alpha = (norma^2)/(grad'*H*grad);

fprintf(1,'    iter    f          norma    \n\n');

% Comenzamos el descenso
while(norma > tol && iter < max_iter)
    x = x - alpha.*grad;
    
    % Actualizamos
    grad = H*x + g;
    norma = norm(grad);
    f = g'*x + (x'*H*x)/2;
    alpha = (norma^2)/(grad'*H*grad);
    iter = iter + 1;
    
    fprintf(1, '    %3i    %1.4f    %1.4f    \n', iter, f, norma);
    
    % Informacion
    if(norma < tol)
        info = 'La norma del gradiente fue menor a la tolerancia en:';
    end

    if(iter >= max_iter)
        info = 'Se alcanzo el maximo de iteraciones en:';
    end
end


fprintf(1, info)
fprintf(1, '\n')
x
