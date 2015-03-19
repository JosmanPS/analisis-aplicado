function [ x ] = RegionConfianza( hatDelta, delta, etha, x, fun, maxiter, tol )

% Errores comunes
if delta <= 0 || delta >= hatDelta
    disp('Error: tu delta0 no está bien definida');
    return;
end

if etha < 0 || etha >= 0.25
    disp('Error: tu etha no está bien definida');
    return;
end

% Definimos el modelo local
m = @(p, g, B) (0.5*p'*B*p + g'*p);

% Valores iniciales
k = 0;
[~, g, B] = fun(x);
n = length(g);

while norm(g) > tol && k < maxiter
    
    % Encontramos p_k del problema con restricciones
    p = grad_conj(B, g, 1e-5, 1e-5, 2 * n);
    
    % Evaluamos el ratio
    ratio = ( fun(x) - fun(x + p) ) / ( m(zeros(n,1), g, B) - m(p, g, B) );
    
    if ratio < 0.25
        delta = 0.25 * norm(p);
    else
        if ratio > 0.75 && norm(p) == delta
            delta = min(2 * delta, hatDelta);
        end
    end
    
    if ratio > etha
        x = x + p;
    end
    
    [~, g, B] = fun(x);
    k = k + 1;
    
end

