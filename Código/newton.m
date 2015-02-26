function [x] = newton ( x0, funder, maxiter, tol );
%
% metodo de Newton global (busqueda lineal con condiciones fuertes
% de Wolfe) para optimizacion de funciones sin restricciones
% j-l morales 19 de febrero, 2015, ITAM
%
%--------------------------------------------------------------------------
%
% ... tolerancias y valores iniciales
%
n = length(x0);
x = x0;  
k = 0; 
%
maxgc = 2*n;
tolgc = tol;
c1 = 1.0e-4; c2 = 0.9;

[ f, g, H ] = feval( funder, x );

norm_g = norm(g);

fprintf(1,'   k        f            ||g||     alfa     nfg     GC           ||r|| \n');
fprintf(1,'-----------------------------------------------------------------------\n');

while  norm_g > tol  &&  k < maxiter
    
    [ p_N, termina, normr ] = grad_conj( H, g, tolgc, maxgc, 1 );
    [ alfa, x, f, g, falla, numfg ] = ...
                           biseccion ( x, f, g, p_N, funder, c1, c2 );

    [ ff, gg, H ] = feval( funder, x );  % obtener la Hessiana en x
    
    norm_g   = norm(g); 
    k = k + 1;
    fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i  %8s  %8.2e \n', ...
                            k, f,  norm_g, alfa, numfg, termina, normr );
end
%
%==========================================================================
%
