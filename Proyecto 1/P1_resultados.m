function P1_resultados

% =========================================================
%
% (Descripcion)
%
% 12 Febrero 2015
%
% Jose Manuel Proudinat Silva
% 130056
%
% =========================================================


% Guardamos los numeros de condicion para los que haremos el
% experimento
conds = [1.0d1 1.0d2 1.0d3];
prints = [5 15 100];

% Por la forma en que construimos A y g esta es la solcuion
x_opt = ones(10,1);
tol = 1.0d-5;

for i = 1:3

	ncond = conds(i);
	nprint = prints(i);

	fprintf(1, '============================================================== \n');
	fprintf(1, 'Resultados para el experimento con condicion: ');
	fprintf(1, '%1.1e \n', ncond);

	% Obtenemos las matrices para el experimento
	[A, g, x] = P1_matriz(10, ncond);

	% Calculamos ambos lados de la desigualdad (1) del proyecto
	% Primero el lado derecho ya que no necesitamos x_{k+1}
	eigValues = eig(A);
	eigMin = eigValues(1);
	eigMax = eigValues(length(eigValues));

	eigValues = ((eigMax - eigMin) / (eigMax + eigMin)) ^ 2;
	l_der = eigValues * ((x - x_opt)' * A * (x - x_opt));

	% Para el lado izquierdo calculamos x_{k+1}
	grad_f = A * x + g; 
	norma = norm(grad_f);
	alpha = norm(grad_f, 2)^2 / (grad_f' * A * grad_f); 
	x = x - alpha * grad_f;

	% Calculamos el lado izquierdo
	l_izq = ((x - x_opt)' * A * (x - x_opt));

	% Predecimos el numero de pasos en que se llegara a que:
	% ||x_{k+1} - x_opt|| <= tol
	predStep = ceil((log(tol) - log(l_der / eigValues)) / log(eigValues));
	fprintf(1, 'Predecimos que llegaremos al resultado en: ');
	fprintf(1, '%4i pasos \n\n', predStep);

	% Imprimimos el primer calculo
	desig = l_izq <= l_der;
	igualdad = norm(l_izq - l_der) < 1e-3;
	iter = 1;
	fprintf(1,'    iter    lado izq      lado der     desigualdad    igualdad    \n\n');
	fprintf(1, '    %4i     %1.3e    %1.3e         %1i             %1i   \n', iter, l_izq, l_der, desig, igualdad);


	% Hacemos el proceso de optimizacion y observamos que pasa
	% con la desigualdad
	while(norma > tol && iter < 10000)
    
	    % Actualizamos
	    l_der = eigValues * ((x - x_opt)' * A * (x - x_opt));
	    grad_f = A * x + g; 
	    norma = norm(grad_f);
	    alpha = norma^2 / (grad_f' * A * grad_f); 
	    x = x - alpha * grad_f;
	    l_izq = ((x - x_opt)' * A * (x - x_opt));
	    desig = l_izq <= l_der;
		igualdad = norm(l_izq - l_der) < 1e-3;
	    iter = iter + 1;

	    if mod(iter, nprint) == 0
	    	% Imprimimos
	    	fprintf(1, '    %4i     %1.3e    %1.3e         %1i             %1i   \n', iter, l_izq, l_der, desig, igualdad);
	    end
	
	end

	% Esto obliga a que se imprima la ultima iteracion
	if mod(iter, nprint) > 0
    	% Imprimimos
    	fprintf(1, '    %4i     %1.3e    %1.3e         %1i             %1i   \n', iter, l_izq, l_der, desig, igualdad);
    end

    fprintf(1, '\n');
    fprintf(1, 'El valor del gradiente de la funcion fue %1.3e \ncon %4i iteraciones \n', norma, iter);
	fprintf(1, '============================================================== \n\n');

end
	