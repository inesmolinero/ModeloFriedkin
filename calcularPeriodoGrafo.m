function d = calcularPeriodoGrafo(adj)
% calcularPeriodoGrafo: calcula el periodo de un grafo dirigido fuerte
n = size(adj,1);
periods = [];
A = adj;
% Se buscan primeros 2n caminos de vuelta a cada nodo
for k = 1:2*n
    if any(diag(A))
        periods(end+1) = k; %#ok<AGROW>
    end
    A = A * adj;
end
if isempty(periods)
    d = NaN;
else
    d = periods(1);
    for val = periods(2:end)
        d = gcd(d,val);
    end
end
end
