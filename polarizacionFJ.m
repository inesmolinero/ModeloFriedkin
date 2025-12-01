function stats = polarizacionFJ(x, varargin)
% POLARIZACIONFJ  Calcula índices de polarización y características básicas
%   stats = polarizacionFJ(x)
%
%   Entradas:
%       x : vector de opiniones finales (Nx1)
%
%   Parámetros opcionales:
%       'alpha'   : parámetro de Esteban & Ray (por defecto = 1.3)
%       'bandwidth' : ancho de banda para detección de modas (por defecto = automático)
%
%   Salidas (estructura stats):
%       .P_ER   -> índice de Esteban & Ray (1994)
%       .P_W    -> índice de Wolfson (1994)
%       .var    -> varianza de opiniones
%       .range  -> rango de opiniones
%       .modes  -> número de modas detectadas
%       .isBimodal -> true si hay al menos 2 modas
%       .interpretacion -> texto explicativo sobre qué índice es más relevante

% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'alpha', 1.3);
addParameter(p, 'bandwidth', []);
parse(p, varargin{:});
alpha = p.Results.alpha;
bw = p.Results.bandwidth;

x = x(:);
x = x(~isnan(x)); % elimina NaN
n = numel(x);

% Normaliza a [0,1] por si acaso
x = (x - min(x)) / (max(x) - min(x) + eps);

% --- 1) Índice de Esteban & Ray (1994)
p_i = ones(n,1)/n; % pesos iguales
k = 1/2; % normalización típica
P_ER = k * sum( (p_i.^(1+alpha)) .* (abs(x - x')).*p_i', 'all' );

% --- 2) Índice de Wolfson (1994)
x_sorted = sort(x);
mu = mean(x_sorted);
mediana = median(x_sorted);
mu_R = mean(x_sorted(x_sorted >= mediana));
P_W = 2 * (mu_R - mu) / (mu + eps);

% --- 3) Varianza y rango
var_final = var(x);
range_final = max(x) - min(x);

% --- 4) Detección de modas (Kernel Density Estimate)
if isempty(bw)
    [f, xi] = ksdensity(x);
else
    [f, xi] = ksdensity(x, 'Bandwidth', bw);
end
modes = islocalmax(f);
numModes = sum(modes);

% --- 5) Interpretación
if numModes <= 1
    interpretacion = sprintf(['Distribución unimodal: consenso parcial o dispersión continua.\n' ...
                              'El índice de Wolfson (%.3f) refleja mejor la polarización global.\n'], P_W);
else
    interpretacion = sprintf(['Distribución multimodal (%d modas): presencia de grupos diferenciados.\n' ...
                              'El índice de Esteban & Ray (%.3f) es más representativo de la polarización grupal.\n'], ...
                              numModes, P_ER);
end

% --- 6) Salida
stats.P_ER = P_ER;
stats.P_W = P_W;
stats.var = var_final;
stats.range = range_final;
stats.modes = numModes;
stats.isBimodal = numModes > 1;
stats.interpretacion = interpretacion;
end
