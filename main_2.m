%% Simulación del Modelo Friedkin-Johnsen en Redes Erdős-Rényi
% Análisis de propagación de opiniones con agentes negativos (trolls)
% en diferentes configuraciones de red y centralidad

clear; close all; clc;

%% Configuración inicial
n = 300;                    % Número de nodos
seed = 123;                 % Semilla para reproducibilidad
max_iter = 500;             % Iteraciones máximas de simulación
tol = 1e-6;                 % Tolerancia de convergencia

% Cargar opiniones iniciales desde archivo
data_file = fullfile(fileparts(mfilename('fullpath')), 'data', 'Twitter_Data_scores.csv');
opiniones_raw = readmatrix(data_file);
opiniones_raw = opiniones_raw(~isnan(opiniones_raw(:)));

% Generar vector de opiniones iniciales
rng(seed);
if length(opiniones_raw) >= n
    x0_vec = opiniones_raw(randperm(length(opiniones_raw), n));
else
    temp = repmat(opiniones_raw, ceil(n/length(opiniones_raw)), 1);
    x0_vec = temp(randperm(length(temp), n));
end

% Normalizar opiniones al rango [0,1]
x0_vec = (x0_vec - min(x0_vec)) / (max(x0_vec) - min(x0_vec));
fprintf('Opiniones iniciales: rango [%.3f, %.3f]\n', min(x0_vec), max(x0_vec));

%% Definir escenarios de simulación
% Tres regímenes de densidad de red
p_base = log(n) / n;
regimenes = {'desconectada', 'umbral', 'fuerte'};
densidades = [0.5, 1.0, 3.0] * p_base;

% Parámetros de variación
prop_trolls = 0:0.06:0.30;      % Proporción de trolls
lambdas = 0.1:0.15:0.85;        % Nivel de prejuicio
bandas = {'bajo', 'medio', 'alto'};  % Bandas de centralidad

num_replicas = 50;              % Réplicas por configuración

%% Crear directorio de resultados
timestamp = datestr(now, 'yyyy-mm-dd_HHMMSS');
dir_resultados = fullfile(pwd, ['resultados_' timestamp]);
mkdir(dir_resultados);
diary(fullfile(dir_resultados, 'log.txt'));

%% Generar banco de redes base
fprintf('Generando %d redes por régimen...\n', num_replicas);
redes = struct();

for i = 1:length(regimenes)
    reg = regimenes{i};
    p = densidades(i);
    
    redes.(reg).A = cell(num_replicas, 1);
    redes.(reg).W = cell(num_replicas, 1);
    redes.(reg).pr = cell(num_replicas, 1);
    
    for r = 1:num_replicas
        rng(seed + 1000*i + r);
        
        % Crear red aleatoria Erdős-Rényi
        A = rand(n) < p;
        A(1:n+1:end) = 0;  % Sin auto-lazos
        
        % Asegurar que todos los nodos tienen al menos una salida
        grados = sum(A, 2);
        for nodo = find(grados == 0)'
            destino = mod(nodo, n) + 1;
            A(nodo, destino) = 1;
        end
        
        % Normalizar matriz de influencia
        W = A ./ sum(A, 2);
        
        % Calcular centralidad PageRank
        G = digraph(A);
        pr = centrality(G, 'pagerank');
        pr = pr / sum(pr);
        
        redes.(reg).A{r} = A;
        redes.(reg).W{r} = W;
        redes.(reg).pr{r} = pr;
    end
end

%% Simulación principal
resultados = {};
fila = 1;

total_configs = length(regimenes) * length(prop_trolls) * length(lambdas) * length(bandas);
fprintf('Simulando %d configuraciones × %d réplicas...\n', total_configs, num_replicas);

for ir = 1:length(regimenes)
    reg = regimenes{ir};
    
    for pt = prop_trolls
        num_trolls = round(pt * n);
        
        for lam = lambdas
            
            for ib = 1:length(bandas)
                banda = bandas{ib};
                
                % Simular réplicas
                for r = 1:num_replicas
                    A = redes.(reg).A{r};
                    W = redes.(reg).W{r};
                    pr = redes.(reg).pr{r};
                    
                    % Seleccionar trolls según banda de centralidad
                    trolls = seleccionar_trolls(pr, num_trolls, banda);
                    normales = setdiff(1:n, trolls);
                    
                    % Configurar opiniones iniciales
                    x0 = x0_vec;
                    x0(trolls) = -1;  % Trolls con opinión negativa
                    
                    % Configurar prejuicios
                    lambda_vec = zeros(n, 1);
                    lambda_vec(normales) = lam;
                    
                    % Calcular punto fijo teórico
                    I = eye(n);
                    Lambda = diag(lambda_vec);
                    x_equilibrio = (I - Lambda*W) \ ((I - Lambda)*x0);
                    
                    % Simular evolución temporal
                    X = simular_modelo(W, x0, lambda_vec, max_iter, x_equilibrio, tol);
                    x_final = X(:, end);
                    
                    % Calcular métricas
                    rango = max(x_final) - min(x_final);
                    desv = std(x_final);
                    mediana = median(x_final);
                    media = mean(x_final);
                    prop_neg = mean(x_final(normales) < 0);
                    
                    % Tiempo de convergencia
                    dist = vecnorm(X - x_equilibrio, inf, 1);
                    t_conv = find(dist < tol, 1);
                    if isempty(t_conv), t_conv = NaN; end
                    
                    % Componentes fuertemente conexas
                    G_dir = digraph(A);
                    [~, sizes] = conncomp(G_dir, 'Type', 'strong');
                    num_comp = length(sizes);
                    
                    % Guardar resultados
                    resultados{fila, 1} = reg;
                    resultados{fila, 2} = banda;
                    resultados{fila, 3} = pt;
                    resultados{fila, 4} = lam;
                    resultados{fila, 5} = r;
                    resultados{fila, 6} = rango;
                    resultados{fila, 7} = desv;
                    resultados{fila, 8} = mediana;
                    resultados{fila, 9} = media;
                    resultados{fila, 10} = prop_neg;
                    resultados{fila, 11} = t_conv;
                    resultados{fila, 12} = num_comp;
                    
                    fila = fila + 1;
                end
            end
        end
    end
end

%% Exportar resultados
headers = {'regimen', 'banda_centralidad', 'prop_trolls', 'lambda', 'replica', ...
           'rango', 'desviacion', 'mediana', 'media', 'prop_negativos', ...
           'tiempo_conv', 'num_componentes'};

tabla = cell2table(resultados, 'VariableNames', headers);
writetable(tabla, fullfile(dir_resultados, 'resultados.csv'));
save(fullfile(dir_resultados, 'datos.mat'), 'tabla');

fprintf('\nSimulación completada. Resultados en: %s\n', dir_resultados);
diary off;

%% ========== FUNCIONES AUXILIARES ==========

function trolls = seleccionar_trolls(pr, num, banda)
    % Selecciona nodos según su centralidad
    n = length(pr);
    [~, idx_ord] = sort(pr, 'descend');
    tercil = floor(n/3);
    
    switch banda
        case 'bajo'
            rango = idx_ord(2*tercil+1:end);
        case 'medio'
            rango = idx_ord(tercil+1:2*tercil);
        case 'alto'
            rango = idx_ord(1:tercil);
    end
    
    if num > length(rango)
        trolls = rango;
    else
        trolls = rango(randperm(length(rango), num));
    end
end

function X = simular_modelo(W, x0, lambda_vec, T, x_eq, tol)
    % Simula la evolución del modelo Friedkin-Johnsen
    n = length(x0);
    X = zeros(n, T+1);
    X(:,1) = x0;
    
    for t = 1:T
        x_actual = X(:,t);
        x_nuevo = (1 - lambda_vec) .* x0 + lambda_vec .* (W * x_actual);
        X(:,t+1) = x_nuevo;
        
        % Detener si converge
        if norm(x_nuevo - x_eq, inf) < tol
            X = X(:, 1:t+1);
            break;
        end
    end
end