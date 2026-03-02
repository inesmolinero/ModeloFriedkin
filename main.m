%% Análisis de sentimiento y simulación Friedkin-Johnsen
clear; close all; clc;
tic;
graficar_redes = false;
%% Configuración inicial
rango_opiniones = [0, 1];
semilla_base = 123;          
max_iter = 5000;
tol_pf = 1e-6;  
n = 300; 

%% Opiniones iniciales
directorio = fileparts(mfilename('fullpath'));
cd(directorio);

% Ejecutar análisis de sentimiento
muestra_py = pyrunfile("analisis_sentimiento.py", "muestra");
% Extraer columnas
textos_py = muestra_py.get('clean_text').values.tolist();
sentimientos_py = muestra_py.get('sentimiento').values.tolist();

% Convertir a tipos de MATLAB
muestra_matlab = table(cell(textos_py)', ...
                       double(sentimientos_py)', ...
                       'VariableNames', {'texto', 'sentimiento'});
x0_raw = muestra_matlab.sentimiento;
assert(~isempty(x0_raw), 'El archivo csv está vacío');

%% Construir x0_base a tamaño n y escalar al rango de opiniones
rng(semilla_base);   

m = numel(x0_raw); % numero de valores en vector de opiniones x0

% Si el número de valores en x0_raw es mayor que n, se toman n de esos m. 
% Si no lo es, usa repmat para superar n y de ese vector extendido se toman n valores.
if m >= n
    x_vector = x0_raw(randperm(m, n)); % n índices aleatorios entre 1 y m
else
    x_vector = repmat(x0_raw, ceil(n/m), 1); %repetir el vector x0_raw (entero>= n/m) veces x 1
    x_vector = x_vector(randperm(numel(x_vector), n));
end

% Escalado al rango rango_opiniones
x_min = min(x_vector);
x_max = max(x_vector);

if x_max > x_min
    x_vector = (x_vector - x_min) / (x_max - x_min); % [0,1]
else
    x_vector = zeros(n,1);
end
x0_base = rango_opiniones(1) + x_vector*diff(rango_opiniones);
fprintf('x0_base generado en rango [%.3f, %.3f]\n', min(x0_base), max(x0_base));

%% Definir escenarios de simulación 

% Tres regimenes ER 
p_base = log(n)/n;
regimen = {'desconectada','umbral','fuerte'};
prob_regimen = [0.5, 1.0, 3.0] * p_base;

% Bandas de centralidad PageRank
bandas_pr = {'bajo_pr','medio_pr','alto_pr'};

% Proporción de trolls
prop_trolls = 0:0.06:0.30;   

% Valores lambda 
lambdas = 0.1:0.15:0.85;        
               
%% Crear directorio de resultados
str_fecha = char(datetime('now','Format','yyyy-MM-dd_HHmmss'));
carpeta_resultados = fullfile(directorio, ['resultados_' str_fecha]);

if ~exist(carpeta_resultados,'dir'), mkdir(carpeta_resultados); end
%diary(fullfile(carpeta_resultados,'log.txt'));  
fprintf('Directorio de resultados: %s\n', carpeta_resultados);


%% Establecer las 50 redes por régimen 
num_replicas = 50; % cantidad de réplicas por cada régimen de densidad
alpha = 0.9; % parámetro de salto aleatorio

fprintf('Generando banco de %d redes por régimen \n', num_replicas);
rng(semilla_base);

% fija para las redes base
redes = struct();   % Banco de réplicas de A, W y pr (vectores pageRank)
for i_reg = 1:numel(regimen)
    tipo_regimen  = regimen{i_reg};
    p_reg = prob_regimen(i_reg);

    redes.(tipo_regimen).replicas = repmat( ...
        struct('A', [], 'W', [], 'pr', []), ...
        num_replicas, 1);

    for r = 1:num_replicas
        rng(semilla_base + i_reg + r); % semilla distinta por régimen y réplica

        % ER dirigido sin autolazos
        A0 = rand(n) < p_reg;        % ER dirigido G(n,p) con probabilidad de arista p_reg
        A0(1:n+1:end) = 0;          % Anular la diagonal para evitar autolazos

        % Garantizar al menos 1 salida por nodo (fila)
        grado = sum(A0,2);         % Suma de grados de cada fila
        nodos_grado_nulo = (grado == 0);       % Filas/nodos con grado = 0
        if any(nodos_grado_nulo)
            for i = find(nodos_grado_nulo).' % Posiciones de los nodos con grado 0
                j = randi(n);       % Elige un destino aleatorio, si es el mismo i se cambia
                if j == i, j = mod(j, n) + 1; end
                A0(i,j) = 1;
            end
        end

        % Normalización por filas
        grado = sum(A0,2);
        W0 = A0 ./ grado;

        % Centralidad PageRank 
        G  = digraph(A0); % Grafo dirigido a partir de la matriz de adyacencia
        pr = centrality(G, 'pagerank', 'FollowProbability', alpha); % Centralidad de cada nodo

        % Establecer matrices de influencia y PageRank
        redes.(tipo_regimen).replicas(r).A  = A0;
        redes.(tipo_regimen).replicas(r).W  = W0;
        redes.(tipo_regimen).replicas(r).pr = pr;

    end
end

%% Establecer los parámetros para cada escenario
n_esc = numel(regimen) * numel(prop_trolls) * numel(lambdas) * numel(bandas_pr);

escenarios(n_esc) = struct('n', [], 'trolls', [], 'frac_trolls', [], ...
                               'loc', '', 'regimen', '', 'p', [], ...
                               'lam', [], 'seed', [], 'tag', '');
c = 1;
for i_reg = 1:numel(regimen)
    tipo_regimen = regimen{i_reg};
    prob_conexion = prob_regimen(i_reg);
    for j_prop_trolls = prop_trolls
        num_trolls = round(j_prop_trolls * n);
        for k_lam = lambdas
            for l_banda = 1:numel(bandas_pr)
                loc_trolls = bandas_pr{l_banda};
                etiqueta = sprintf('n%d_trolls%d_%s_%s_lam%.2f', n, num_trolls, tipo_regimen, loc_trolls, k_lam);
                escenarios(c) = struct( ...
                    'n',    n, ...
                    'trolls',  num_trolls, ...
                    'frac_trolls', j_prop_trolls, ...
                    'loc', loc_trolls, ...
                    'regimen', tipo_regimen, ...
                    'p',    prob_conexion, ...
                    'lam',  k_lam, ...
                    'seed', semilla_base + c, ...
                    'tag',  etiqueta ...
                );
                c = c + 1;
            end
        end
    end
end
%N = numel(escenarios);
fprintf('Generados %d escenarios\n', n_esc);

%% Bucle de todos los escenarios y réplicas
total_simulaciones = n_esc * num_replicas;

resultados = table('Size', [total_simulaciones, 16], ...
    'VariableTypes', {'string','double','string','string','double','double',...
                      'double','double','double','double','double','double',...
                      'double','double','double','double'}, ...
    'VariableNames', {'tag','replica','regimen','loc','fracTrolls','lambda',...
                      'rango_final','std_final','median_final','mean_final',...
                      'prop_neg','median_norm','resid','t_conv','rho_lambdaW','nSCC'});


idx_fila = 1;

for i = 1:n_esc
    s = escenarios(i);

    for r = 1:num_replicas
        
        % Red de cada réplica (indice r)        
        A  = redes.(s.regimen).replicas(r).A;
        W  = redes.(s.regimen).replicas(r).W;
        pr = redes.(s.regimen).replicas(r).pr;

        % Seleccionar nodos para trolls
        if num_trolls > 0
            idx_banda = seleccionar_trolls(pr, s.trolls, s.loc, semilla_base + i);
            trolls = idx_banda(1:min(s.trolls, length(idx_banda)));
        else
            trolls = [];
        end
        normales = setdiff((1:s.n).', trolls(:));

        % Establecer x0 y lambdas
        x0 = x0_base; 
        x0(trolls) = -1;
        lambdas = zeros(s.n,1);
        lambdas(normales) = s.lam;
        
        % Cálculo de las componentes fuertemente conexas 
        Gd = digraph(A);
        [componentes, sizes] = conncomp(Gd, 'Type', 'strong');
        nSCC = numel(sizes);    

        % Punto fijo teórico 
        I = speye(s.n);
        Lam = spdiags(lambdas,0,s.n,s.n);
        M = (I - Lam*W);
        b = (I - Lam)*x0;
        x_obj = M \ b;
    
        % Radio espectral
        rho_LW = abs(eigs(Lam*W, 1, 'lm'));
        rho_LW = min(rho_LW, 1 - 1e-12);
        
        % Simulación
        [X, converge] = simular_friedkin(W, x0, lambdas, max_iter, x_obj, tol_pf); 
        x_final = X(:,end);
        x_final_normales = x_final(normales);

        % Tiempo de convergencia
        if converge
            t_conv = size(X,2) - 1;
        else
            t_conv = NaN; 
        end
        % Métricas
        rango_f  = max(x_final) - min(x_final);
        sd_f  = std(x_final);
        mediana_f = median(x_final);
        mean_f  = mean(x_final);

        prop_neg_f = mean(x_final_normales < 0);
        mediana_norm_f = median(x_final_normales);

        resid_f = norm(x_final - x_obj, Inf);

        % Guardar resultados de fila
        resultados(idx_fila, :) = {s.tag, r, s.regimen, s.loc, ...
                                      s.frac_trolls, s.lam, rango_f, sd_f, ...
                                      mediana_f, mean_f, prop_neg_f, mediana_norm_f, ...
                                      resid_f, t_conv, rho_LW, nSCC};
        
        idx_fila = idx_fila + 1;
      
    end
end

%% Comparación de topologías con opiniones iniciales
if graficar_redes 
    % ER ya generada
    A_er_d = redes.desconectada.replicas(1).A;
    A_er_u = redes.umbral.replicas(1).A;
    A_er_f = redes.fuerte.replicas(1).A;
    
    % Generar Small-World
    A_sw = generar_small_world(n, 4, 0.1);
    
    % Generar Scale-Free
    A_sf = generar_scale_free(n, 5, 3);
    
    figure;
    
    % ER
    subplot(1,5,1);
    G1 = digraph(A_er_d);
    h1 = plot(G1,'Layout','force');
    h1.NodeCData = x0;
    colormap(jet);
    clim([-1 1]);
    title('Erdős–Rényi- Desconectada');

    colorbar;
    subplot(1,5,2);
    G1 = digraph(A_er_u);
    h1 = plot(G1,'Layout','force');
    h1.NodeCData = x0;
    colormap(jet);
    clim([-1 1]);
    title('Erdős–Rényi -Umbral');
    colorbar;

    subplot(1,5,3);
    G1 = digraph(A_er_f);
    h1 = plot(G1,'Layout','force');
    h1.NodeCData = x0;
    colormap(jet);
    clim([-1 1]);
    title('Erdős–Rényi- Fuerte');
    colorbar;

    % Small-world
    subplot(1,5,4);
    G2 = digraph(A_sw);
    h2 = plot(G2,'Layout','force');
    h2.NodeCData = x0;
    clim([-1 1]);
    title('Small-World');
    
    % Scale-free
    subplot(1,5,5);
    G3 = digraph(A_sf);
    h3 = plot(G3,'Layout','force');
    h3.NodeCData = x0;
    clim([-1 1]);
    title('Scale-Free');
    
    sgtitle('Comparación de topologías con opiniones iniciales');
end


%% Exportar resultados
writetable(resultados, fullfile(carpeta_resultados,'resultados.csv'));
save(fullfile(carpeta_resultados,'resultados.mat'),'resultados');

fprintf('Guardado: resultados.csv y resultados.mat\n');

tiempo_total = toc;
fprintf('Tiempo total de simulación: %.2f segundos (%.2f minutos)\n', ...
        tiempo_total, tiempo_total/60);

%% ========================================================================
% FUNCIONES AUXILIARES

%%  Seleccionar trolls según banda de centralidad
function idx_banda = seleccionar_trolls(pr, n_trolls, nivel_banda, seed)
    rng(seed, 'twister');
    n = length(pr);
    [~, idx] = sort(pr, 'descend');
    
    tercil = floor(n/3);
    
    switch nivel_banda
        case 'bajo_pr'
            banda = idx(2*tercil+1:end);
        case 'medio_pr'
            banda = idx(tercil+1:2*tercil);
        case 'alto_pr'
            banda = idx(1:tercil);
        otherwise
            error('Banda desconocida: %s', nivel_banda);
    end
    
    if n_trolls > length(banda)
        warning('numTrolls (%d) > tamaño banda (%d). Se toman todos los de la banda.', n_trolls, length(banda));
        idx_banda = banda;
    else
        perm = randperm(length(banda), n_trolls);
        idx_banda = banda(perm);
    end
end

%% Función de simulación del modelo de FJ

function [X,converge] = simular_friedkin(W, x0, lambdas, n_iter, x_star, tol)
    n = length(x0);
    X = zeros(n, n_iter+1);
    X(:,1) = x0;
    converge = false;

    for t = 1:n_iter
        x_prev = X(:,t);
        x_actual = (1 - lambdas) .* x0 + lambdas .* (W * x_prev);
        X(:,t+1) = x_actual;
        
        % Convergencia anticipada
        if norm(x_actual - x_star, Inf) < tol
            X = X(:, 1:t+1);
            converge = true;
            break;
        end
    end
end
%% Generar red scale-free tipo Barabási–Albert

function A = generar_scale_free(n, m0, m)

    % m0 nodos iniciales totalmente conectados
    A = zeros(n);
    
    % Red inicial completa
    A(1:m0,1:m0) = 1;
    A(1:m0+1:end) = 0;
    
    grados = sum(A,2);
    
    for i = m0+1:n
        
        prob = grados(1:i-1) / sum(grados(1:i-1));
        conexiones = randsample(i-1, m, true, prob);
        
        for j = conexiones'
            A(i,j) = 1;
        end
        
        grados = sum(A,2);
    end
end
%% Generar red small-world tipo Watts–Strogatz

function A = generar_small_world(n, k, beta)

    % Red regular en anillo
    A = zeros(n);
    
    for i = 1:n
        for j = 1:k
            vecino = mod(i+j-1,n) + 1;
            A(i,vecino) = 1;
        end
    end
    
    % Rewire con probabilidad beta
    for i = 1:n
        for j = find(A(i,:))
            if rand < beta
                nuevo = randi(n);
                A(i,j) = 0;
                A(i,nuevo) = 1;
            end
        end
    end
end
