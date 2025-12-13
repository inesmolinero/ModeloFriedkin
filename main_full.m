%% Simulación del modelo de Friedkin-Johnsen

clear; close all; clc;

%% Configuración inicial
rango_opiniones = [0, 1];
semilla_base = 123;          
max_iteraciones = 500;
tol_consenso = 0.05;
tol_punto_fijo = 1e-6;  % tolerancia al punto fijo
n = 300; % número de agentes

%% Opiniones iniciales
% Lectura
Dir = fileparts(mfilename('fullpath'));
x0_raw_all  = readmatrix(fullfile(Dir, 'data', 'Twitter_Data_scores.csv'));
x0_raw_all  = x0_raw_all(:);        % convierte a vector columna
x0_raw_all  = x0_raw_all(~isnan(x0_raw_all));   % elimina nulos
assert(~isempty(x0_raw_all), 'Twitter_Data_scores.csv no tiene datos numéricos válidos.');

%% Construir x0_base a tamaño n y escalar al rango de opiniones
rng(semilla_base);   

m = numel(x0_raw_all); % numero de valores en vector de opiniones x0

% Si el número de valores es mayor, se toman n de esos m. En otro caso,se
% usa repmat para superar n y de ese vector extendido se toman n valores.
if m >= n
    xv = x0_raw_all(randperm(m, n));
else
    xv = repmat(x0_raw_all, ceil(n/m), 1);
    xv = xv(randperm(numel(xv), n));
end

% Escalado al rango param.rango_opiniones
xmin = min(xv); xmax = max(xv);
if xmax > xmin
    xv = (xv - xmin) / (xmax - xmin);  % [0,1]
else
    xv = zeros(param.n,1);
end
x0_base = rango_opiniones(1) + xv*diff(rango_opiniones);
fprintf('x0_base generado: rango [%.3f, %.3f]\n', min(x0_base), max(x0_base));

%% Definir escenarios de simulación 

% Tres regímenes ER 
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
if ~isempty(mfilename)
    directorio = fileparts(mfilename('fullpath'));
else
    directorio = pwd;  
end
str_fecha       = char(datetime('now','Format','yyyy-MM-dd_HHmmss'));
carpeta_resultados  = fullfile(directorio, ['resultados_' str_fecha]);

if ~exist(carpeta_resultados,'dir'), mkdir(carpeta_resultados); end
diary(fullfile(carpeta_resultados,'log.txt'));  
fprintf('Directorio de resultados: %s\n', carpeta_resultados);





%% Banco de R redes por régimen para reproducibilidad
num_replicas = 50;          % R réplicas por cada régimen de densidad
damp = 0.9;      % damping factor: parámetro de salto aleatorio

fprintf('Generando banco de %d redes por régimen...\n', num_replicas);
rng(semilla_base);

% fija para las redes base
redes = struct();   % Banco de R réplicas de A, W y pr (vectores pageRank)
for i_reg = 1:numel(regimen)
    reg  = regimen{i_reg};
    pReg = prob_regimen(i_reg);

    redes.(reg).A  = cell(num_replicas,1);
    redes.(reg).W  = cell(num_replicas,1);
    redes.(reg).pr = cell(num_replicas,1);

    for r = 1:num_replicas
        rng(param.seed_base + 1000*i_reg + r);  % semilla distinta por regimen y réplica

        % --- ER dirigido sin autolazos ---
        A0 = rand(param.n) < pReg;        % ER dirigido G(n,p) con probabilidad de arista pReg
        A0(1:n+1:end) = 0;          % Anular la diagonal para evitar autolazos

        % --- Garantiza al menos 1 salida por nodo ---
        grado = sum(A0,2);         % Suma de grados de cada fila
        nodos_grado_nulo = (grado == 0);       % Marca nodos con grado = 0
        if any(nodos_grado_nulo)
            for i = find(nodos_grado_nulo).'          % Devuelve posiciones de esos nodos con grado 0
                j = randi(n);       % Elige un destino aleatorio, si es él mismo lo cambia
                if j == i, j = mod(j, n) + 1; end
                A0(i,j) = 1;
            end
        end

        % --- Normalización por filas ---
       % grado = sum(A0,2);
        W0 = A0 ./ grado;

        % --- PageRank ---
        G  = digraph(A0);            % Se construye el grafo dirigido con la matriz de adyacencia
        pr = centrality(G, 'pagerank', 'FollowProbability', damp);      % Cálculo de la centralidad de cada nodo
        pr = pr / sum(pr);          % Normalización a 1.

        % Recogida n Bank por cada tipo de régimen y réplica
        redes.(reg).A{r}  = A0;
        redes.(reg).W{r}  = W0;
        redes.(reg).pr{r} = pr;
    end
end

%% Generar banco de redes base
esc = struct('n',{},'trolls',{},'p',{}, 'lam',{}, ...
             'seed',{}, 'tag',{}, 'regimen',{}, 'loc',{}, 'fracTrolls',{});
c = 1;
for i_reg = 1:numel(regimen)
    reg = regimen{i_reg};
    p   = prob_regimen(i_reg);
    for j_prop_trolls = prop_trolls
        trolls_total = round(j_prop_trolls * n);
        for k_lam = lambdas
            for l_banda = 1:numel(bandas_pr)
                loc = bandas_pr{l_banda};
                tag_str = sprintf('n%d_trolls%d_%s_%s_lam%.2f', n, trolls_total, reg, loc, k_lam);
                esc(c) = struct( ...
                    'n',    n, ...
                    'trolls',  trolls_total, ...
                    'p',    p, ...
                    'lam',  k_lam, ...
                    'seed', semilla_base + c, ...
                    'tag',  tag_str, ...
                    'regimen', reg, ...
                    'loc', loc, ...
                    'fracTrolls', j_prop_trolls ...
                );
                c = c + 1;
            end
        end
    end
end
N = numel(esc);
fprintf('Generados %d escenarios\n', N);

%% Bucle principal
% Preasignación de memoria
totalFilas = N * num_replicas;
rows = cell(totalFilas, 17);   % Preasignación de filas y columnas
row_id = 1;

for i = 1:N
    s = esc(i);

    for r = 1:num_replicas
        
        % --- Red de la réplica r ---
        A  = redes.(s.regimen).A{r};
        W  = redes.(s.regimen).W{r};
        pr = redes.(s.regimen).pr{r};

        % --- Selección trolls ---
        trolls_total = s.trolls;
        if trolls_total > 0
            idxBand = seleccionar_trolls(pr, trolls_total, s.loc, 42);
            trolls = idxBand(1:min(trolls_total, length(idxBand)));
        else
            trolls = [];
        end
        normales = setdiff((1:s.n).', trolls(:));

        % --- x0 y lambdas ---
        x0 = x0_base; 
        x0(trolls) = -1;
        lambdas = zeros(s.n,1);
        lambdas(normales) = s.lam;
        
        % Numero de componentes conexas 
        A_bin  = spones(A) ~= 0;
        Gd = digraph(A_bin);
        [compStrong, sizes] = conncomp(Gd, 'Type', 'strong');
        nSCC = numel(sizes);    

        % --- Punto fijo FJ ---
        I = speye(s.n);
        Lam = spdiags(lambdas,0,s.n,s.n);
        M = (I - Lam*W);
        b = (I - Lam)*x0;
        x_star = M \ b;
    
        % --- Simulación ---
        rho_LW = abs(eigs(Lam*W, 1, 'lm'));
        rho = min(rho_LW, 1 - 1e-12);
        C0 = norm(x0 - x_star, Inf);
        k_needed = ceil(log(param.tol_fp / max(C0,eps)) / log(rho));
        k_cap = 10000;
        Tsim = max(param.iteraciones, min(k_needed, k_cap));

        X = simular_friedkin(W, x0, lambdas, Tsim, x_star, param.tol_fp); 
        z = X(:,end);
        zn = z(normales);

        % --- Tiempo de convergencia ---
        d_inf = vecnorm(X - x_star, Inf, 1);
        tconv = find(d_inf < param.tol_fp, 1, 'first');
        if isempty(tconv), tconv = NaN; end

        % --- Métricas ---
        rf  = max(z) - min(z);
        sd  = std(z);
        medianaFinal = median(z);
        meanFinal  = mean(z);

        propNeg = mean(zn < 0);
        medianaNorm = median(zn);

        Delta = z - z.';
        NDI = sum(sum(W .* (Delta.^2)));

        resid = norm(z - x_star, Inf);

        % --- Guardar fila ---
        rows{row_id,1} = s.tag;
        rows{row_id,2} = r;
        rows{row_id,3} = s.regimen;
        rows{row_id,4} = s.loc;
        rows{row_id,5} = s.fracTrolls;
        rows{row_id,6} = s.lam;

        rows{row_id,7} = rf;
        rows{row_id,8} = sd;
        rows{row_id,9} = medianaFinal;
        rows{row_id,10} = meanFinal;
        rows{row_id,11} = propNeg;
        rows{row_id,12} = medianaNorm;

        rows{row_id,13} = NDI;
        rows{row_id,14} = resid;
        rows{row_id,15} = tconv;
        rows{row_id,16}= rho;
        rows{row_id,17}= nSCC;
        row_id = row_id + 1;
    end
end

%% Exportar resultados
headers = {'tag','replica','regimen','loc','fracTrolls','lambda',...
           'rangoFinal','stdFinal','medianFinal','meanFinal',...
           'propNeg','medianNorm','NDI','resid','tconv','rho_lambdaW','nSCC'};
tabla = cell2table(rows,'VariableNames',headers);
writetable(tabla, fullfile(carpeta_resultados,'resultados.csv'));
save(fullfile(carpeta_resultados,'resultados.mat'),'tabla');

fprintf('Guardado: resultados.csv\n');
diary off;

%% ========================================================================
%  FUNCIONES AUXILIARES
%  ========================================================================

function idxBand = seleccionar_trolls(pr, n_trolls, nivel_banda, seed)
% Selecciona trolls según banda de centralidad (bajoPR, medioPR, altoPR)
    rng(seed, 'twister');
    n = length(pr);
    [~, sortIdx] = sort(pr, 'descend');
    
    tercil = floor(n/3);
    
    switch nivel_banda
        case 'bajoPR'
            band = sortIdx(2*tercil+1:end);
        case 'medioPR'
            band = sortIdx(tercil+1:2*tercil);
        case 'altoPR'
            band = sortIdx(1:tercil);
        otherwise
            error('Banda desconocida: %s', nivel_banda);
    end
    
    if n_trolls > length(band)
        warning('numTrolls (%d) > tamaño banda (%d). Se toman todos los de la banda.', ...
                n_trolls, length(band));
        idxBand = band;
    else
        perm = randperm(length(band), n_trolls);
        idxBand = band(perm);
    end
end

function X = simular_friedkin(W, x0, lambdas, n_iter, x_star, tol)
% Simula el modelo de Friedkin-Johnsen
% W: matriz de influencia (n x n)
% x0: opiniones iniciales (n x 1)
% lambdas: vector de prejuicios (n x 1)
% n_iter: número de iteraciones
% x_star: punto fijo teórico (para convergencia)
% tol: tolerancia de convergencia

    n = length(x0);
    X = zeros(n, n_iter+1);
    X(:,1) = x0;
    
    for t = 1:n_iter
        x_prev = X(:,t);
        x_actual = (1 - lambdas) .* x0 + lambdas .* (W * x_prev);
        X(:,t+1) = x_actual;
        
        % Convergencia anticipada
        if norm(x_actual - x_star, Inf) < tol
            X = X(:, 1:t+1);
            break;
        end
    end
end
