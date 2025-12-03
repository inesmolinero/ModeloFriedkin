%% main.m – Experimentos Friedkin–Johnsen
%===========================================================================
% Redes ER, simulación Friedkin–Johnsen, conectividad, centralidad,
% convergencia, polarización, clústeres y exportación de resultados.
% - Trolls: negativos
% - Colocación: determinística por centralidad (bajoPR/medioPR/altoPR)
% - x0 fijo en todos los escenarios
% - W mínima variación y 50 réplicas de redes por régimen
%===========================================================================

%% LIMPIEZA
clear; close all; clc;

%% 0.  OPINIONES INICIALES
Dir = fileparts(mfilename('fullpath'));
x0_raw_all  = readmatrix(fullfile(Dir, 'data', 'Twitter_Data_scores.csv'));
x0_raw_all  = x0_raw_all(:);        % convierte a vector columna
x0_raw_all  = x0_raw_all(~isnan(x0_raw_all));   % elimina nulos
assert(~isempty(x0_raw_all), 'Twitter_Data_scores.csv no tiene datos numéricos válidos.');

fprintf('Cargadas opiniones iniciales\n');

%% 1 · PARÁMETROS BASE
param = struct();
param.rango_opiniones  = [0, 1];
param.seed_base        = 123;          % semilla base
param.iteraciones      = 500;
param.tol_consenso     = 0.05;
param.tol_fp           = 1e-6;         % tolerancia a punto fijo x*

param.n = 300;

%% Construir x0_base a tamaño n y escalar al rango de opiniones
rng(param.seed_base+777,'twister');   

m = numel(x0_raw_all); % numero de valores en vector de opiniones x0

% Si el número de valores es mayor, se toman n de esos m. En otro caso,se
% usa repmat para superar n y de ese vector extendido se toman n valores.
if m >= param.n
    xv = x0_raw_all(randperm(m, param.n));
else
    xv = repmat(x0_raw_all, ceil(param.n/m), 1);
    xv = xv(randperm(numel(xv), param.n));
end

% Escalado al rango param.rango_opiniones
xmin = min(xv); xmax = max(xv);
if xmax > xmin
    xv = (xv - xmin) / (xmax - xmin);  % [0,1]
else
    xv = zeros(param.n,1);
end
x0_base = param.rango_opiniones(1) + xv*diff(param.rango_opiniones);
fprintf('x0_base generado: rango [%.3f, %.3f]\n', min(x0_base), max(x0_base));

%% 2 · ESCENARIOS (fracción total de trolls, lambda, regimen de unión y bandas de centralidad)
% Regímenes ER por densidad p = c * log(n)/n
pStar  = log(param.n)/param.n;
regimenes   = {'desconectada','umbral','fuerte'};
pPorRegimen = [0.5, 1.0, 3.0] * pStar;                 % 3 regímenes
pPorRegimen = min(0.99, max(1e-6, pPorRegimen));       % seguridad
locBands       = {'bajoPR','medioPR','altoPR'};        % terciles de centralidad

% Proporción de trolls
propTrollsVals = [0.00, 0.10, 0.20, 0.30];             % 0%, 10%, 20% y 30% de trolls totales

% Valores lambda (prejuicio)
lamVals        = [0.2, 0.35, 0.5, 0.65, 0.8];          % lambda para NO tercos
               
% Construcción de escenarios
esc = struct('n',{},'trolls',{},'p',{}, 'lam',{}, ...
             'seed',{}, 'tag',{}, 'regimen',{}, 'loc',{}, 'fracTrolls',{});
c = 1;

for ireg = 1:numel(regimenes)
    reg = regimenes{ireg};
    p   = pPorRegimen(ireg);
    for fracTrolls = propTrollsVals
        tTot = round(fracTrolls * param.n);
        for lam = lamVals
            for ib = 1:numel(locBands)
                loc = locBands{ib};
                tagStr = sprintf('n%d_trolls%d_%s_%s_lam%.2f', param.n, tTot, reg, loc, lam);
                esc(c) = struct( ...
                    'n',    param.n, ...
                    'trolls',  tTot, ...
                    'p',    p, ...
                    'lam',  lam, ...
                    'seed', param.seed_base + c, ...
                    'tag',  tagStr, ...
                    'regimen', reg, ...
                    'loc', loc, ...
                    'fracTrolls', fracTrolls ...
                );
                c = c + 1;
            end
        end
    end
end
N = numel(esc);
fprintf('Generados %d escenarios\n', N);

%% 3 · RESULTADOS: carpeta y log
if ~isempty(mfilename)
    scriptDir = fileparts(mfilename('fullpath'));
else
    scriptDir = pwd;  
end
resultsRoot = scriptDir;
tsStr       = char(datetime('now','Format','yyyy-MM-dd_HHmmss'));
resultsDir  = fullfile(resultsRoot, ['resultados_' tsStr]);

if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
diary(fullfile(resultsDir,'log.txt'));  
fprintf('Directorio de resultados: %s\n', resultsDir);

%% 4 · Banco de R redes por régimen para reproducibilidad
R = 50;          % R réplicas por cada régimen de densidad
damp = 0.9;      % damping factor: parámetro de salto aleatorio

fprintf('Generando banco de %d redes por régimen...\n', R);
rng(param.seed_base, 'twister');

% fija para las redes base
Bank = struct();   % Banco de R réplicas de A, W y pr (vectores pageRank)
for ireg = 1:numel(regimenes)
    reg  = regimenes{ireg};
    pReg = pPorRegimen(ireg);

    Bank.(reg).A  = cell(R,1);
    Bank.(reg).W  = cell(R,1);
    Bank.(reg).pr = cell(R,1);

    for r = 1:R
        rng(param.seed_base + 1000*ireg + r, 'twister');  % semilla distinta por regimen y réplica

        % --- ER dirigido sin autolazos ---
        A0 = rand(param.n) < pReg;        % ER dirigido G(n,p) con probabilidad de arista pReg
        A0(1:param.n+1:end) = 0;          % Anular la diagonal para evitar autolazos

        % --- Garantiza al menos 1 salida por nodo ---
        outdeg = sum(A0,2);         % Suma de grados de cada fila
        dang = (outdeg == 0);       % Marca nodos con grado = 0
        if any(dang)
            for i = find(dang).'          % Devuelve posiciones de esos nodos con grado 0
                j = randi(param.n);       % Elige un destino aleatorio, si es él mismo lo cambia
                if j == i, j = mod(j, param.n) + 1; end
                A0(i,j) = 1;
            end
            outdeg = sum(A0,2);
        end

        % --- Normalización por filas ---
        W0 = A0 ./ outdeg;

        % --- PageRank ---
        G  = digraph(A0);            % Se construye el grafo dirigido con la matriz de adyacencia
        pr = centrality(G, 'pagerank', 'FollowProbability', damp);      % Cálculo de la centralidad de cada nodo
        pr = pr / sum(pr);          % Normalización a 1.

        % --- Recogida n Bank por cada tipo de régime y réplica ---
        Bank.(reg).A{r}  = A0;
        Bank.(reg).W{r}  = W0;
        Bank.(reg).pr{r} = pr;
    end
end

%% 5 · BUCLE PRINCIPAL — UNA FILA POR RÉPLICA (SIN MEDIAS)
% Preasignación de memoria
totalFilas = N * R;
rows = cell(totalFilas, 17);   % Preasignamos todas las filas y columnas
row_id = 1;

for i = 1:N
    s = esc(i);

    for r = 1:R
        
        % --- Red de la réplica r ---
        A  = Bank.(s.regimen).A{r};
        W  = Bank.(s.regimen).W{r};
        pr = Bank.(s.regimen).pr{r};

        % --- Selección trolls ---
        tTot = s.trolls;
        if tTot > 0
            idxBand = pickByCentralityBands(pr, tTot, s.loc, 42);
            trolls = idxBand(1:min(tTot, length(idxBand)));
        else
            trolls = [];
        end
        normales = setdiff((1:s.n).', trolls(:));

        % --- x0 y lambdas ---
        x0 = x0_base; 
        x0(trolls) = -1;
        lambdas = zeros(s.n,1);
        lambdas(normales) = s.lam;

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

        X = simularFriedkin(W, x0, lambdas, Tsim, x_star, param.tol_fp); 
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
        P2 = mean(z.^2);
        P4 = sum(abs(z));

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
        rows{row_id,14} = P2;
        rows{row_id,15} = P4;
        rows{row_id,16} = resid;
        rows{row_id,17} = tconv;

        row_id = row_id + 1;
    end
end

%% 6 · Exportar CSV FULL REPLICAS
headers = {'tag','replica','regimen','loc','fracTrolls','lambda',...
           'rangoFinal','stdFinal','medianFinal','meanFinal',...
           'propNeg','medianNorm','NDI','P2','P4','resid','tconv'};
T = cell2table(rows,'VariableNames',headers);
writetable(T, fullfile(resultsDir,'resultados.csv'));
save(fullfile(resultsDir,'resultados.mat'),'T');

fprintf('Guardado: todas_las_replicas.csv\n');
diary off;

%% ========================================================================
%  FUNCIONES AUXILIARES
%  ========================================================================

function idxBand = pickByCentralityBands(pr, numTrolls, bandName, seed)
% Selecciona trolls según banda de centralidad (bajoPR, medioPR, altoPR)
    rng(seed, 'twister');
    n = length(pr);
    [~, sortIdx] = sort(pr, 'descend');
    
    tercil = floor(n/3);
    
    switch bandName
        case 'bajoPR'
            band = sortIdx(2*tercil+1:end);
        case 'medioPR'
            band = sortIdx(tercil+1:2*tercil);
        case 'altoPR'
            band = sortIdx(1:tercil);
        otherwise
            error('Banda desconocida: %s', bandName);
    end
    
    if numTrolls > length(band)
        warning('numTrolls (%d) > tamaño banda (%d). Se toman todos los de la banda.', ...
                numTrolls, length(band));
        idxBand = band;
    else
        perm = randperm(length(band), numTrolls);
        idxBand = band(perm);
    end
end

function X = simularFriedkin(W, x0, lambdas, T, x_star, tol)
% Simula el modelo de Friedkin-Johnsen
% W: matriz de influencia (n x n)
% x0: opiniones iniciales (n x 1)
% lambdas: vector de prejuicios (n x 1)
% T: número de iteraciones
% x_star: punto fijo teórico (para convergencia)
% tol: tolerancia de convergencia

    n = length(x0);
    X = zeros(n, T+1);
    X(:,1) = x0;
    
    for t = 1:T
        x_prev = X(:,t);
        x_new = (1 - lambdas) .* x0 + lambdas .* (W * x_prev);
        X(:,t+1) = x_new;
        
        % Convergencia anticipada
        if norm(x_new - x_star, Inf) < tol
            X = X(:, 1:t+1);
            break;
        end
    end
end

function [isStrong, period, scc_periods, nSCC] = digraph_props(A)
% Propiedades dirigidas: SCCs (fuertemente conexo), irreducible, periodo.

    % A como adyacencia 0/1
    A  = spones(A) ~= 0;
    Gd = digraph(A);

    % compStrong: etiqueta 1..K por nodo
    % sizes: tamaños de cada SCC (vector de longitud K)
    [compStrong, sizes] = conncomp(Gd, 'Type', 'strong');
    nSCC = numel(sizes);            % nº de SCC (ESCALAR)

    isStrong = (nSCC == 1);

    % Preasigna vector de períodos por SCC
    scc_periods = NaN(nSCC,1);
    for k = 1:nSCC
        nodes = find(compStrong == k);
        if isempty(nodes)
            scc_periods(k) = NaN;
            continue
        end
        Ak = A(nodes, nodes);
        scc_periods(k) = scc_period(Ak);
    end

    if isStrong
        period = scc_periods(1);
    else
        period = NaN;
    end
end

function d = scc_period(Ak)
% Periodo d de una SCC dirigida (mcd de longitudes de ciclos).
    Ak = spones(Ak) ~= 0;
    n  = size(Ak,1);
    if n == 0
        d = NaN; 
        return
    end

    % BFS desde un nodo (1) para obtener distancias topológicas
    dist = inf(n,1); 
    dist(1) = 0; 
    Q = 1;
    while ~isempty(Q)
        u = Q(1); 
        Q(1) = [];
        nb = find(Ak(u,:));
        for v = nb
            if isinf(dist(v))
                dist(v) = dist(u) + 1; 
                Q(end+1) = v; %#ok<AGROW>
            end
        end
    end

    % mcd de (dist(u)+1 - dist(v)) sobre todas las aristas u->v alcanzables
    g = 0; 
    [u,v] = find(Ak);
    for e = 1:numel(u)
        du = dist(u(e)); 
        dv = dist(v(e));
        if isfinite(du) && isfinite(dv)
            g = gcd(g, abs(du + 1 - dv));
        end
    end
    d = max(g,1);
end