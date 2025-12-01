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
x0_raw_all  = x0_raw_all(:);        %convierte a vector columna
x0_raw_all  = x0_raw_all(~isnan(x0_raw_all));   %elimina nulos
assert(~isempty(x0_raw_all), 'Twitter_Data_scores.csv no tiene datos numéricos válidos.'); %si el vector se queda en nulo -> error

fprintf('Cargadas opiniones iniciales\n');

%% 1 · PARÁMETROS BASE
param = struct();
param.rango_opiniones  = [0, 1];
param.seed_base        = 123;          % semilla base
param.iteraciones      = 500;
param.tol_consenso     = 0.05;
param.tol_fp        = 1e-6;    % tolerancia a punto fijo x*

param.n = 300;

%% Construir x0_base a tamaño n y escalar al rango de opiniones
rng(param.seed_base+777,'twister');   

m = numel(x0_raw_all); % numero de valores en vector de opiniones x0

%Si el número de valores es mayor, se toman n de esos m. En otro caso,se
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
lamVals        = [ 0.2, 0.35, 0.5, 0.65, 0.8];          % lambda para NO tercos
               
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

%% 3 · RESULTADOS: carpeta y log (carpeta del script)
if ~isempty(mfilename)
    scriptDir = fileparts(mfilename('fullpath'));
else
    scriptDir = pwd;  
end
resultsRoot = scriptDir;
tsStr       = datestr(now,'yyyy-mm-dd_HHMMSS');
resultsDir  = fullfile(resultsRoot, ['resultados_' tsStr]);

if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
diary(fullfile(resultsDir,'log.txt'));  
fprintf('Directorio de resultados: %s\n', resultsDir);

%% === Config de guardado por réplica ===
saveRuns   = true;                                      % activar/desactivar guardado
runsDir    = fullfile(resultsDir, 'runs');              % carpeta raíz de runs
if saveRuns && ~exist(runsDir, 'dir'), mkdir(runsDir); end

snap_count = 12;    % nº de snapshots (tiempos log-espaciados)
sample_m   = 60;    % nº de agentes en la muestra de trayectorias

% índice de runs (lo rellenamos dentro del bucle y lo exportamos al final)
runsIndex = table('Size',[0 3], ...
    'VariableTypes', {'string','double','string'}, ...
    'VariableNames', {'tag','rep','path'}); 

%% 4 · Banco de R redes por régimen para reproducibilidad
R= 50;          % R réplicas por cada régimen de densidad
damp = 0.9;    % damping factor: parámetro de salto aleatorio

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
        A0 = rand(param.n) < pReg;        %ER dirigido G(n,p) con probabilidad de arista pReg
        A0(1:param.n+1:end) = 0;          % Anular la diagonal para evitar autolazos

        % --- Garantiza al menos 1 salida por nodo ---
        outdeg = sum(A0,2);         %Suma de grados de cada fila
        dang = (outdeg == 0);       %Marca nodos con grado = 0
        if any(dang)
            for i = find(dang).'          %Devuelve posiciones de esos nodos con grado 0
                j = randi(param.n);       %Elige un destino aleatorio, si es él mismo lo cambia
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

%% 5 · BUCLE PRINCIPAL con R=50 réplicas
clear MAgg;
cnt = 1;

% === Métricas a promediar (numéricas) ===
numFields = {'rho_W', 'rho_LambdaW', ...
             'rangoFinal', 'stdFinal', ...
             'convTime', ...
             'P_norm_vs_trolls', ...
             'periodLW', 'periodW', 'nSCC_W', 'nSCC_LW', ...
             'NDI', 'P2', 'P4', ...
             'avgCentTrolls', 'avgCentNormals', ...
             'tconv', 'kneed', 'propNegNormals', 'medianFinalNormals', 'medianFinal', 'resid_fp'};

boolFields = {'stochastic', 'substochastic', ...
              'isStrongLW', 'isPrimLW', 'isStrongW', ...
              'consenso'};

% === Estructura agregada (una fila por escenario) ===
MAgg(N,1) = struct( ...
  'tag','', 'n',[], 'trolls',[],  'p',[], 'lam',[], ...
  'regimen','', 'loc','', 'fracTrolls',[], ...
  'mean_rho_W',[], 'mean_rho_LW',[],  ...
  'mean_rangoFinal',[], 'mean_stdFinal',[], ...
  'mean_convTime',[],  ...
  'mean_P_norm_vs_trolls',[], ...
  'mean_propNegNormals',[], 'mean_medianFinalNormals',[], 'mean_resid_fp',[], ...
  'mean_medianFinal', [], 'mean_meanFinal', [], ...
  'std_convTime',[], 'std_rangoFinal',[], 'std_stdFinal',[], ...
  'prop_stochastic',[], 'prop_substochastic',[], ...
  'prop_isStrongLW',[],  ...
  'mean_nSCC_W',[], 'mean_nSCC_LW',[], ...                 
  'mean_periodW',[], 'mean_periodLW',[], ...                
  'prop_isStrongW',[], 'prop_isPrimLW', [], ...                                   
  'mean_nSCC_W_over_n',[], 'mean_nSCC_LW_over_n',[] ...     
);


for i = 1:N
    s = esc(i);

    % Acumuladores por réplica
    acc = struct();
    for f = numFields,  acc.(f{1})  = nan(R,1);   end
    for f = boolFields, acc.(f{1})  = false(R,1); end

    for r = 1:R
        % --- Red de la réplica r ---
        A  = Bank.(s.regimen).A{r};
        W  = Bank.(s.regimen).W{r};
        pr = Bank.(s.regimen).pr{r};

        % --- Comprobación estocasticidad ---
        is_stoch = all(abs(sum(W,2) - 1) < 1e-12);

        % --- Selección trolls por banda ---
        tTot = s.trolls;
        if tTot > 0
            idxBand = pickByCentralityBands(pr, tTot, s.loc, 42);
        else
            idxBand = [];
        end
        trolls = idxBand(1:tTot);
        normales   = setdiff((1:s.n).', trolls.');
        % --- x0 y lambdas ---
        x0 = x0_base;  x0(trolls) = -1; 
        lambdas = zeros(s.n,1);  lambdas(normales) = s.lam;

        % --- ΛW subestocástica ---
        LambdaW  = diag(lambdas) * W;
        is_subst = all(sum(LambdaW,2) <= 1 + 1e-12);

        % --- Espectros ---
        rho_W        = max(abs(eig(full(W))));
        spLambdaW    = spdiags(lambdas,0,s.n,s.n) * W;
        rho_LambdaW  = abs(eigs(spLambdaW, 1, 'lm'));

        % --- Propiedades de ΛW ---
        spLW = spones(spLambdaW);
        [isStrongLW, periodLW, scc_periods_LW, nSCC_LW] = digraph_props(spLW);
        isPrimLW = isStrongLW && (periodLW==1);

        % --- Propiedades de W ---
        spW = spones(W);
        [isStrongW, periodW, scc_periods_W, nSCC_W] = digraph_props(spW);
        isPrimW = isStrongW && (periodW==1);

        % --- Punto fijo analítico de FJ  ---
        I   = speye(s.n);
        Lam = spdiags(lambdas,0,s.n,s.n);
        M   = (I - Lam*W);
        b   = (I - Lam)*x0;
        x_star = M \ b;
        
        % --- Estimación de iteraciones necesarias ---
        rho = min(rho_LambdaW, 1 - 1e-12);           % margen numérico
        C0  = norm(x0 - x_star, Inf);
        k_needed = ceil( log(param.tol_fp / max(C0,eps)) / log(rho) );
        k_cap    = 10000;                             % tope de seguridad
        Tsim     = max(param.iteraciones, min(k_needed, k_cap));
        
        % --- Simulación con parada temprana respecto a x* ---
        X = simularFriedkin(W, x0, lambdas, Tsim, x_star, param.tol_fp); 
        z = X(:,end);
        zn = z(normales);

        % --- Tiempo de convergencia a x* ---
        d_inf = vecnorm(X - x_star, Inf, 1);
        tconv = find(d_inf < param.tol_fp, 1, 'first');
        converge = ~isempty(tconv);
        if isempty(tconv), tconv = NaN; end
        
        % --- Residuos ---
        resid = norm(z - x_star, Inf);


        % --- Métricas básicas ---
        rf = max(z) - min(z);
        sd = std(z);
        medianFinal = median(z);
        meanFinal = mean(z);
        propNegNorm = mean(zn < 0);            
        medianFinalNorm  = median(zn); 

        % --- Métricas de polarización ---
        x = z(:);                 % vector columna

        % Diferencias por pares (matriz n x n)
        Delta = x - x.';          % (x_i - x_j)

        % 1) NDI
        %    Sum_{i,j} w_ij (x_i - x_j)^2
        NDI = sum(sum(W .* (Delta.^2)));

        % 4) P2 (1/|V|) * sum_i x_i^2
        P2 = mean(x.^2);              

        % 6) P4 sum_i |x_i|
        P4 = sum(abs(x));       

       
        % --- Centralidades medias ---

        idxT= trolls;
        if isempty(idxT), avgCentT = NaN; else, avgCentT = mean(pr(idxT),'omitnan'); end
        avgCentN = mean(pr(normales),'omitnan');

        labelAgentes = zeros(s.n,1);
        labelAgentes(normales)=1;
        labelAgentes(trolls)=2;
       
        is_cons = range(z) < param.tol_consenso;


 %% === Guardado por réplica: finales, stats temporales, snapshots y muestra ===
        if saveRuns
            % 1) Estadísticos temporales
            T = size(X,2);
            rango_it = max(X,[],1) - min(X,[],1);
            mean_it  = mean(X,1);
            mediana_it = median(X, 1);
            std_it   = std(X,0,1);
            q05_it   = quantile(X,0.05,1);
            q50_it   = quantile(X,0.50,1);
            q95_it   = quantile(X,0.95,1);
           
            % 2) Snapshots log-espaciados
            tgrid  = unique( max(1, min(T, round(logspace(0, log10(max(T,1)), snap_count)))) );
            X_snap = X(:, tgrid);
            K_snap = nan(size(tgrid));
        
            % 3) Muestra de agentes (estratificada por tipo y centralidad)
            idx_sample = sample_by_strata(pr, labelAgentes, sample_m);
            X_sample   = X(idx_sample, :);
        
            % 4) Empaquetar
            out = struct();
        
            % Parámetros del escenario
            out.params = struct('tag', s.tag, 'n', s.n, 'trolls', s.trolls, ...
                                'p', s.p, 'lam', s.lam, 'regimen', s.regimen, ...
                                'loc', s.loc, 'fracTrolls', s.fracTrolls);
        
            % Índices 
            out.indices = struct('normales', normales,  ...
                                 'trolls', trolls, 'idx_sample', idx_sample);
            
            % Resultados finales
            out.final = struct( ...
                        'z', z, ...
                        'x_star', x_star, ...              % punto fijo
                        'resid_fp', resid, ...             % ||z - x*||inf
                        'tconv', tconv, ...
                        'consenso', is_cons, ...          % = is_cons
                        'rf', rf, ...      % rango final
                        'sd', sd, ...      % std final
                        'NDI', NDI, ...    % polarización 
                        'P2', P2, ...     
                        'P4', P4 ...       
                    );
            %Propiedades algebraicas
            out.spectral   = struct('rho_LW', rho_LambdaW, 'rho_W', rho_W);

            % Estadísticos temporales
            out.time_stats = struct('t', 1:T, 'rango', rango_it, 'mean', mean_it, 'median',mediana_it, 'std', std_it, ...
                                    'q05', q05_it, 'q50', q50_it, 'q95', q95_it);
        
            % Snapshots
            out.snapshots = struct('tgrid', tgrid, 'X', X_snap, 'K_snap', K_snap);
        
            % Muestra de trayectorias
            out.sample = struct('idx', idx_sample, 'X', X_sample);
        
            % Propiedades de grafos
            out.graph = struct('isStrongW', isStrongW, 'periodW', periodW, 'nSCC_W', nSCC_W, ...
                               'isStrongLW', isStrongLW, 'periodLW', periodLW, 'nSCC_LW', nSCC_LW, ...
                               'isPrimLW', isPrimLW);
           
            % Centralidades
            out.centrality = struct('avgCentTrolls', avgCentT, 'avgCentNormals', avgCentN, 'pr', pr);
        
            % 5) Guardar archivo
            scenDir = fullfile(runsDir, s.tag);
            if ~exist(scenDir, 'dir'), mkdir(scenDir); end
            runFile = fullfile(scenDir, sprintf('rep_%03d.mat', r));
            save(runFile, 'out', '-v7.3');
        
            % Añadir al índice
            runsIndex = [runsIndex; {string(s.tag), r, string(runFile)}]; 
        end
        % ======= Acumuladores =======
        acc.rho_W(r)    = rho_W;
        acc.rho_LambdaW(r)  = rho_LambdaW;

        acc.rangoFinal(r)   = rf;
        acc.stdFinal(r)     = sd;
        acc.resid_fp(r) = resid;
        acc.medianFinal(r) = medianFinal;
        acc.meanFinal(r)= meanFinal;
        acc.NDI(r)   = NDI;
        acc.P2(r)    = P2;
        acc.P4(r)    = P4;

        acc.propNegNormals(r)      = propNegNorm;
        acc.medianFinalNormals(r)  = medianFinalNorm;


        acc.avgCentTrolls(r)  = avgCentT;
        acc.avgCentNormals(r) = avgCentN;

        acc.convTime(r)     = ifelse(converge, tconv, NaN); % guarda NaN si no hay consenso

        
        acc.consenso(r)         = is_cons;
        acc.P_norm_vs_trolls(r)   = P_nt;
        acc.stochastic(r)    = is_stoch;
        acc.substochastic(r) = is_subst;

        acc.isStrongW(r)  = isStrongW;
        acc.periodW(r)    = periodW;
        acc.nSCC_W(r)     = nSCC_W;
        
        acc.isStrongLW(r) = isStrongLW;
        acc.periodLW(r)   = periodLW;
        acc.nSCC_LW(r)    = nSCC_LW;
        acc.isPrimLW(r)   = isPrimLW;
        
        acc.tconv(r)  = tconv;
        acc.kneed(r)  = k_needed;
        
    end

    % ======= Agregación por escenario =======
    MAgg(cnt).tag        = s.tag;     MAgg(cnt).n   = s.n;
    MAgg(cnt).trolls        = s.trolls;     
    MAgg(cnt).p          = s.p;       MAgg(cnt).lam = s.lam;
    MAgg(cnt).regimen    = s.regimen; MAgg(cnt).loc = s.loc;
    MAgg(cnt).fracTrolls = s.fracTrolls;

    % Medias numéricas
    MAgg(cnt).mean_rho_emp_W   = mean(acc.rho_W,'omitnan');
    MAgg(cnt).mean_rho_LW      = mean(acc.rho_LambdaW,'omitnan');

    MAgg(cnt).mean_rangoFinal  = mean(acc.rangoFinal,'omitnan');
    MAgg(cnt).mean_stdFinal    = mean(acc.stdFinal,'omitnan');
    MAgg(cnt).mean_convTime    = mean(acc.convTime,'omitnan');
    MAgg(cnt).mean_NDI = mean(acc.NDI, 'omitnan');
    MAgg(cnt).mean_P2  = mean(acc.P2,  'omitnan');
    MAgg(cnt).mean_P4  = mean(acc.P4,  'omitnan');

    MAgg(cnt).mean_P_norm_vs_trolls = mean(acc.P_norm_vs_trolls,'omitnan');
    
    MAgg(cnt).std_convTime   = std(acc.convTime, 'omitnan');
    MAgg(cnt).std_rangoFinal = std(acc.rangoFinal, 'omitnan');
    MAgg(cnt).std_stdFinal   = std(acc.stdFinal, 'omitnan');
    MAgg(cnt).mean_propNegNormals      = mean(acc.propNegNormals,     'omitnan');
    MAgg(cnt).mean_medianFinalNormals  = mean(acc.medianFinalNormals, 'omitnan');
    MAgg(cnt).mean_resid_fp            = mean(acc.resid_fp,           'omitnan');
    MAgg(cnt).mean_medianFinal = mean(acc.medianFinal, 'omitnan');
    MAgg(cnt).mean_meanFinal =mean(acc.meanFinal, 'omitnan');

    % Proporciones booleanas
    MAgg(cnt).prop_stochastic      = mean(acc.stochastic);
    MAgg(cnt).prop_substochastic   = mean(acc.substochastic);
    MAgg(cnt).prop_isStrongLW      = mean(acc.isStrongLW);
    MAgg(cnt).prop_isPrimLW        = mean(acc.isPrimLW);
    MAgg(cnt).prop_consenso        = mean(acc.consenso);
    
    % Periodo y SCC (ΛW)
    MAgg(cnt).mean_nSCC_W  = mean(acc.nSCC_W, 'omitnan');
    MAgg(cnt).mean_nSCC_LW  = mean(acc.nSCC_LW,  'omitnan');

    idxStrongW  = acc.isStrongW  == true;
    idxStrongLW = acc.isStrongLW == true;
    MAgg(cnt).mean_periodW  = ifelse(any(idxStrongW),  mean(acc.periodW(idxStrongW),  'omitnan'), NaN);
    MAgg(cnt).mean_periodLW = ifelse(any(idxStrongLW), mean(acc.periodLW(idxStrongLW),'omitnan'), NaN);
    
    % Por si quieres ver proporciones (útil para interpretar NaN)
    MAgg(cnt).prop_isStrongW  = mean(acc.isStrongW);
    MAgg(cnt).prop_isStrongLW = mean(acc.isStrongLW);
    MAgg(cnt).prop_isPrimLW   = mean(acc.isPrimLW);
    MAgg(cnt).mean_periodLW = mean(acc.periodLW, 'omitnan');
    MAgg(cnt).mean_nSCC_W_over_n  = mean(acc.nSCC_W  / s.n, 'omitnan');
    MAgg(cnt).mean_nSCC_LW_over_n = mean(acc.nSCC_LW / s.n, 'omitnan');
    
    MAgg(cnt).mean_convTime    = mean(acc.tconv, 'omitnan');
    MAgg(cnt).mean_k_needed    = mean(acc.kneed, 'omitnan');

    cnt = cnt + 1;
end

% Helper
function y = ifelse(cond, a, b)
    if cond, y = a; else, y = b; end
end
%% 7 · Exportar índice de runs
if saveRuns && exist('runsIndex','var') && ~isempty(runsIndex)
    writetable(runsIndex, fullfile(runsDir, 'runs_index.csv'));
end

%% 8 · EXPORTAR RESULTADOS
TAgg = struct2table(MAgg);
writetable(TAgg, fullfile(resultsDir,'summary_metrics_AGG.csv'));
save(fullfile(resultsDir,'summary_metrics.mat'),'TAgg');
disp(TAgg);


%% Funciones auxiliares
function [isStrong, period, scc_periods, nSCC] = digraph_props(A)
% Propiedades dirigidas: SCCs (fuertemente conexo), irreducible, periodo.

    % A como adyacencia 0/1
    A  = spones(A) ~= 0;
    Gd = digraph(A);

    % compStrong: etiqueta 1..K por nodo
    % binsizes  : tamaños de cada SCC (vector de longitud K)
    [compStrong, sizes] = conncomp(Gd, 'Type', 'strong');
    nSCC = numel(sizes);            % nº de SCC (ESCALAR)

    isStrong = (nSCC == 1);
    %isIrred  = isStrong;            % para A>=0: irreducible ⇔ una SCC

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

function idx = sample_by_strata(pr, labelAgentes, m)
% Devuelve 'm' índices estratificados por tipo (1 normal, 2 t+, 3 t−)
% y por centralidad (baja/media/alta) dentro de cada tipo.

    n   = numel(pr);
    lbl = labelAgentes(:);   % columna
    pr  = pr(:);             % columna

    tipos = [1,2,3];
    nbins = 3;
    target_per_bin = max(1, round(m / (numel(tipos)*nbins)));

    idx = zeros(0,1);        % *** columna vacía (0x1) ***

    for t = tipos
        S = find(lbl == t);
        if isempty(S), continue; end

        % Quantiles para estratos
        q = quantile(pr(S), linspace(0,1,nbins+1));

        for b = 1:nbins
            % intervalo cerrado en extremos -> luego hacemos unique
            inbin = S(pr(S) >= q(b) & pr(S) <= q(b+1));
            if isempty(inbin), continue; end

            k = min(target_per_bin, numel(inbin));
            pick = randsample(inbin, k);   % Kx1
            idx  = [idx; pick];            % concatena columnas
        end
    end

    % Quitar duplicados y asegurar columna
    idx = unique(idx, 'stable');
    idx = idx(:);

    % Ajustar a tamaño m
    if numel(idx) > m
        idx = idx(1:m);
    elseif numel(idx) < m
        rest = setdiff((1:n).', idx);      % *** columna ***
        addk = min(m - numel(idx), numel(rest));
        if addk > 0
            idx = [idx; randsample(rest, addk)];   % concat columna con columna
        end
    end
end
