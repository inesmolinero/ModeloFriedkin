%% main.m – Experimentos Friedkin–Johnsen
%===========================================================================
% Generación de redes ER, simulación Friedkin–Johnsen,
% métricas espectrales, conectividad, primitividad,
% centralidad, convergencia, figuras y tablas.
%===========================================================================

%% 0 · LIMPIEZA
clear; close all; clc;

%% 1 · PARÁMETROS BASE
param.rango_opiniones  = [-0.7, 0.7];
param.seed_base        = 123;
param.iteraciones      = 1000;
param.tol_consenso     = 0.05;
param.tol_polarizacion = 0.60;

%% 2 · DEFINIR ESCENARIOS
sizes       = [20, 50, 100];
propNegVals = [0.0, 0.15, 0.3];  % proporción de trolls negativos
propPosVals = [0.0, 0.15, 0.3];  % proporción de trolls positivos
pVals       = [0.1, 0.3, 0.7];

% Construir lista de escenarios
esc = struct('n', {}, 'neg', {}, 'pos', {}, 'p', {}, 'seed', {}, 'tag', {});
c = 1;
for n = sizes
    for propNeg = propNegVals
        for propPos = propPosVals
            neg = round(propNeg * n);
            pos = round(propPos * n);
            if neg + pos >= n, continue; end
            for p = pVals
                tagStr = sprintf('n%d_neg%d_pos%d_p%.2f', n, neg, pos, p);
                esc(c) = struct(...
                    'n',    n, ...
                    'neg',  neg, ...
                    'pos',  pos, ...
                    'p',    p, ...
                    'seed', param.seed_base + c, ...
                    'tag',  tagStr ...
                );
                c = c + 1;
            end
        end
    end
end

N = numel(esc);

%% 3 · PREPARAR RESULTADOS
resultsRoot = 'C:/Users/HP/Documents/4 mates/TFG/modelomatlab';
tsStr       = datestr(now,'yyyy-mm-dd_HHMMSS');
resultsDir  = fullfile(resultsRoot, ['resultados_' tsStr]);
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
diary(fullfile(resultsDir,'log.txt'));

%% 3.1 · PREALOCAR STRUCT DE RESULTADOS
M(N,1) = struct(...
  'tag',           '', ...
  'n',             [],  ...
  'neg',           [],  ...
  'pos',           [],  ...
  'p',             [],  ...
  'irreducible',   false,...
  'periodo',       NaN, ...
  'primitiva',     false,...
  'rho_emp',       [],  ...
  'rho_teor',      [],  ...
  'norma2',        [],  ...
  'rangoFinal',    [],  ...
  'stdFinal',      [],  ...
  'polarizado',    false,...
  'avgCentTrolls', [],  ...
  'avgCentNormals',[],  ...
  'convTime',      NaN, ...
  'total_trolls',  [],  ...
  'prop_trolls_pos',[]  ...
);
cnt = 1;

%% 4 · BUCLE PRINCIPAL
for i = 1:N
    s = esc(i);
    rng(s.seed, 'twister');

    % 4.1 · Generar ER no dirigido
    A = rand(s.n) < s.p;
    A = A | A.';

    % 4.2 · Normalizar filas → W
    W = A ./ max(sum(A,2),1);

    % 4.3 · Conectividad e irreducibilidad
    A_test = double((W~=0) | (W.'~=0));
    irr    = (max(conncomp(graph(A_test))) == 1);

    % 4.4 · Periodicidad y Primitividad
    d    = NaN;
    prim = irr;

    % 4.5 · Radio espectral y norma
    ev      = eig(W);
    rho_emp = max(abs(ev));
    norma2  = norm(W);
    rho_teor = s.p * (s.n - 1);

    % 4.6 · Opiniones iniciales
    x0    = zeros(s.n,1);
    ids   = randperm(s.n);
    trolls_neg = ids(1:s.neg);
    trolls_pos = ids(s.neg+1:s.neg+s.pos);
    normales   = setdiff(1:s.n, [trolls_neg, trolls_pos]);
    x0(trolls_neg) = -1;
    x0(trolls_pos) =  1;
    x0(normales)   = param.rango_opiniones(1) + diff(param.rango_opiniones) * rand(numel(normales),1);

    % 4.7 · Simulación Friedkin–Johnsen
    lambdas = zeros(s.n,1);
    lambdas(normales) = 0.5 + 0.3*rand(numel(normales),1);
    X = simularFriedkin(W, x0, lambdas, param.iteraciones);

    % 4.8 · Métricas finales
    rf  = max(X(:,end)) - min(X(:,end));
    sd  = std(X(:,end));
    pol = (rf > 0.8) && (sd > 0.3);

    % 4.9 · Tiempo de convergencia
    rango_it = max(X,[],1) - min(X,[],1);
    tconv    = find(rango_it < param.tol_consenso, 1);
    if isempty(tconv), tconv = NaN; end

    % 4.10 · Centralidad
    [V,D]   = eig(W.');
    [~, idx]= max(abs(diag(D)));
    pi_vec  = abs(real(V(:,idx)));
    pi_vec  = pi_vec / sum(pi_vec);
    avgCentT = mean(pi_vec([trolls_neg, trolls_pos]));
    avgCentN = mean(pi_vec(normales));

    % 4.11 · Guardar .mat individual
    save(fullfile(resultsDir, [safeName(s.tag) '.mat']), 'W','x0','lambdas','X');

    % 4.12 · Almacenar en M
    M(cnt) = struct(...
      'tag',            s.tag,     ...
      'n',              s.n,       ...
      'neg',            s.neg,     ...
      'pos',            s.pos,     ...
      'p',              s.p,       ...
      'irreducible',    irr,       ...
      'periodo',        d,         ...
      'primitiva',      prim,      ...
      'rho_emp',        rho_emp,   ...
      'rho_teor',       rho_teor,  ...
      'norma2',         norma2,    ...
      'rangoFinal',     rf,        ...
      'stdFinal',       sd,        ...
      'polarizado',     pol,       ...
      'avgCentTrolls',  avgCentT,  ...
      'avgCentNormals', avgCentN,  ...
      'convTime',       tconv,     ...
      'total_trolls',   s.neg + s.pos,...
      'prop_trolls_pos', s.pos/(s.neg+s.pos) ...
    );
    cnt = cnt + 1;
end

%% 5 · EXPORTAR RESULTADOS
T = struct2table(M);
writetable(T, fullfile(resultsDir,'summary_metrics.csv'));
save(fullfile(resultsDir,'summary_metrics.mat'),'T');
disp(T);

analisis2(resultsDir);
