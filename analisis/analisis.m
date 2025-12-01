%% 8. ANÁLISIS FACTORIAL
TT = TAgg;  
% copia para análisis
% Tipos/categorización
TT.regimen = string(TT.regimen);
TT.loc     = categorical(string(TT.loc), ["low","mid","high"], 'Ordinal', true);
TT.p_num   = TT.p;   % si te resulta cómodo, alias numérico de p

% Outcome ya agregado: usa proporciones o moda
% Si quieres una etiqueta "dominante" consistente:
[~, idxMax] = max([TT.prop_outcome_consenso, TT.prop_outcome_mixto, TT.prop_outcome_polarizado], [], 2);
lab = strings(height(TT),1);
lab(idxMax==1) = "consenso";
lab(idxMax==2) = "mixto";
lab(idxMax==3) = "polarizado";
TT.outcome_mode = lab;   % no usa distCentros


%% hjjj
gscatter3(TT.lam, TT.p, TT.Ppol, TT.fracTrolls) % o usa scatter3 y colorbar
% ver si va

%% 9. antigua-----
%%catsLocAll = {'low','mid','high'};   % lista global de categorías para loc (ORDINAL)

TT.regimen = string(TT.regimen);
TT.loc     = categorical(string(TT.loc), ["low","mid","high"], "Ordinal",true);
TT.tot_trolls = TT.neg + TT.pos;
TT.fracTrolls = TT.fracTrolls;       % ya viene (0 ó 0.3)
TT.meanDist = zeros(height(TT),1);
hasD = cellfun(@(c) ~isempty(c), TT.distCentros);
TT.meanDist(hasD) = cellfun(@(c) mean(c), TT.distCentros(hasD));

% outcome categórico (por si te ayuda):
tol = 0.05;
TT.outcome = repmat("mixto", height(TT), 1);
TT.outcome(TT.rangoFinal < tol) = "consenso";
TT.outcome(TT.polarizado)       = "polarizado";

% Colapsa duplicados sin trolls (loc irrelevante):
mask0 = TT.tot_trolls==0;
if any(mask0)
    T0 = TT(mask0,:);
    [~, ia] = unique([T0.regimen, string(T0.lam)], 'rows');
    TT = [TT(~mask0,:); T0(ia,:)];
end

% "Intensidad de conectividad" como numérica (usa p directo)
TT.p_num = TT.p;                     % ya lo tienes numérico

%% 9. cpy
% Efecto conjunto por (regimen, lam, loc, frac) en métricas clave (agregadas)
G = groupsummary(TT, ["regimen","p_num","lam","loc","fracTrolls"], "mean", ...
    ["mean_convTime","mean_K_clusters","mean_rangoFinal","mean_stdFinal","prop_polarizado"]);

disp(G(:,["regimen","lam","loc","fracTrolls", ...
          "mean_mean_convTime","mean_mean_K_clusters","mean_prop_polarizado"]));



%% 10. MODELOS (efectos principales) con loc ORDINAL por cada lambda) — AGREGADA
% TODO: da error con distCentros
TT = TAgg;                                % <-- trabajar con la agregada
TT.loc = categorical(string(TT.loc), {'low','mid','high'}, 'Ordinal', true);
if ~ismember('p_num', TT.Properties.VariableNames)
    TT.p_num = TT.p;                      % alias por comodidad
end

% Usamos solo escenarios con trolls (loc sólo tiene sentido con trolls)
D_loc_all = TT(TT.fracTrolls>0 & ~isnan(TT.mean_meanDist), :);
uLam = unique(D_loc_all.lam);

CoefTabs = table(); 
ContrTabs = table();

for il = 1:numel(uLam)
    lam0 = uLam(il);
    D = D_loc_all(D_loc_all.lam==lam0, :);

    % Chequeo de mínimos: variación en p y al menos 2 niveles de loc
    if height(D) < 6 || numel(unique(D.p_num))<2 || numel(categories(D.loc))<2
        warning('Insuficientes observaciones para lam=%.2f; salto modelo.', lam0);
        continue
    end

    % Modelo estable (sin interacciones) sobre la media agregada de distancias
    mdl = fitlm(D, 'mean_meanDist ~ p_num + loc');

    % --- guarda coeficientes etiquetados
    C = mdl.Coefficients;
    C.Term = string(C.Properties.RowNames);   % pasa rownames a columna
    C.Properties.RowNames = {};
    C.lam = repmat(lam0, height(C), 1);
    CoefTabs = [CoefTabs; C]; %#ok<AGROW>

    % --- CONTRASTE: HIGH vs LOW a p = min/med/max
    pvals = [min(D.p_num), median(D.p_num), max(D.p_num)];
    catsLoc = {'low','mid','high'};
    rows = table();
    for p0 = pvals
        Xlow  = table(p0, categorical("low",  catsLoc, 'Ordinal', true), ...
                      'VariableNames', {'p_num','loc'});
        Xhigh = table(p0, categorical("high", catsLoc, 'Ordinal', true), ...
                      'VariableNames', {'p_num','loc'});
        yL = predict(mdl, Xlow);
        yH = predict(mdl, Xhigh);
        rows = [rows; table(lam0, p0, yH-yL, ...
            'VariableNames', {'lam','p_num','delta_mean_meanDist_high_minus_low'})]; %#ok<AGROW>
    end
    ContrTabs = [ContrTabs; rows]; %#ok<AGROW>

    % exporta por-lambda
    writetable(C,    fullfile(resultsDir, sprintf('AGG_mdl_meanMeanDist_coeffs_lam%.2f.csv', lam0)));
    writetable(rows, fullfile(resultsDir, sprintf('AGG_contraste_loc_high_vs_low_meanMeanDist_lam%.2f.csv', lam0)));
end

% exporta concatenados si hay
if ~isempty(CoefTabs)
    writetable(CoefTabs, fullfile(resultsDir,'AGG_mdl_meanMeanDist_coeffs_ALL.csv'));
end
if ~isempty(ContrTabs)
    writetable(ContrTabs, fullfile(resultsDir,'AGG_contraste_loc_high_vs_low_meanMeanDist_ALL.csv'));
end

%% 10.bis Modelo para tiempo de convergencia (AGREGADA)

% Opción A: usar tiempo CONDICIONAL (donde hubo alguna convergencia)
if ismember('mean_convTime_cond', TT.Properties.VariableNames)
    D_ct = TT(TT.fracTrolls>0 & ~isnan(TT.mean_convTime_cond), :);
    if ~isempty(D_ct) && numel(unique(D_ct.p_num))>=2 && numel(categories(D_ct.loc))>=2
        D_ct.loc = categorical(string(D_ct.loc), {'low','mid','high'}, 'Ordinal', true);
        mdl_time_cond = fitlm(D_ct, 'mean_convTime_cond ~ p_num + fracTrolls + loc');
        writetable(mdl_time_cond.Coefficients, fullfile(resultsDir,'AGG_mdl_convTimeCOND_coeffs.csv'));
    else
        warning('Datos insuficientes para mean\\_convTime\\_cond: no ajusto mdl_time_cond.');
    end
end

% Opción B: usar tiempo IMPUTADO (si lo calculaste; útil para tener siempre respuesta)
if ismember('mean_convTime_imputed', TT.Properties.VariableNames)
    D_cti = TT(TT.fracTrolls>0 & ~isnan(TT.mean_convTime_imputed), :);
    if ~isempty(D_cti) && numel(unique(D_cti.p_num))>=2 && numel(categories(D_cti.loc))>=2
        D_cti.loc = categorical(string(D_cti.loc), {'low','mid','high'}, 'Ordinal', true);
        mdl_time_imp = fitlm(D_cti, 'mean_convTime_imputed ~ p_num + fracTrolls + loc');
        writetable(mdl_time_imp.Coefficients, fullfile(resultsDir,'AGG_mdl_convTimeIMP_coeffs.csv'));
    else
        warning('Datos insuficientes para mean\\_convTime\\_imputed: no ajusto mdl_time_imp.');
    end
end

%% 10bis. MODELO GLOBAL para meanDist (para usar en CONTRASTES)
D_all = TT(TT.tot_trolls>0, :);                          % loc solo con trolls
D_all.loc = categorical(string(D_all.loc), {'low','mid','high'}, 'Ordinal', true);

if height(D_all) >= 6 && numel(unique(D_all.p_num))>=2 && numel(categories(D_all.loc))>=2
    % Modelo sencillo y estable (sin interacciones para evitar rank deficiency)
    mdl_dist = fitlm(D_all, 'meanDist ~ p_num + fracTrolls + loc');
else
    warning('Datos insuficientes para meanDist global: no ajusto mdl_dist.');
    mdl_dist = [];   % marcador para saltar contrastes
end
%% 11. CONTRASTES (solo si existe mdl_dist) — AGREGADA
if exist('mdl_dist','var') && ~isempty(mdl_dist)
    catsLocAll = {'low','mid','high'};   % categorías ordenadas

    % Rejilla de p (usa la del propio TT)
    p_grid = [min(TT.p_num), median(TT.p_num), max(TT.p_num)];
    % Fracciones de trolls a contrastar
    f_grid = [0.3];   % puedes añadir 0.15 si quieres

    % 11.1 HIGH vs LOW a p y frac fijos (predice mean_meanDist)
    Contr = table(); k = 1;
    for p0 = p_grid
        for f0 = f_grid
            Xlow  = table(p0, f0, categorical("low",  catsLocAll, 'Ordinal', true), ...
                          'VariableNames', {'p_num','fracTrolls','loc'});
            Xhigh = table(p0, f0, categorical("high", catsLocAll, 'Ordinal', true), ...
                          'VariableNames', {'p_num','fracTrolls','loc'});

            yL = predict(mdl_dist, Xlow);    % predicción de mean_meanDist
            yH = predict(mdl_dist, Xhigh);

            Contr.k(k,1)             = k;
            Contr.p_num(k,1)         = p0;
            Contr.fracTrolls(k,1)    = f0;
            Contr.delta_meanMeanDist(k,1) = yH - yL;   % >0: high separa más que low
            k = k + 1;
        end
    end
    writetable(Contr, fullfile(resultsDir,'AGG_contraste_loc_high_vs_low_meanMeanDist.csv'));

    % 11.2 EFECTO de subir p (desconectada→umbral→fuerte) con loc=mid, frac=0.3
    Xp = table(p_grid', repmat(0.3,numel(p_grid),1), ...
               categorical(repmat("mid",numel(p_grid),1), catsLocAll, 'Ordinal', true), ...
               'VariableNames', {'p_num','fracTrolls','loc'});
    yp = predict(mdl_dist, Xp);
    Eff_p = table(p_grid', yp, 'VariableNames', {'p_num','pred_meanMeanDist'});
    writetable(Eff_p, fullfile(resultsDir,'AGG_efecto_p_meanMeanDist.csv'));

    % 11.3 EFECTO de aumentar fracTrolls (0 → 0.3) con loc=mid
    Xf = table(repmat(p_grid',2,1), [zeros(numel(p_grid),1); 0.3*ones(numel(p_grid),1)], ...
               categorical(repmat("mid",2*numel(p_grid),1), catsLocAll, 'Ordinal', true), ...
               'VariableNames', {'p_num','fracTrolls','loc'});
    yf = predict(mdl_dist, Xf);
    Ef = table(Xf.p_num, Xf.fracTrolls, yf, 'VariableNames', ...
               {'p_num','fracTrolls','pred_meanMeanDist'});
    writetable(Ef, fullfile(resultsDir,'AGG_efecto_fracTrolls_meanMeanDist.csv'));
else
    warning('Sin mdl_dist: se omiten los contrastes de mean\\_meanDist.');
end



%% 11. CONTRASTES 
%11.1. HIGH VS LOW ( A P Y FRACTROLLS FIJADOS)
% Grid de predicción para comparar loc a p y frac fijos
p_grid   = [min(TT.p_num), median(TT.p_num), max(TT.p_num)];
f_grid   = [0.15 0.3];  % con trolls
lam_grid = unique(TT.lam)';  % por si quieres estratificar por lambda

Contr = table(); k=1;
for p0 = p_grid
  for f0 = f_grid
    for lam0 = lam_grid
      XcatsLocAll = {'low','mid','high'};
      Xlow  = table(p0, f0, categorical("low",  catsLocAll, 'Ordinal', true), lam0, ...
              'VariableNames', {'p_num','fracTrolls','loc','lam'});
      Xhigh = table(p0, f0, categorical("high", catsLocAll, 'Ordinal', true), lam0, ...
              'VariableNames', {'p_num','fracTrolls','loc','lam'});

      yL = predict(mdl_dist, Xlow);
      yH = predict(mdl_dist, Xhigh);
      Contr.k(k,1)      = k;
      Contr.p_num(k,1)  = p0;
      Contr.frac(k,1)   = f0;
      Contr.lam(k,1)    = lam0;
      Contr.deltaDist(k,1) = yH - yL;  % >0: high separa más que low
      k = k+1;
    end
  end
end
writetable(Contr, fullfile(resultsDir,'contraste_loc_high_vs_low_meanDist.csv'));

%11.2. EFECTO DE INCREMENTAR P (DESCONECTADA -> FUERTE)
% Mantén loc=mid y frac=0.3, observa cómo cambia la predicción al subir p
Xp = table(p_grid', repmat(0.3,numel(p_grid),1), ...
           categorical(repmat("mid",numel(p_grid),1), catsLocAll, 'Ordinal', true), ...
           repmat(median(TT.lam),numel(p_grid),1), ...
           'VariableNames', {'p_num','fracTrolls','loc','lam'});

yp = predict(mdl_dist, Xp);
Eff_p = table(p_grid', yp, 'VariableNames', {'p_num','pred_meanDist'});
writetable(Eff_p, fullfile(resultsDir,'efecto_p_meanDist.csv'));

% EFECTO DE INCREMENTAR PROPTROLLS
% Con p en cada régimen (p_grid) y loc=mid
Xf = table(repmat(p_grid',2,1), [zeros(numel(p_grid),1); 0.3*ones(numel(p_grid),1)], ...
           categorical(repmat("mid",2*numel(p_grid),1), catsLocAll, 'Ordinal', true), ...
           repmat(median(TT.lam),2*numel(p_grid),1), ...
           'VariableNames', {'p_num','fracTrolls','loc','lam'});
yf = predict(mdl_dist, Xf);
Ef = table(Xf.p_num, Xf.fracTrolls, yf, 'VariableNames', {'p_num','fracTrolls','pred_meanDist'});
writetable(Ef, fullfile(resultsDir,'efecto_fracTrolls_meanDist.csv'));

%GRAFICOS 
% 4.1. MeanDist vs lambda (strat por loc) con trolls
Tc = TT(TT.tot_trolls>0, :);
figure; hold on;
uLoc = categories(Tc.loc);
for i=1:numel(uLoc)
    m = Tc.loc==uLoc{i};
    plot(Tc.lam(m), Tc.meanDist(m), '.', 'DisplayName', uLoc{i});
end
xlabel('\lambda (no tercos)'); ylabel('Distancia media entre centros');
legend('Location','best'); title('Separación por ubicación de trolls');

% 4.2. K_clusters por régimen (histos)
figure; tiledlayout(1,3);
uReg = unique(Tc.regimen);
for i=1:numel(uReg)
    nexttile; histogram(Tc.K_clusters(Tc.regimen==uReg(i)));
    title("K por " + uReg(i)); xlabel('K'); ylabel('freq');
end

% 4.3. ConvTime vs p (p_num) con facet por loc (marcas por lambda)
figure; hold on;
scatter(Tc.p_num, Tc.convTime, 10, Tc.lam, 'filled'); % color por lam
xlabel('p'); ylabel('Tiempo de convergencia'); title('Convergencia vs p (color=\lambda)');
colormap turbo; colorbar;

writetable(G,    fullfile(resultsDir,'agg_regimen_p_lam_loc_frac.csv'));
writetable(Gc,   fullfile(resultsDir,'agg_loc_con_trolls.csv'));
%% RESULTADOS OBTENIDOS TRAS ANÁLISIS

% Llamada a tu análisis posterior (si aplica)
%if exist('analisis2','file')==2
%    analisis2(resultsDir);
%end

% === Tabla 1: convTime cuando NO hay trolls (fracTrolls=0)
T0 = TT(TT.neg+TT.pos==0,:);
T_tab1 = groupsummary(T0, ["regimen","lam"], "median", "convTime");
T_tab1 = removevars(T_tab1, setdiff(T_tab1.Properties.VariableNames, ...
           ["regimen","lam","median_convTime"]));
writetable(T_tab1, fullfile(resultsDir,'TABLA1_convTime_sin_trolls.csv'));
disp('TABLA1 – Convergencia sin trolls (mediana convTime):');
disp(T_tab1);

% === Tabla 2: con trolls (0.15 y 0.30) – clústeres y separaciones
Tc = TT(TT.neg+TT.pos>0,:);
% meanDist = media de distCentros (debe ser >=0 si clusters1D devuelve distancias correctas)
T_tab2 = groupsummary(Tc, ["regimen","fracTrolls","loc"], "mean", ["K_clusters","meanDist"]);
T_tab2 = sortrows(T_tab2, ["regimen","fracTrolls","loc"]);
writetable(T_tab2, fullfile(resultsDir,'TABLA2_fragmentacion_con_trolls.csv'));
disp('TABLA2 – Fragmentación con trolls (media K y meanDist):');
disp(T_tab2(:,["regimen","fracTrolls","loc","mean_K_clusters","mean_meanDist"]));
% Mapea regimen a p teórico medio para el eje x
regKey = groupsummary(TT, "regimen", "mean", "p");
reg2p  = containers.Map(regKey.regimen, regKey.mean_p);

T0 = TT(TT.neg+TT.pos==0,:);
x  = arrayfun(@(r) reg2p(char(r)), T0.regimen);  % p en eje x
figure; hold on;
ms = {'o','s','^'}; uLam = unique(T0.lam);
for i=1:numel(uLam)
    m = T0.lam==uLam(i);
    scatter(x(m), T0.convTime(m), 45, 'filled', ms{i});
end
xlabel('p (intensidad ER)'); ylabel('Tiempo de convergencia (iter.)');
legend(cellstr("λ="+string(uLam)),'Location','northeastoutside');
title('Sin trolls: p reduce el tiempo de convergencia');
grid on; box on;
saveas(gcf, fullfile(resultsDir,'FIG_A_convTime_sin_trolls.png'));
print(gcf, fullfile(resultsDir,'FIG_A_convTime_sin_trolls.pdf'), '-dpdf', '-bestfit');
Tc = TT(TT.neg+TT.pos>0,:);
Tc.loc = categorical(string(Tc.loc), ["low","mid","high"], "Ordinal",true);
% Media de K por (regimen, loc) colapsando fracTrolls (si quieres, separa por frac)
Kgrid = groupsummary(Tc, ["regimen","loc"], "mean", "K_clusters");
% Pivot
ureg = unique(Kgrid.regimen);
uloc = categories(Tc.loc);
Mat  = nan(numel(ureg), numel(uloc));
for i=1:numel(ureg)
  for j=1:numel(uloc)
    m = Kgrid.regimen==ureg(i) & string(Kgrid.loc)==uloc{j};
    if any(m), Mat(i,j) = Kgrid.mean_K_clusters(m); end
  end
end
figure;
imagesc(Mat); colorbar; caxis([min(Mat,[],'all') max(Mat,[],'all')]);
xticks(1:numel(uloc)); xticklabels(uloc);
yticks(1:numel(ureg)); yticklabels(ureg);
xlabel('Ubicación de trolls (PageRank)'); ylabel('Régimen (p)');
title('Media de K (nº de clústeres) con trolls');
for i=1:numel(ureg)
  for j=1:numel(uloc)
    if ~isnan(Mat(i,j))
      text(j,i, sprintf('%.1f',Mat(i,j)), 'HorizontalAlignment','center','Color','w','FontWeight','bold');
    end
  end
end
set(gca,'YDir','normal'); box on;
saveas(gcf, fullfile(resultsDir,'FIG_B_Kclusters_heatmap.png'));
print(gcf, fullfile(resultsDir,'FIG_B_Kclusters_heatmap.pdf'), '-dpdf', '-bestfit');
Tc = TT(TT.neg+TT.pos>0 & TT.fracTrolls==0.30,:);
Tc.loc = categorical(string(Tc.loc), ["low","mid","high"], "Ordinal",true);
G = groupsummary(Tc, ["regimen","loc"], "mean", "meanDist");

figure; tiledlayout(1, numel(unique(G.regimen)));
ureg = unique(G.regimen);
for i=1:numel(ureg)
    nexttile; Gi = G(G.regimen==ureg(i),:);
    bar(categorical(string(Gi.loc), ["low","mid","high"]), Gi.mean_meanDist);
    title("Régimen: " + string(ureg(i)));
    ylabel('Distancia media entre centros');
    xlabel('Ubicación de trolls');
    box on; grid on;
end
sgtitle('Con 30% trolls: separación por ubicación (PageRank)');
saveas(gcf, fullfile(resultsDir,'FIG_C_meanDist_bar_30pct.png'));
print(gcf, fullfile(resultsDir,'FIG_C_meanDist_bar_30pct.pdf'), '-dpdf', '-bestfit');
Td = TT(TT.regimen=="desconectada" & TT.neg+TT.pos>0,:);
Td.loc = categorical(string(Td.loc), ["low","mid","high"], "Ordinal",true);
C = groupsummary(Td, ["fracTrolls","loc"], "mean", "meanDist");
disp('CONTRASTE simple (desconectada): meanDist por loc y fracción de trolls');
disp(C(:,["fracTrolls","loc","mean_meanDist"]));
writetable(C, fullfile(resultsDir,'CONTRASTE_desconectada_loc_frac_meanDist.csv'));


%% A) Sensibilidad tipo Sobol (eta^2 por factor) ============================
% Usamos varianza explicada por cada factor (main effects) sobre tu malla
% de escenarios. Resulta muy interpretable y sin ajustar modelos.

OUTS = {'K_clusters','meanDist','convTime'};
FACT = {'regimen','loc','fracTrolls','lam'};   % p_num lo representas con regimen; si quieres, añade 'p_num'

SobolTable = table();
for io = 1:numel(OUTS)
    y = TT.(OUTS{io});
    good = ~isnan(y);           % por si convTime tiene NaN
    y  = y(good);
    VY = var(y,1);              % Var total (poblacional)
    if VY < eps, continue; end

    S = table('Size',[0 3], 'VariableTypes',{'string','string','double'}, ...
              'VariableNames',{'output','factor','S1'});
    for jf = 1:numel(FACT)
        f = FACT{jf};
        g  = TT.(f)(good);
        % media por nivel del factor y varianza de esa media (E[Var] / Var total)
        [lev,~,idx] = unique(string(g));    % robusto para num/cat
        m = accumarray(idx, y, [], @mean);
        w = accumarray(idx, 1, [], @sum);
        % Var entre-niveles ponderada:
        mu = sum(w.*m)/sum(w);
        V_between = sum(w.*(m - mu).^2) / sum(w);
        S1 = V_between / VY;
        S = [S; {OUTS{io}, f, S1}]; %#ok<AGROW>
    end
    SobolTable = [SobolTable; S]; %#ok<AGROW>
end
disp('--- Sobol-like (S1 por factor y output) ---');
disp(SobolTable);

% Barplots por output
uOut = unique(SobolTable.output);
for i=1:numel(uOut)
    sub = SobolTable(SobolTable.output==uOut(i),:);
    figure; bar(categorical(sub.factor), sub.S1);
    ylim([0 1]); ylabel('S1 (varianza explicada)'); title("Sobol-like: " + uOut(i));
end

%% B) Morris post-hoc (efectos elementales observacionales)
% Objetivo: estimar importancia de cada factor a partir de diferencias
% locales entre combinaciones observadas (vecinos en niveles adyacentes).

% ---- 1) Definir niveles con orden fijo para factores categóricos
ord_regimen = ["desconectada","umbral","fuerte"];
ord_loc     = ["low","mid","high"];

% Limpiar/forzar tipos
TTm = TT;
TTm.regimen = categorical(string(TTm.regimen), ord_regimen, 'Ordinal', true);
TTm.loc     = categorical(string(TTm.loc),     ord_loc,     'Ordinal', true);

% Outputs a evaluar (usa sólo los que existan)
OUTS_all = {'K_clusters','meanDist','convTime'};
OUTS = OUTS_all(ismember(OUTS_all, TTm.Properties.VariableNames));

% Tabla de resultados final (cabecera fija para evitar errores de concatenación)
Morris = table('Size',[0 4], ...
    'VariableNames', {'output','factor','mu_star','sigma'}, ...
    'VariableTypes', {'string','string','double','double'});

for io = 1:numel(OUTS)
    yname = OUTS{io};

    % ---- 2) Agregar: media por combinación de factores (evita duplicados)
    G = groupsummary(TTm, {'regimen','loc','fracTrolls','lam'}, 'mean', yname);
    % El nombre sale como 'mean_<yname>'; renómbralo a 'ymean'
    colMean = "mean_" + string(yname);
    G.ymean = G.(colMean);
    G.(colMean) = [];  % ya no lo necesitamos

    % Niveles presentes (ordenados)
    lev.regimen    = categories(TTm.regimen);
    lev.loc        = categories(TTm.loc);
    lev.fracTrolls = unique(G.fracTrolls);  lev.fracTrolls = sort(lev.fracTrolls);
    lev.lam        = unique(G.lam);         lev.lam        = sort(lev.lam);

    % Almacenes de efectos elementales
    EE_vals = [];  EE_lab = {};

    % ---- 3) Recorre cada fila del grid agregado y busca vecinos por factor
    for i = 1:height(G)
        r  = G.regimen(i);  l  = G.loc(i);
        f  = G.fracTrolls(i);
        la = G.lam(i);
        yi = G.ymean(i);
        if isnan(yi), continue; end

        % --- 3a) Vecinos en 'regimen' (nivel previo/siguiente)
        kR = find(lev.regimen == string(r));
        for step = [-1 1]
            kN = kR + step;
            if kN>=1 && kN<=numel(lev.regimen)
                r2 = categorical(lev.regimen{kN}, lev.regimen, 'Ordinal', true);
                j = G.regimen==r2 & G.loc==l & G.fracTrolls==f & G.lam==la;
                if any(j)
                    yj = G.ymean(j);
                    if ~isnan(yj)
                        EE_vals(end+1) = (yj - yi); %#ok<AGROW>
                        EE_lab{end+1}  = 'regimen'; %#ok<AGROW>
                    end
                end
            end
        end

        % --- 3b) Vecinos en 'loc'
        kL = find(lev.loc == string(l));
        for step = [-1 1]
            kN = kL + step;
            if kN>=1 && kN<=numel(lev.loc)
                l2 = categorical(lev.loc{kN}, lev.loc, 'Ordinal', true);
                j = G.regimen==r & G.loc==l2 & G.fracTrolls==f & G.lam==la;
                if any(j)
                    yj = G.ymean(j);
                    if ~isnan(yj)
                        EE_vals(end+1) = (yj - yi); %#ok<AGROW>
                        EE_lab{end+1}  = 'loc'; %#ok<AGROW>
                    end
                end
            end
        end

        % --- 3c) Vecinos en 'fracTrolls'
        kF = find(abs(lev.fracTrolls - f) < 1e-12);
        for step = [-1 1]
            kN = kF + step;
            if kN>=1 && kN<=numel(lev.fracTrolls)
                f2 = lev.fracTrolls(kN);
                j = G.regimen==r & G.loc==l & G.fracTrolls==f2 & abs(G.lam - la)<1e-12;
                if any(j)
                    yj = G.ymean(j);
                    if ~isnan(yj)
                        EE_vals(end+1) = (yj - yi); %#ok<AGROW>
                        EE_lab{end+1}  = 'fracTrolls'; %#ok<AGROW>
                    end
                end
            end
        end

        % --- 3d) Vecinos en 'lam'
        kA = find(abs(lev.lam - la) < 1e-12);
        for step = [-1 1]
            kN = kA + step;
            if kN>=1 && kN<=numel(lev.lam)
                la2 = lev.lam(kN);
                j = G.regimen==r & G.loc==l & abs(G.fracTrolls - f)<1e-12 & abs(G.lam - la2)<1e-12;
                if any(j)
                    yj = G.ymean(j);
                    if ~isnan(yj)
                        EE_vals(end+1) = (yj - yi); %#ok<AGROW>
                        EE_lab{end+1}  = 'lam'; %#ok<AGROW>
                    end
                end
            end
        end
    end

    % ---- 4) Resumen por factor: mu* = mean(|EE|), sigma = std(EE)
    if ~isempty(EE_vals)
        Ttmp = table(categorical(EE_lab(:)), abs(EE_vals(:)), EE_vals(:), ...
                     'VariableNames', {'factor','absEE','EE'});
        Tm = groupsummary(Ttmp, 'factor', {'mean','std'}, {'absEE','EE'});
        % Renombra y añade etiqueta de output
        Tm.Properties.VariableNames(end-1:end) = {'mu_star','sigma'};
        Tm.output = repmat(string(yname), height(Tm), 1);
        % Selecciona columnas finales en orden
        Tm = Tm(:, {'output','factor','mu_star','sigma'});

        Morris = [Morris; Tm]; %#ok<AGROW>
    end
end

disp('--- Morris post-hoc (mu* y sigma) ---');
disp(Morris);

% ---- 5) Visualización: diagrama mu* vs sigma por output
uOut = unique(Morris.output);
for i = 1:numel(uOut)
    sub = Morris(Morris.output==uOut(i),:);
    figure; scatter(sub.mu_star, sub.sigma, 60, 'filled'); grid on;
    text(sub.mu_star, sub.sigma, string(sub.factor), ...
        'VerticalAlignment','bottom','HorizontalAlignment','left');
    xlabel('\mu^* (tamaño de efecto medio |EE|)');
    ylabel('\sigma (variabilidad → no linealidad/interacción)');
    title("Morris post-hoc: " + uOut(i));
end


%% C) Correlaciones parciales y MI (heatmap) ================================
% Matriz de inputs numéricos (codificamos factores)
X = table();
X.p_num      = TT.p;                 % ya numérico
X.fracTrolls = TT.fracTrolls;
% codificación ordinal/entera para categoricals:
X.loc_num     = double(grp2idx(categorical(TT.loc, {'low','mid','high'})));
X.regimen_num = double(grp2idx(categorical(TT.regimen)));

Y = table();
Y.K_clusters = TT.K_clusters;
Y.meanDist   = TT.meanDist;
Y.convTime   = TT.convTime;

% Correlaciones de Spearman (robusto) entre cada input y cada output
inNames  = X.Properties.VariableNames;
outNames = Y.Properties.VariableNames;

Rho = nan(numel(outNames), numel(inNames));
for i=1:numel(outNames)
    for j=1:numel(inNames)
        v1 = Y.(outNames{i});
        v2 = X.(inNames{j});
        m  = ~isnan(v1) & ~isnan(v2);
        if nnz(m)>3
            Rho(i,j) = corr(v1(m), v2(m), 'Type','Spearman');
        end
    end
end

% Información mutua (discretizando a 10 bins)
nb=10;
MI = nan(numel(outNames), numel(inNames));
for i=1:numel(outNames)
    for j=1:numel(inNames)
        v1 = Y.(outNames{i}); v2 = X.(inNames{j});
        m = ~isnan(v1) & ~isnan(v2);
        if nnz(m)>10
            MI(i,j) = mutual_info_disc(v1(m), v2(m), nb);
        end
    end
end

% Heatmaps
figure; imagesc(Rho); colorbar; colormap parula;
set(gca,'XTick',1:numel(inNames),'XTickLabel',inNames, ...
        'YTick',1:numel(outNames),'YTickLabel',outNames, 'XTickLabelRotation',45);
title('Spearman (inputs vs outputs)');

figure; imagesc(MI); colorbar; colormap parula;
set(gca,'XTick',1:numel(inNames),'XTickLabel',inNames, ...
        'YTick',1:numel(outNames),'YTickLabel',outNames, 'XTickLabelRotation',45);
title('Información mutua (discretizada)');
function mi = mutual_info_disc(x, y, nb)
    % Discretiza y calcula I(X;Y) con estimador plug-in
    x = x(:); y = y(:);
    [~,~,xb] = histcounts(x, nb);
    [~,~,yb] = histcounts(y, nb);
    valid = xb>0 & yb>0;
    xb = xb(valid); yb = yb(valid);
    if isempty(xb), mi = NaN; return; end
    N = numel(xb);
    Pxy = accumarray([xb,yb], 1, [nb,nb]) / N;
    Px  = sum(Pxy,2);
    Py  = sum(Pxy,1);
    [ii,jj,v] = find(Pxy);
    mi = 0;
    for k=1:numel(v)
        mi = mi + v(k) * log( v(k) / (Px(ii(k))*Py(jj(k))) );
    end
    mi = mi / log(2);  % bits
end

%% D) Efectos marginales ====================================================
% 1) meanDist vs p (por loc) con y sin trolls
figure; hold on;
for L = {'low','mid','high'}
    m = TT.loc==string(L);
    scatter(TT.p(m), TT.meanDist(m), 18, 'filled', 'DisplayName', string(L));
end
xlabel('p'); ylabel('meanDist'); legend('Location','best');
title('Efecto marginal de p sobre meanDist (coloreado por loc)');

% 2) K_clusters vs fracTrolls (por regimen)
figure; hold on;
for R = unique(string(TT.regimen))'
    m = string(TT.regimen)==R;
    scatter(TT.fracTrolls(m), TT.K_clusters(m), 18, 'filled', 'DisplayName', R);
end
xlabel('fracTrolls'); ylabel('K\_clusters'); legend('Location','best');
title('Efecto marginal de fracTrolls sobre K');

% 3) convTime vs p por regimen (solo filas con convTime válido)
C = TT(~isnan(TT.convTime),:);
if ~isempty(C)
    figure; hold on;
    for R = unique(string(C.regimen))'
        m = string(C.regimen)==R;
        scatter(C.p(m), C.convTime(m), 25, 'filled', 'DisplayName', R);
    end
    xlabel('p'); ylabel('convTime'); legend('Location','best');
    title('Efecto marginal de p sobre convTime (por régimen)');
end
%%%%


