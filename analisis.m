function analisis(resultsDir)
% analisis.m – Análisis centrado en trolls con contraste de hipótesis y p-valores

if nargin<1
    resultsDir = uigetdir('C:/Users/HP/Documents/4 mates/TFG/modelomatlab', ...
                         'Selecciona la carpeta de resultados generada por main.m');
    if isequal(resultsDir,0)
        error('No seleccionaste ninguna carpeta.');
    end
end

% 1 · Carga la tabla
T = readtable(fullfile(resultsDir,'summary_metrics.csv'));

%% 2 · Variables derivadas
T.total_trolls = T.neg + T.pos;
T.prop_trolls_pos = T.pos ./ (T.total_trolls + eps);
T.prop_trolls_neg = T.neg ./ (T.total_trolls + eps);

%% 2.1 · Cargar métricas adicionales si existen
piFile = fullfile(resultsDir,'piTx0.csv');
if isfile(piFile)
    T_pi = readtable(piFile);
    T = join(T, T_pi, 'Keys','tag');
    for nval = unique(T.n)'
        Tn = T(T.n==nval,:);
        simVals = unique(Tn.sim);
        pVals   = unique(Tn.p);

        for simval = simVals'
            for pval = pVals'
                subset = Tn(Tn.sim==simval & abs(Tn.p-pval)<1e-6,:);
                if isempty(subset), continue; end

                % Malla de proporciones
                [NEG, POS] = meshgrid(0:0.1:1, 0:0.1:1);
                ZZ_pi    = nan(size(NEG));
                ZZ_cons  = nan(size(NEG));
                ZZ_ctime = nan(size(NEG));
                
                for i = 1:size(NEG,1)
                    for j = 1:size(NEG,2)
                        mask = abs(subset.prop_trolls_neg - NEG(i,j)) < 0.05 & ...
                               abs(subset.prop_trolls_pos - POS(i,j)) < 0.05;
                        if any(mask)
                            ZZ_pi(i,j)    = mean(subset.piTx0(mask), 'omitnan');
                            ZZ_cons(i,j)  = mean(~isnan(subset.convTime(mask)));   % P(consenso)
                            ZZ_ctime(i,j) = mean(subset.convTime(mask),  'omitnan');
                        end
                    end
                end

                ttlSim = tern(simval,'sim','dir');
                tag    = sprintf('n%d_p%.2f_%s', nval, pval, ttlSim);

                % ===== 
                % Superficie pi^T x0
                fig = figure('Visible','off');
                surf(NEG, POS, ZZ_pi, 'EdgeColor','none');
                shading interp; colormap turbo; view(45,30);
                xlabel('Prop trolls neg'); ylabel('Prop trolls pos');
                zlabel('\pi^T x_0');
                title(['\pi^T x_0 — ' tag]);
                colorbar;
                saveas(fig, fullfile(resultsDir, ['piTx0_3D_' tag '.png']));
                close(fig);

                % Probabilidad de consenso
                fig = figure('Visible','off');
                surf(NEG, POS, ZZ_cons, 'EdgeColor','none');
                shading interp; colormap parula; view(45,30);
                xlabel('Prop trolls neg'); ylabel('Prop trolls pos');
                zlabel('P(consenso)');
                title(['Prob Consenso — ' tag]);
                colorbar;
                saveas(fig, fullfile(resultsDir, ['consenso3D_' tag '.png']));
                close(fig);

                % Tiempo medio de convergencia
                fig = figure('Visible','off');
                surf(NEG, POS, ZZ_ctime, 'EdgeColor','none');
                shading interp; colormap hot; view(45,30);
                xlabel('Prop trolls neg'); ylabel('Prop trolls pos');
                zlabel('Tiempo conv');
                title(['ConvTime — ' tag]);
                colorbar;
                saveas(fig, fullfile(resultsDir, ['convTime3D_' tag '.png']));
                close(fig);
            end
        end
    end
end


%% 3 · Correlaciones globales
R_ct_tt = corrcoef(T.convTime, T.total_trolls, 'Rows','complete');
fprintf('\nCorrelación convTime vs. total_trolls: %.3f\n', R_ct_tt(1,2));
R_ct_prop = corrcoef(T.convTime, T.prop_trolls_pos, 'Rows','complete');
fprintf('Correlación convTime vs. prop_trolls_pos: %.3f\n', R_ct_prop(1,2));

%% 4 · Regresiones globales (logística y lineal)
Xlog = [T.total_trolls, T.prop_trolls_pos, T.total_trolls .* T.prop_trolls_pos];
Ylog = double(T.polarizado);
[b_log, ~, stats_log] = glmfit(Xlog, Ylog, 'binomial','link','logit');
Tlog = table({'Intercept';'total_trolls';'prop_pos';'interacción'}, b_log, stats_log.p, ...
      'VariableNames', {'Predictor', 'Estimate', 'pValue'});
fprintf('\nRegresión logística (logit):\n');
disp(Tlog);

mask = ~isnan(T.convTime);
Xlin = [T.total_trolls(mask), T.prop_trolls_pos(mask), T.total_trolls(mask) .* T.prop_trolls_pos(mask)];
tbl_lin = table(Xlin(:,1), Xlin(:,2), Xlin(:,3), T.convTime(mask), ...
    'VariableNames', {'total_trolls','prop_pos','interaccion','convTime'});
modelo_lin = fitlm(tbl_lin, 'convTime ~ total_trolls + prop_pos + interaccion');
fprintf('\nRegresión lineal convTime ~ total_trolls + prop_pos + interacción:\n');
disp(modelo_lin);

%% 5 · Análisis por tamaño de red (n) con gráficos de coeficientes
sizes = unique(T.n);
colors = lines(numel(sizes));
coefs_log = zeros(4, numel(sizes));
pvals_log = zeros(4, numel(sizes));
coefs_lin = zeros(4, numel(sizes));
pvals_lin = zeros(4, numel(sizes));

resumen_n = table();

for i = 1:numel(sizes)
    fprintf('\n==== Análisis para n = %d ====\n', sizes(i));
    Tn = T(T.n == sizes(i), :);

    % Logística
    Xlog_n = [Tn.total_trolls, Tn.prop_trolls_pos, Tn.total_trolls .* Tn.prop_trolls_pos];
    Ylog_n = double(Tn.polarizado);
    [b_log_n, ~, stats_log_n] = glmfit(Xlog_n, Ylog_n, 'binomial','link','logit');
    coefs_log(:,i) = b_log_n;
    pvals_log(:,i) = stats_log_n.p;

    % Lineal
    mask_n = ~isnan(Tn.convTime);
    Xlin_n = [Tn.total_trolls(mask_n), Tn.prop_trolls_pos(mask_n), Tn.total_trolls(mask_n).*Tn.prop_trolls_pos(mask_n)];
    tbl_lin_n = table(Xlin_n(:,1), Xlin_n(:,2), Xlin_n(:,3), Tn.convTime(mask_n), ...
        'VariableNames', {'total_trolls','prop_pos','interaccion','convTime'});
    modelo_lin_n = fitlm(tbl_lin_n, 'convTime ~ total_trolls + prop_pos + interaccion');
    coefs_lin(:,i) = modelo_lin_n.Coefficients.Estimate;
    pvals_lin(:,i) = modelo_lin_n.Coefficients.pValue;

    % === Gráfico 3D de polarización por total_trolls y prop_pos ===
    tt_vals = unique(Tn.total_trolls);
    pp_vals = linspace(0,1,11);
    [TT, PP] = meshgrid(tt_vals, pp_vals);
    ZZ_pol = nan(size(TT));
    ZZ_ct  = nan(size(TT));
    for r = 1:size(TT,1)
        for c = 1:size(TT,2)
            mask3D = Tn.total_trolls == TT(r,c) & abs(Tn.prop_trolls_pos - PP(r,c)) < 0.05;
            ZZ_pol(r,c) = mean(Tn.polarizado(mask3D), 'omitnan');
            ZZ_ct(r,c)  = mean(Tn.convTime(mask3D), 'omitnan');
        end
    end
    fig3d = figure('Visible','off');
    surf(TT, PP, ZZ_pol);
    xlabel('Total trolls'); ylabel('Proporción trolls positivos'); zlabel('% polarización');
    title(sprintf('Polarización esperada para n = %d', sizes(i)));
    saveas(fig3d, fullfile(resultsDir, sprintf('polarizacion3D_n%d.png', sizes(i))));
    close(fig3d);

    fig3d_ct = figure('Visible','off');
    surf(TT, PP, ZZ_ct);
    xlabel('Total trolls'); ylabel('Proporción trolls positivos'); zlabel('Tiempo medio de convergencia');
    title(sprintf('Convergencia esperada para n = %d', sizes(i)));
    saveas(fig3d_ct, fullfile(resultsDir, sprintf('convergencia3D_n%d.png', sizes(i))));
    close(fig3d_ct);

    % === Tabla resumen
    resumen_n = [resumen_n; table(sizes(i), mean(Tn.polarizado), mean(Tn.convTime,'omitnan'), ...
                   mean(Tn.total_trolls), mean(Tn.prop_trolls_pos), mean(Tn.prop_trolls_neg), ...
                   'VariableNames', {'n','pct_polarizado','media_convTime','media_totalTrolls', ...
                                      'media_prop_pos','media_prop_neg'})];
end

% Guardar resumen
writetable(resumen_n, fullfile(resultsDir,'resumen_por_n.csv'));

% Gráficos: regresión logística
fig = figure; clf;
bar(coefs_log');
hold on;
legend({'Intercept','total\_trolls','prop\_pos','interacción'}, 'Location','best');
title('Coeficientes Regresión Logística por tamaño de red');
xlabel('Tamaño de red (n)');
xlim([0.5 numel(sizes)+0.5]);
set(gca,'XTickLabel',string(sizes));
ylabel('Estimación');
saveas(fig, fullfile(resultsDir,'coef_logistic_por_n.png'));

fig = figure; clf;
bar(pvals_log');
legend({'Intercept','total\_trolls','prop\_pos','interacción'}, 'Location','best');
title('p-valores Regresión Logística por tamaño de red');
xlabel('Tamaño de red (n)');
xlim([0.5 numel(sizes)+0.5]);
set(gca,'XTickLabel',string(sizes));
ylabel('p-valor');
saveas(fig, fullfile(resultsDir,'pvals_logistic_por_n.png'));

% Gráficos: regresión lineal
fig = figure; clf;
bar(coefs_lin');
hold on;
legend({'Intercept','total\_trolls','prop\_pos','interacción'}, 'Location','best');
title('Coeficientes Regresión Lineal por tamaño de red');
xlabel('Tamaño de red (n)');
xlim([0.5 numel(sizes)+0.5]);
set(gca,'XTickLabel',string(sizes));
ylabel('Estimación');
saveas(fig, fullfile(resultsDir,'coef_lineal_por_n.png'));

fig = figure; clf;
bar(pvals_lin');
legend({'Intercept','total\_trolls','prop\_pos','interacción'}, 'Location','best');
title('p-valores Regresión Lineal por tamaño de red');
xlabel('Tamaño de red (n)');
xlim([0.5 numel(sizes)+0.5]);
set(gca,'XTickLabel',string(sizes));
ylabel('p-valor');
saveas(fig, fullfile(resultsDir,'pvals_lineal_por_n.png'));

%% 6 · Heatmaps
tt = unique(T.total_trolls);
pv = linspace(0,1,6);
Z = nan(numel(tt), numel(pv)-1);
for i = 1:numel(tt)
  for j = 1:numel(pv)-1
    m = T.total_trolls==tt(i) & T.prop_trolls_pos>=pv(j) & T.prop_trolls_pos<pv(j+1);
    Z(i,j) = mean(T.convTime(m), 'omitnan');
  end
end

figHeat = figure('Visible','off');
h = heatmap(sprintfc('%.1f', (pv(1:end-1)+pv(2:end))/2), string(tt), Z);
h.XLabel = 'Proporción trolls positivos';
h.YLabel = 'Total trolls';
h.Title  = 'Tiempo medio de convergencia';
saveas(figHeat, fullfile(resultsDir,'heatmap_convTime_prop.png'));
close(figHeat);

fprintf('\nAnálisis completo con contraste de hipótesis por tamaño de red guardado en %s\n', resultsDir);
end
