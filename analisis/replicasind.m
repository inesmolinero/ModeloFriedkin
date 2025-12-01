%% analyze_replicas.m - Análisis Detallado de Réplicas Individuales
% 
% Este script complementa el análisis en Python analizando réplicas 
% individuales guardadas en archivos .mat
%
% Uso:
%   1. Configurar resultsDir (directorio de resultados)
%   2. Ejecutar secciones para diferentes análisis

%% LIMPIEZA Y CONFIGURACIÓN
clear; close all; clc;

fprintf('=================================================================\n');
fprintf('     ANÁLISIS DE RÉPLICAS INDIVIDUALES - FRIEDKIN-JOHNSEN\n');
fprintf('=================================================================\n\n');

%% 1. CONFIGURACIÓN
% Solicitar directorio de resultados
resultsDir = input('Directorio de resultados: ', 's');

if isempty(resultsDir) || ~exist(resultsDir, 'dir')
    error('Directorio no válido: %s', resultsDir);
end

runsDir = fullfile(resultsDir, 'runs');
if ~exist(runsDir, 'dir')
    error('No se encuentra el directorio de runs: %s', runsDir);
end

% Crear directorio para análisis de réplicas
replicaAnalysisDir = fullfile(resultsDir, 'replica_analysis');
if ~exist(replicaAnalysisDir, 'dir')
    mkdir(replicaAnalysisDir);
end

fprintf('✓ Directorio de resultados: %s\n', resultsDir);
fprintf('✓ Directorio de análisis: %s\n\n', replicaAnalysisDir);

%% 2. CARGAR ÍNDICE DE RÉPLICAS
indexFile = fullfile(runsDir, 'runs_index.csv');
if ~exist(indexFile, 'file')
    error('No se encuentra runs_index.csv');
end

indexTable = readtable(indexFile);
fprintf('✓ Cargado índice con %d réplicas\n\n', height(indexTable));

% Obtener lista de escenarios únicos
scenarios = unique(indexTable.tag);
nScenarios = numel(scenarios);

fprintf('Total de escenarios: %d\n', nScenarios);
fprintf('Réplicas por escenario: ~%.0f\n\n', height(indexTable) / nScenarios);

%% 3. SELECCIONAR ESCENARIO PARA ANÁLISIS DETALLADO
fprintf('=================================================================\n');
fprintf('SELECCIÓN DE ESCENARIO PARA ANÁLISIS DETALLADO\n');
fprintf('=================================================================\n\n');

fprintf('Primeros 10 escenarios disponibles:\n');
for i = 1:min(10, nScenarios)
    fprintf('  %2d) %s\n', i, scenarios{i});
end
fprintf('\n');

scenIdx = input('Seleccione número de escenario (1-10): ');
if isempty(scenIdx) || scenIdx < 1 || scenIdx > min(10, nScenarios)
    scenIdx = 1;
    fprintf('Usando escenario por defecto: %d\n', scenIdx);
end

selectedScenario = scenarios{scenIdx};
fprintf('\n✓ Escenario seleccionado: %s\n\n', selectedScenario);

% Cargar todas las réplicas de este escenario
scenReplicas = indexTable(strcmp(indexTable.tag, selectedScenario), :);
nReps = height(scenReplicas);

fprintf('Cargando %d réplicas...\n', nReps);

replicas = cell(nReps, 1);
for r = 1:nReps
    load(scenReplicas.path{r}, 'out');
    replicas{r} = out;
end

fprintf('✓ Réplicas cargadas\n\n');

%% 4. ANÁLISIS DE TRAYECTORIAS
fprintf('=================================================================\n');
fprintf('ANÁLISIS DE TRAYECTORIAS\n');
fprintf('=================================================================\n\n');

% Seleccionar algunas réplicas representativas
repsToPlot = min(5, nReps);
repIndices = round(linspace(1, nReps, repsToPlot));

figure('Position', [100, 100, 1400, 900]);
sgtitle(sprintf('Trayectorias de Opiniones: %s', strrep(selectedScenario, '_', '\_')), ...
        'FontSize', 14, 'FontWeight', 'bold');

for idx = 1:repsToPlot
    r = repIndices(idx);
    out = replicas{r};
    
    % Snapshots
    subplot(2, 3, idx);
    X_snap = out.snapshots.X;
    tgrid = out.snapshots.tgrid;
    
    % Plotear trayectorias de la muestra
    X_sample = out.sample.X;
    
    plot(1:size(X_sample, 2), X_sample', 'LineWidth', 0.5, 'Color', [0.7 0.7 0.7 0.3]);
    hold on;
    
    % Resaltar normales vs trolls
    idx_norm = out.indices.normales;
    idx_troll = out.indices.trolls;
    
    if ~isempty(idx_troll)
        % Trayectorias de algunos trolls
        troll_sample = intersect(out.sample.idx, idx_troll);
        if ~isempty(troll_sample)
            for ti = troll_sample(1:min(3, length(troll_sample))).'
                loc = find(out.sample.idx == ti);
                plot(1:size(X_sample, 2), X_sample(loc, :), 'r-', 'LineWidth', 2);
            end
        end
    end
    
    % Línea de consenso
    if out.final.consenso
        yline(mean(out.final.z), 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('Iteración');
    ylabel('Opinión');
    title(sprintf('Réplica %d (Consenso: %s)', r, ...
          ifelse(out.final.consenso, 'Sí', 'No')));
    grid on;
    xlim([1, size(X_sample, 2)]);
    ylim([0, 1]);
    hold off;
end

% Gráfico combinado de estadísticos
subplot(2, 3, 6);
hold on;

for r = 1:min(10, nReps)
    out = replicas{r};
    ts = out.time_stats;
    
    plot(ts.t, ts.mean, 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5 0.3]);
end

xlabel('Iteración');
ylabel('Opinión Media');
title('Evolución de Opinión Media (10 réplicas)');
grid on;
hold off;

saveas(gcf, fullfile(replicaAnalysisDir, sprintf('trajectories_%s.png', selectedScenario)));
fprintf('✓ Guardado: trajectories_%s.png\n\n', selectedScenario);

%% 5. ANÁLISIS DE DISTRIBUCIONES FINALES
fprintf('=================================================================\n');
fprintf('ANÁLISIS DE DISTRIBUCIONES FINALES\n');
fprintf('=================================================================\n\n');

figure('Position', [100, 100, 1400, 600]);
sgtitle(sprintf('Distribuciones Finales: %s', strrep(selectedScenario, '_', '\_')), ...
        'FontSize', 14, 'FontWeight', 'bold');

% Recopilar todas las opiniones finales
all_z = [];
for r = 1:nReps
    all_z = [all_z; replicas{r}.final.z]; %#ok<AGROW>
end

% Histograma global
subplot(1, 3, 1);
histogram(all_z, 30, 'Normalization', 'probability', 'FaceColor', [0.3 0.6 0.9]);
xlabel('Opinión Final');
ylabel('Probabilidad');
title('Distribución Agregada (Todas las Réplicas)');
grid on;

% Boxplots por tipo de agente
subplot(1, 3, 2);
out1 = replicas{1};
n = length(out1.final.z);

z_by_type = [];
labels_by_type = [];

for r = 1:nReps
    out = replicas{r};
    z_norm = out.final.z(out.indices.normales);
    z_troll = out.final.z(out.indices.trolls);
    
    z_by_type = [z_by_type; z_norm; z_troll]; %#ok<AGROW>
    labels_by_type = [labels_by_type; ...
                     repmat({'Normal'}, length(z_norm), 1); ...
                     repmat({'Troll'}, length(z_troll), 1)]; %#ok<AGROW>
end

if ~isempty(z_by_type)
    boxplot(z_by_type, labels_by_type);
    ylabel('Opinión Final');
    title('Distribución por Tipo de Agente');
    grid on;
end

% Evolución de varianza
subplot(1, 3, 3);
var_evolution = zeros(nReps, 1);
for r = 1:nReps
    out = replicas{r};
    T = length(out.time_stats.t);
    var_evolution(r) = out.time_stats.std(end)^2;
end

histogram(var_evolution, 20, 'Normalization', 'probability', 'FaceColor', [0.9 0.4 0.4]);
xlabel('Varianza Final');
ylabel('Probabilidad');
title('Distribución de Varianza Final');
grid on;

saveas(gcf, fullfile(replicaAnalysisDir, sprintf('distributions_%s.png', selectedScenario)));
fprintf('✓ Guardado: distributions_%s.png\n\n', selectedScenario);

%% 6. ANÁLISIS DE CONVERGENCIA
fprintf('=================================================================\n');
fprintf('ANÁLISIS DE CONVERGENCIA\n');
fprintf('=================================================================\n\n');

figure('Position', [100, 100, 1400, 600]);
sgtitle(sprintf('Convergencia: %s', strrep(selectedScenario, '_', '\_')), ...
        'FontSize', 14, 'FontWeight', 'bold');

% Recopilar datos de convergencia
conv_times = nan(nReps, 1);
consensus_achieved = false(nReps, 1);

for r = 1:nReps
    out = replicas{r};
    conv_times(r) = out.final.tconv;
    consensus_achieved(r) = out.final.consenso;
end

% Histograma de tiempos de convergencia
subplot(1, 3, 1);
histogram(conv_times(~isnan(conv_times)), 20, 'FaceColor', [0.4 0.7 0.4]);
xlabel('Tiempo de Convergencia');
ylabel('Frecuencia');
title(sprintf('Distribución (%.1f%% convergen)', 100*mean(~isnan(conv_times))));
grid on;

% Proporción de consenso
subplot(1, 3, 2);
pie([sum(consensus_achieved), sum(~consensus_achieved)], ...
    {'Consenso', 'No Consenso'});
title(sprintf('Consenso Alcanzado: %.1f%%', 100*mean(consensus_achieved)));

% Evolución del rango temporal
subplot(1, 3, 3);
hold on;
for r = 1:min(20, nReps)
    out = replicas{r};
    ts = out.time_stats;
    
    color = ifelse(out.final.consenso, [0.2 0.7 0.2 0.3], [0.7 0.2 0.2 0.3]);
    plot(ts.t, ts.rango, 'LineWidth', 1, 'Color', color);
end

xlabel('Iteración');
ylabel('Rango de Opiniones');
title('Evolución del Rango (20 réplicas)');
legend({'Consenso', 'No Consenso'}, 'Location', 'best');
grid on;
hold off;

saveas(gcf, fullfile(replicaAnalysisDir, sprintf('convergence_%s.png', selectedScenario)));
fprintf('✓ Guardado: convergence_%s.png\n\n', selectedScenario);

%% 7. ANÁLISIS ESPECTRAL
fprintf('=================================================================\n');
fprintf('ANÁLISIS ESPECTRAL\n');
fprintf('=================================================================\n\n');

figure('Position', [100, 100, 1200, 500]);
sgtitle(sprintf('Propiedades Espectrales: %s', strrep(selectedScenario, '_', '\_')), ...
        'FontSize', 14, 'FontWeight', 'bold');

% Recopilar radios espectrales
rho_W_all = zeros(nReps, 1);
rho_LW_all = zeros(nReps, 1);

for r = 1:nReps
    out = replicas{r};
    rho_W_all(r) = out.spectral.rho_W;
    rho_LW_all(r) = out.spectral.rho_LW;
end

% Scatter ρ(W) vs ρ(ΛW)
subplot(1, 2, 1);
scatter(rho_W_all, rho_LW_all, 50, conv_times, 'filled', 'MarkerEdgeColor', 'k');
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
xlabel('ρ(W)');
ylabel('ρ(ΛW)');
title('Radio Espectral: W vs ΛW');
colorbar('Label', 'Tiempo Convergencia');
grid on;
axis equal;
xlim([0 1]);
ylim([0 1]);
hold off;

% Histograma de ρ(ΛW)
subplot(1, 2, 2);
histogram(rho_LW_all, 20, 'FaceColor', [0.6 0.4 0.8]);
xlabel('ρ(ΛW)');
ylabel('Frecuencia');
title(sprintf('Radio Espectral de ΛW (μ=%.4f)', mean(rho_LW_all)));
grid on;

saveas(gcf, fullfile(replicaAnalysisDir, sprintf('spectral_%s.png', selectedScenario)));
fprintf('✓ Guardado: spectral_%s.png\n\n', selectedScenario);

%% 8. ANÁLISIS DE POLARIZACIÓN
fprintf('=================================================================\n');
fprintf('ANÁLISIS DE POLARIZACIÓN\n');
fprintf('=================================================================\n\n');

figure('Position', [100, 100, 1400, 500]);
sgtitle(sprintf('Métricas de Polarización: %s', strrep(selectedScenario, '_', '\_')), ...
        'FontSize', 14, 'FontWeight', 'bold');

% Recopilar métricas
NDI_all = zeros(nReps, 1);
P2_all = zeros(nReps, 1);
P4_all = zeros(nReps, 1);

for r = 1:nReps
    out = replicas{r};
    NDI_all(r) = out.final.NDI;
    P2_all(r) = out.final.P2;
    P4_all(r) = out.final.P4;
end

% NDI
subplot(1, 3, 1);
histogram(NDI_all, 20, 'FaceColor', [0.9 0.5 0.3]);
xlabel('NDI');
ylabel('Frecuencia');
title(sprintf('NDI (μ=%.4f, σ=%.4f)', mean(NDI_all), std(NDI_all)));
grid on;

% P2
subplot(1, 3, 2);
histogram(P2_all, 20, 'FaceColor', [0.5 0.3 0.9]);
xlabel('P2 (Energía)');
ylabel('Frecuencia');
title(sprintf('P2 (μ=%.4f, σ=%.4f)', mean(P2_all), std(P2_all)));
grid on;

% P4
subplot(1, 3, 3);
histogram(P4_all, 20, 'FaceColor', [0.3 0.7 0.9]);
xlabel('P4 (Extremismo)');
ylabel('Frecuencia');
title(sprintf('P4 (μ=%.1f, σ=%.1f)', mean(P4_all), std(P4_all)));
grid on;

saveas(gcf, fullfile(replicaAnalysisDir, sprintf('polarization_%s.png', selectedScenario)));
fprintf('✓ Guardado: polarization_%s.png\n\n', selectedScenario);

%% 9. COMPARACIÓN ENTRE ESCENARIOS (MUESTRA)
fprintf('=================================================================\n');
fprintf('COMPARACIÓN ENTRE ESCENARIOS\n');
fprintf('=================================================================\n\n');

% Seleccionar 4 escenarios diferentes para comparar
nScenariosCompare = min(4, nScenarios);
scenIndices = round(linspace(1, nScenarios, nScenariosCompare));

figure('Position', [100, 100, 1400, 900]);
sgtitle('Comparación entre Escenarios', 'FontSize', 14, 'FontWeight', 'bold');

for idx = 1:nScenariosCompare
    scenTag = scenarios{scenIndices(idx)};
    scenReps = indexTable(strcmp(indexTable.tag, scenTag), :);
    
    % Cargar primera réplica de este escenario
    load(scenReps.path{1}, 'out');
    
    subplot(2, 2, idx);
    
    % Plotear snapshots de opiniones
    X_snap = out.snapshots.X;
    imagesc(out.snapshots.tgrid, 1:size(X_snap, 1), X_snap);
    colormap(jet);
    colorbar;
    caxis([0 1]);
    
    xlabel('Iteración');
    ylabel('Agente');
    title(strrep(scenTag, '_', '\_'), 'FontSize', 10);
end

saveas(gcf, fullfile(replicaAnalysisDir, 'scenario_comparison.png'));
fprintf('✓ Guardado: scenario_comparison.png\n\n', selectedScenario);

%% 10. RESUMEN ESTADÍSTICO
fprintf('=================================================================\n');
fprintf('RESUMEN ESTADÍSTICO DEL ESCENARIO SELECCIONADO\n');
fprintf('=================================================================\n\n');

% Crear tabla resumen
summary = struct();
summary.Escenario = selectedScenario;
summary.NumReplicas = nReps;
summary.NDI_mean = mean(NDI_all);
summary.NDI_std = std(NDI_all);
summary.P2_mean = mean(P2_all);
summary.P2_std = std(P2_all);
summary.P4_mean = mean(P4_all);
summary.P4_std = std(P4_all);
summary.ConvTime_mean = mean(conv_times, 'omitnan');
summary.ConvTime_std = std(conv_times, 'omitnan');
summary.ProporcionConsenso = mean(consensus_achieved);
summary.RhoLW_mean = mean(rho_LW_all);
summary.RhoLW_std = std(rho_LW_all);

% Mostrar resumen
fprintf('Escenario: %s\n', summary.Escenario);
fprintf('Réplicas: %d\n\n', summary.NumReplicas);

fprintf('POLARIZACIÓN:\n');
fprintf('  NDI:  %.4f ± %.4f\n', summary.NDI_mean, summary.NDI_std);
fprintf('  P2:   %.4f ± %.4f\n', summary.P2_mean, summary.P2_std);
fprintf('  P4:   %.2f ± %.2f\n\n', summary.P4_mean, summary.P4_std);

fprintf('CONVERGENCIA:\n');
fprintf('  Tiempo medio: %.2f ± %.2f iteraciones\n', ...
        summary.ConvTime_mean, summary.ConvTime_std);
fprintf('  Consenso: %.1f%%\n\n', 100*summary.ProporcionConsenso);

fprintf('ESPECTRAL:\n');
fprintf('  ρ(ΛW): %.4f ± %.4f\n\n', summary.RhoLW_mean, summary.RhoLW_std);

% Guardar resumen
summaryFile = fullfile(replicaAnalysisDir, sprintf('summary_%s.txt', selectedScenario));
fid = fopen(summaryFile, 'w');
fprintf(fid, 'RESUMEN ESTADÍSTICO - %s\n', selectedScenario);
fprintf(fid, '=================================================================\n\n');
fprintf(fid, 'Réplicas: %d\n\n', summary.NumReplicas);
fprintf(fid, 'POLARIZACIÓN:\n');
fprintf(fid, '  NDI:  %.4f ± %.4f\n', summary.NDI_mean, summary.NDI_std);
fprintf(fid, '  P2:   %.4f ± %.4f\n', summary.P2_mean, summary.P2_std);
fprintf(fid, '  P4:   %.2f ± %.2f\n\n', summary.P4_mean, summary.P4_std);
fprintf(fid, 'CONVERGENCIA:\n');
fprintf(fid, '  Tiempo medio: %.2f ± %.2f iteraciones\n', ...
        summary.ConvTime_mean, summary.ConvTime_std);
fprintf(fid, '  Consenso: %.1f%%\n\n', 100*summary.ProporcionConsenso);
fprintf(fid, 'ESPECTRAL:\n');
fprintf(fid, '  ρ(ΛW): %.4f ± %.4f\n', summary.RhoLW_mean, summary.RhoLW_std);
fclose(fid);

fprintf('✓ Resumen guardado: %s\n\n', summaryFile);

%% FINALIZACIÓN
fprintf('=================================================================\n');
fprintf('✓ ANÁLISIS DE RÉPLICAS COMPLETADO\n');
fprintf('=================================================================\n\n');

fprintf('Archivos generados en: %s\n\n', replicaAnalysisDir);

% Función auxiliar
function y = ifelse(cond, a, b)
    if cond, y = a; else, y = b; end
end