%% === 0) Ruta fija + lectura robusta ===
resultsDir = 'C:\Users\inesm\OneDrive\Documentos\tfg\ModeloFriedkin\resultados_2025-10-10_190104';

% Nombre esperado (ajústalo si es otro)
fname = fullfile(resultsDir, 'summary_metrics_AGG.csv');

% Si no existe, intenta localizar cualquier CSV parecido
if ~exist(fname,'file')
    d = dir(fullfile(resultsDir, '*.csv'));
    if isempty(d)
        error('No se encontró ningún CSV en: %s', resultsDir);
    else
        % coge el primero (o cambia el criterio si quieres uno concreto)
        fname = fullfile(resultsDir, d(1).name);
        fprintf('Usando CSV encontrado: %s\n', fname);
    end
end

% Intento 1: inferir con detectImportOptions
try
    opts = detectImportOptions(fname);
    T = readtable(fname, opts);
catch
    % Intento 2: probar coma y luego punto y coma
    try
        T = readtable(fname, 'Delimiter', ',');
    catch
        T = readtable(fname, 'Delimiter', ';');
    end
end

fprintf('Leído: %s (filas=%d, columnas=%d)\n', fname, height(T), width(T));

%% === Carpeta de resultados (tu ruta fija) ===
resultsDir = 'C:\Users\inesm\OneDrive\Documentos\tfg\ModeloFriedkin\resultados_2025-10-10_190104';
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
diary(fullfile(resultsDir,'log.txt'));

%% === 1) Agregado por lambda (barras) ===
reqLam = {'lam','prop_consenso','prop_polarizado','prop_fragmentacion','mean_K_clusters','mean_gapMin','mean_centersRange'};
if all(ismember(reqLam, T.Properties.VariableNames))
    G_lam = groupsummary(T, "lam", "mean", reqLam(2:end));
    G_lam = sortrows(G_lam, "lam");

    fig1 = figure('Name','Barras por lambda');
    bar(G_lam.lam, [G_lam.mean_prop_consenso, G_lam.mean_prop_polarizado, G_lam.mean_prop_fragmentacion], 'grouped');
    xlabel('\lambda'); ylabel('Proporción media');
    legend({'consenso','polarización','fragmentación'}, 'Location','best');
    title('Resultados medios por \lambda'); grid on;
    exportgraphics(fig1, fullfile(resultsDir,'grafico_barras_lambda.png'), 'Resolution',300);
else
    warning('Faltan columnas para gráfico por lambda. Omitido.');
end

%% === 2) Agregado por fracTrolls (líneas) ===
reqFrac = {'fracTrolls','prop_consenso','prop_polarizado','prop_fragmentacion','mean_K_clusters','mean_gapMin','mean_centersRange'};
if all(ismember(reqFrac, T.Properties.VariableNames))
    G_frac = groupsummary(T, "fracTrolls", "mean", reqFrac(2:end));
    G_frac = sortrows(G_frac, "fracTrolls");

    fig2 = figure('Name','Líneas por fracTrolls');
    plot(G_frac.fracTrolls, G_frac.mean_prop_consenso, '-o'); hold on;
    plot(G_frac.fracTrolls, G_frac.mean_prop_polarizado, '-o');
    plot(G_frac.fracTrolls, G_frac.mean_prop_fragmentacion, '-o'); hold off;
    xlabel('fracTrolls'); ylabel('Proporción media');
    legend({'consenso','polarización','fragmentación'}, 'Location','best');
    title('Proporciones por fracción de trolls'); grid on;
    exportgraphics(fig2, fullfile(resultsDir,'grafico_lineas_fractrolls.png'), 'Resolution',300);
else
    warning('Faltan columnas para gráfico por fracTrolls. Omitido.');
end

%% === 3) Dispersión gap vs rango de centros (color por outcome dominante) ===
needScatter = {'prop_consenso','prop_polarizado','prop_fragmentacion','mean_gapMin','mean_centersRange'};
if all(ismember(needScatter, T.Properties.VariableNames))
    props = [T.prop_consenso, T.prop_polarizado, T.prop_fragmentacion];
    [~, idxMax] = max(props, [], 2);
    labels = strings(height(T),1);
    labels(idxMax==1) = "consenso";
    labels(idxMax==2) = "polarización";
    labels(idxMax==3) = "fragmentación";

    fig3 = figure('Name','Scatter gap vs range');
    gscatter(T.mean_gapMin, T.mean_centersRange, labels);
    xlabel('mean\_gapMin'); ylabel('mean\_centersRange');
    title('Separación entre centros vs rango (outcome dominante)'); grid on;
    exportgraphics(fig3, fullfile(resultsDir,'scatter_gap_vs_range.png'), 'Resolution',300);
else
    warning('Faltan columnas para scatter gap vs range. Omitido.');
end

%% === 4) Top 5 escenarios consenso / polarización (tablas y export a Excel) ===
cols = intersect({'tag','regimen','loc','n','neg','pos','fracTrolls','lam', ...
                  'mean_K_clusters','prop_polarizado','prop_fragmentacion', ...
                  'prop_consenso','mean_gapMin','mean_centersRange','mean_rangoFinal','mean_stdFinal'}, ...
                  T.Properties.VariableNames);

if any(strcmp('prop_polarizado', cols))
    TopPolar = sortrows(T(:,cols), 'prop_polarizado', 'descend');
    TopPolar = TopPolar(1:min(5,height(TopPolar)), :);
else
    TopPolar = table();
end

if any(strcmp('prop_consenso', cols))
    TopCons  = sortrows(T(:,cols), 'prop_consenso', 'descend');
    TopCons  = TopCons(1:min(5,height(TopCons)), :);
else
    TopCons = table();
end

% Exporta tablas agregadas y tops a Excel dentro de tu carpeta
outXLSX = fullfile(resultsDir,'resumen_analisis.xlsx');
if exist('G_lam','var'),  writetable(G_lam,  outXLSX, 'Sheet','by_lambda'); end
if exist('G_frac','var'), writetable(G_frac, outXLSX, 'Sheet','by_fracTrolls'); end
if ~isempty(TopPolar),    writetable(TopPolar, outXLSX, 'Sheet','top_polarizados'); end
if ~isempty(TopCons),     writetable(TopCons,  outXLSX, 'Sheet','top_consenso');   end

disp('Listo: gráficos exportados (.png) y Excel creado en tu carpeta de resultados.');
diary off
