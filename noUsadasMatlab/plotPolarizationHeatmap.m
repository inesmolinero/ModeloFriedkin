function plotPolarizationHeatmap(csvPath)
    % Lee la tabla de resultados
    T = readtable(csvPath);
    % Número de trolls
    T.NumTrolls = T.Neg + T.Pos;

    % Vectores únicos
    ps = unique(T.p);
    ts = unique(T.NumTrolls);

    % Matriz con proporción de simulaciones polarizadas
    M = nan(numel(ts), numel(ps));
    for i = 1:numel(ts)
      for j = 1:numel(ps)
        idx = T.NumTrolls==ts(i) & T.p==ps(j);
        M(i,j) = mean(T.Polarizado(idx));  % true→1 / false→0
      end
    end

    % Heatmap
    figure;
    heatmap(ps, ts, M, ...
        'Colormap', parula, ...
        'ColorLimits',[0 1], ...
        'XLabel','p (conectividad)', ...
        'YLabel','# Trolls', ...
        'Title','Frac. de escenarios polarizados');
end
