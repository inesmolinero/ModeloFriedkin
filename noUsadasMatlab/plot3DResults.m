function plot3DResults(csvPath)
    % Lee la tabla de resultados
    T = readtable(csvPath);
    % Calcula #trolls = neg + pos
    T.NumTrolls = T.Neg + T.Pos;

    % Vectores únicos de p y trolls
    ps = unique(T.p);
    ts = unique(T.NumTrolls);

    % Matriz con RangoFinal promedio
    Z = nan(numel(ts), numel(ps));
    for i = 1:numel(ts)
      for j = 1:numel(ps)
        idx = T.NumTrolls==ts(i) & T.p==ps(j);
        Z(i,j) = mean(T.RangoFinal(idx));
      end
    end

    % 3D surf
    [PGrid,TGrid] = meshgrid(ps, ts);
    figure;
    surf(PGrid, TGrid, Z, 'EdgeColor','none');
    xlabel('p (conectividad)'); ylabel('# Trolls'); zlabel('RangoFinal');
    title('RangoFinal medio vs p y nº de trolls');
    colorbar;
    view(45,30);
end
