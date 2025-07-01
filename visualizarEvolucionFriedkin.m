%===========================================================================  
% visualizarEvolucionFriedkin.m  
%===========================================================================  
function fig = visualizarEvolucionFriedkin(x_hist, trolls_neg, trolls_pos, normales, opinionRange)
    % Visualiza evolución de opiniones diferenciando trolls negativos,
    % trolls positivos y agentes normales con colores distintos.
    %
    % INPUTS:
    %   x_hist       - matriz [n × iters] de opiniones a lo largo del tiempo
    %   trolls_neg   - vector de índices de trolls negativos
    %   trolls_pos   - vector de índices de trolls positivos
    %   normales     - vector de índices de agentes normales
    %   opinionRange - vector [min, max] del rango de opiniones
    %
    % OUTPUT:
    %   fig          - handle de la figura generada (oculta)

    iters = size(x_hist, 2);
    t     = 1:iters;
    
    fig = figure('Visible','off'); hold on;
    
    % Normales en gris claro
    if ~isempty(normales)
        plot(t, x_hist(normales,:)', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0);
    end
    
    % Trolls negativos en rojo
    if ~isempty(trolls_neg)
        plot(t, x_hist(trolls_neg,:)', 'r--', 'LineWidth', 1.5);
    end
    
    % Trolls positivos en azul
    if ~isempty(trolls_pos)
        plot(t, x_hist(trolls_pos,:)', 'b-.', 'LineWidth', 1.5);
    end
    
    % Líneas de referencia en los extremos del rango y en 0
    yline(opinionRange(1), 'k:');
    yline(0,                 'k:');
    yline(opinionRange(2), 'k:');
    
    xlabel('Iteraciones');
    ylabel('Opinión');
    
    % Leyenda
    legendEntries = {};
    if ~isempty(normales),    legendEntries{end+1} = 'Normales';     end
    if ~isempty(trolls_neg),  legendEntries{end+1} = 'Trolls neg';   end
    if ~isempty(trolls_pos),  legendEntries{end+1} = 'Trolls pos';   end
    legend(legendEntries, 'Location', 'eastoutside');
    
    grid on; hold off;
end

