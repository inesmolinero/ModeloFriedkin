function idx = pickByCentralityBands(pr, k, loc, seed)
% Selecciona k nodos por bandas de centralidad.
% Entradas: 
%    pr: vector con puntuaciones de centralidad PageRank (n x 1), 
%    k: nº total de trolls a escoger, 
%    loc: banda a elegir según el nivel de centralidad que queremos(low,
%    mid o high)
%    seed: semilla (si no pasa el argumento, se establece en 42)
% Salidas:
%   idx: vector de índices 
    if nargin < 4, seed = 42;  end            % valor por defecto si no lo pasanend
    
    n = numel(pr);
    [~, ord] = sort(pr, 'ascend');  % ordenación de nodos por orden ascendente
    terc = floor(n/3); % se establecen 3 bandas de centralidad.
    switch char(loc)
        case 'bajoPR'
            pool = ord(1:terc);
        case 'medioPR'
            pool = ord(terc+1:2*terc);
        case 'altoPR'
            pool = ord(2*terc+1:end);
        otherwise
            pool = ord; % por si acaso
    end

    %  barajado reproducible del pool 
    st = rng;                    % guarda estado global
    rng(seed,'twister');         % fija semilla/algoritmo reproducible
    pool = pool(randperm(numel(pool)));  % baraja estable por semilla
    rng(st);                     % restaura estado global


    % selección determinística: toma los primeros k del pool
    k = min(k, numel(pool));
    idx = pool(1:k);
end
