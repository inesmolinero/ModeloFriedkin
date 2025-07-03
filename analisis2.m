function analisis2(resultsDir)
% analisis2.m – EDA, reducción dimensional y análisis supervisado para trolls

if nargin<1
    resultsDir = uigetdir(pwd, 'Selecciona carpeta con summary_metrics.csv');
    if isequal(resultsDir,0)
        error('No seleccionaste ninguna carpeta.');
    end
end

%% 1 · Carga y preprocesado
T = readtable(fullfile(resultsDir,'summary_metrics.csv'));
% Unificar nombres de columnas
if ismember('total_trolls', T.Properties.VariableNames)
    T.totalTrolls = T.total_trolls;
end
T.ratioPos = T.pos ./ (T.totalTrolls + eps);
T.ratioNeg = T.neg ./ (T.totalTrolls + eps);

%% 2 · Cargar opiniones finales y fragmentación
tags = cellstr(string(T.tag));
Xend = loadOpinions(resultsDir, tags);
silh = nan(height(T),1);
for i = 1:height(T)
    x = Xend{i};
    if numel(x) < 2
        continue;
    end
    % k-means robusto con iteraciones aumentadas y un solo replicate
    try
        idx = kmeans(x(:), 2, 'Replicates',1, 'MaxIter',500, 'Options', statset('Display','off'));
        s = silhouette(x(:), idx);
        silh(i) = mean(s);
    catch ME
        warning('Clustering falló para %s: %s', tags{i}, ME.message);
    end
end
T.fragment = silh;

%% 3 · EDA básico sobre Z-scores
vars = {'rangoFinal','stdFinal','convTime','totalTrolls','ratioPos','fragment'};
Zraw = T{:,vars};
mu = nanmean(Zraw,1);
sigma = nanstd(Zraw,0,1);
Z = (Zraw - mu) ./ sigma;
valid = all(~isnan(Z),2);

if sum(valid) >= 2
    % 3.1 Matriz de correlaciones
    C = corrcoef(Z(valid,:), 'Rows','complete');
    figure; heatmap(vars,vars,C); colormap(parula);
    title('Matriz de Correlaciones Z-score');

    % 3.2 Boxplots
    figure; boxplot(Z(valid,:), 'Labels', vars);
    title('Boxplots de Z-score');

    % 3.3 Clustering jerárquico
    distM = pdist(Z(valid,:), 'euclidean');
    linkM = linkage(distM, 'ward');
    figure; dendrogram(linkM);
    title('Dendrograma de escenarios');
else
    warning('No hay suficientes datos válidos para EDA.');
end

%% 4 · Reducción de dimensionalidad
if sum(valid) >= 2
    % 4.1 PCA
    [~,~,~,~,expl] = pca(Z(valid,:), 'Rows','complete');
    figure; pareto(expl); title('Scree plot PCA');

    % 4.2 EFA exploratorio
    K = min(3, sum(valid)-1);
    if K >= 1
        try
            [Load, ~] = factoran(Z(valid,:), K, 'rotate','varimax');
            figure; imagesc(Load); colorbar;
            set(gca, 'XTick', 1:K, 'XTickLabel', strcat('F', string(1:K)), ...
                     'YTick', 1:numel(vars), 'YTickLabel', vars);
            title('Cargas factoriales');
        catch ME
            warning('EFA falló: %s');
        end
    end
end

%% 5 · Análisis supervisado
features = [T.totalTrolls, T.ratioPos, T.p];
featNames = {'totalTrolls','ratioPos','p'};

% 5.1 GLM polarización
if any(~isnan(T.polarizado))
    glmPol = fitglm(features, T.polarizado, ...
        'Distribution','binomial', 'VarNames',[featNames,'polarizado']);
    disp(glmPol);
else
    warning('Sin datos para GLM polarizado.');
end

% 5.2 Regresión convTime y Random Forest
validTime = ~isnan(T.convTime);
if sum(validTime) >= 2
    % Regresión lineal
    mdlTime = fitlm(features(validTime,:), T.convTime(validTime), ...
        'VarNames',[featNames,'convTime']);
    disp(mdlTime);

    % Random Forest con importancia OOB habilitada
        % Random Forest con importancia OOB habilitada y nombres de predictores
    RF = TreeBagger(200, features(validTime,:), T.convTime(validTime), ...
        'Method','regression','OOBPrediction','On','OOBPredictorImportance','On', ...
        'PredictorNames', featNames);
    impOOB = RF.OOBPermutedPredictorDeltaError;
    figure; bar(impOOB);
    set(gca,'XTick',1:numel(featNames),'XTickLabel',featNames);
    title('Importancia OOB - convTime');
    impOOB = RF.OOBPermutedPredictorDeltaError;
    figure; bar(impOOB);
    set(gca,'XTick',1:numel(featNames),'XTickLabel',featNames);
    title('Importancia OOB - convTime');
end

%% 6 · Visualizaciones adicionales
if exist('RF','var')
    figure; plotPartialDependence(RF,{'ratioPos'});
    title('PDP convTime vs ratioPos');

    TT = unique(T.totalTrolls(validTime));
    RP = linspace(0,1,20);
    [G1,G2] = meshgrid(TT,RP);
    Zpred = nan(size(G1));
    for ii = 1:size(G1,1)
        for jj = 1:size(G1,2)
            Zpred(ii,jj) = predict(RF,[G1(ii,jj),G2(ii,jj),mean(T.p)]);
        end
    end
    figure; surf(G1,G2,Zpred,'EdgeColor','none');
    xlabel('TotalTrolls'); ylabel('ratioPos'); zlabel('convTime');
    title('ConvTime esperado');
end

%% 7 · Guardar resultados
writetable(T, fullfile(resultsDir,'summary_with_features.csv'));
if exist('RF','var')
    saveas(gcf, fullfile(resultsDir,'surfaceConvTime.png'));
end
fprintf('Análisis guardado en %s\n', resultsDir);
end


function Xend = loadOpinions(resultsDir, tags)
% Carga X(:,end) de cada archivo .mat usando safeName
n = numel(tags);
Xend = cell(n,1);
for i = 1:n
    pref = safeName(tags{i});
    files = dir(fullfile(resultsDir, [pref '*.mat']));
    if isempty(files)
        warning('No se encontró MAT (%s*.mat)', pref);
        Xend{i} = [];
    else
        data = load(fullfile(resultsDir, files(1).name));
        if isfield(data,'X')
            Xend{i} = data.X(:,end);
        else
            warning('MAT sin X: %s', files(1).name);
            Xend{i} = [];
        end
    end
end
end
