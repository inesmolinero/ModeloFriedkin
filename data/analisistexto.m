% Script para analizar sentimiento con VADER en MATLAB y guardar a CSV
% Requiere Text Analytics Toolbox (R2021b+ aprox.)

% 1) Carga tu CSV original
csvFile = fullfile("data","Twitter_Data.csv");
T = readtable(csvFile, 'TextType','string');

% Asegúrate de que exista la columna 'clean_text'
if ~ismember("clean_text", T.Properties.VariableNames)
    error("No se encontró la columna 'clean_text' en el CSV.");
end

% 2) (Opcional) Preprocesado mínimo
docs = tokenizedDocument(T.clean_text, 'Language','en');

% 3) Sentimiento con método VADER
% Según versión, 'sentiment' puede devolver etiqueta y/o puntuaciones.
% Con VADER, podemos pedir las puntuaciones detalladas:
[sentLabel, sentScores] = sentiment(docs, 'SentimentMethod','vader');

% 'sentScores' suele incluir Negativity, Neutrality, Positivity y Compound
% Normalizamos a una sola columna 'sentiment_cont' usando Compound
if ismember('Compound', sentScores.Properties.VariableNames)
    sentiment_cont = sentScores.Compound;
else
    % Fallback: aproximar tipo VADER si no viene 'Compound'
    % (pos - neg), acotado [-1,1]
    if all(ismember({'Positivity','Negativity'}, sentScores.Properties.VariableNames))
        raw = sentScores.Positivity - sentScores.Negativity;
        sentiment_cont = max(-1, min(1, raw));
    else
        error("La salida de 'sentiment' no contiene 'Compound' ni Pos/Neg.");
    end
end

% 4) Quedarnos con la columna numérica
out = table(sentiment_cont);

% 5) Guardar a CSV
output_path = fullfile("data","Twitter_Data_scores.csv");
writetable(out, output_path);

% 6) Mostrar primer vistazo
head(out)
