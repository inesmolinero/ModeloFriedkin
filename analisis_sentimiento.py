
import pandas as pd
from vaderSentiment.vaderSentiment import SentimentIntensityAnalyzer
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

base_path = os.getcwd() 
csv_path = os.path.join(base_path, "data", "Twitter_data.csv")

# Cargar los datos
if not os.path.exists(csv_path):
    raise FileNotFoundError(f"No se encuentra el archivo en: {csv_path}")

df = pd.read_csv(csv_path)
# Inicializar el analizador de sentimiento
analyzer = SentimentIntensityAnalyzer()

# Limpiar y analizar sentimiento
def calcular_sentimiento(texto):
    if pd.isna(texto) or texto.strip() == "":
        return None
    puntuacion = analyzer.polarity_scores(texto)
    return puntuacion["compound"]

df["sentimiento"] = df["clean_text"].apply(calcular_sentimiento)

# Quitar filas nulas
df = df.dropna(subset=["sentimiento"])

# Guardar las puntuaciones en un csv
df[["sentimiento"]].to_csv("data/Twitter_data_scores.csv", index=False)

# Crear muestra estratificada de 300 tweets diviendo en deciles
n_muestra= 300
df["decil"] = pd.qcut(df["sentimiento"], q=10, labels=False, duplicates="drop")

# Sacar 30 tweets de cada grupo 
opiniones_objetivo_por_decil = int(n_muestra/ 10)

muestra = (
    df
    .groupby("decil", group_keys=False)
    .apply(lambda x: x.sample(n=min(opiniones_objetivo_por_decil, len(x)), random_state=9),
    include_groups=False)
)
# Si no llegamos a 300, rellenamos con tweets aleatorios
if len(muestra) < n_muestra :
    faltan = 300 - len(muestra)
    resto = df[~df.index.isin(muestra.index)]
    extra = resto.sample(n=min(faltan, len(resto)), random_state=9)
    muestra = pd.concat([muestra, extra], ignore_index=True)

# Guardar muestra 
muestra[["clean_text", "sentimiento"]].to_csv(
    "data/Twitter_data_sample300.csv", 
    index=False
)

# Mostrar resultados
print(f"Total de opiniones: {len(df)}")
print(f" Total de muestra: {len(muestra)}")


# # Histograma superpuesto de proporciones relativas
plt.hist(
    df["sentimiento"],
    bins=50,
    weights=np.ones(len(df)) / len(df),
    alpha=0.4,
    color="#6A1B9A",
    edgecolor="purple",
    label="Conjunto completo"
)

plt.hist(
    muestra["sentimiento"],
    bins=50,
    weights=np.ones(len(muestra)) / len(muestra),
    alpha=0.4,
    color="#CE93D8",
    edgecolor="black",
    label="Muestra (n = 300)"
)

plt.xlabel("Valor de sentimiento")
plt.ylabel("Frecuencia relativa")
plt.title("Distribución del sentimiento: conjunto completo vs muestra")
plt.legend()
plt.grid(axis="y", alpha=0.3)

plt.show()