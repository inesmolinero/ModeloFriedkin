# Modelo de Friedkin–Johnsen con agentes intransigentes

Este repositorio contiene la implementación en MATLAB de un modelo de dinámica de opiniones basado en el modelo de Friedkin–Johnsen. El proyecto permite estudiar cómo la introducción de una minoría de agentes con opiniones fijas y opuestas a la población puede alterar la evolución de las opiniones colectivas.

---

## Contexto académico

Este proyecto ha sido desarrollado en el marco del **Trabajo de Fin de Grado en Matemáticas**, con un enfoque en sistemas dinámicos comlejos y modelización matemática de fenómenos sociales, concretamente en la evolución de opiniones políticas.

---

## Descripción del modelo

Se considera una sociedad compuesta por:
- **Agentes normales o susceptibles**, cuyas opiniones pueden evolucionar en función de la interacción social.
- **Agentes intransigentes o trolls**, caracterizados por mantener una opinión fija y opuesta, actuando como una fuente externa persistente dentro del sistema.

---

## Estructura del repositorio

- `main.m`: script principal para la ejecución de las simulaciones.
- `analisis_final.ipynb`: análisis estadístico en Python de los resultados obtenidos en MATLAB de las simulaciones del modelo
- `analisis_sentimiento.py`: análisis de sentimiento realizado en Python con el objetivo de obtener las opiniones iniciales a partir de un dataset público de Kaggle. 
- `resultados_YYYY-MM-DD_HHMMSS\`: resultados de la simulación para la fecha en la que se ha ejecutado el script
---

## Ejecución
Para ejecutar las simulaciones y el análisis de sentimiento, abrir MATLAB y ejecutar `main.mat`.

Los resultados de las simulaciones se guardarán en la carpeta `resultados_YYYY-MM-DD_HHMMSS\` 

Para obtener los análisis correspondientes de estos resultados, abrir Python, cambiar las rutas correspondientes y ejecutar `analisis_final.ipynb`. 

