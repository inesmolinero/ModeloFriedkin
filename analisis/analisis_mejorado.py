"""
Análisis Estadístico y Visualización de Experimentos Friedkin-Johnsen

Este script realiza:
1. Carga y limpieza de datos
2. Análisis exploratorio de datos (EDA)
3. Tests de hipótesis estadísticas
4. Visualizaciones comprehensivas
5. Generación de informe automático

Autor: Análisis mejorado
Fecha: Noviembre 2024
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from scipy.io import loadmat
import warnings
warnings.filterwarnings('ignore')

# Configuración de estilo
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 11
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['axes.labelsize'] = 12

#==============================================================================
# SECCIÓN 1: CONFIGURACIÓN Y CARGA DE DATOS
#==============================================================================

class FriedkinAnalyzer:
    """Clase principal para análisis de experimentos Friedkin-Johnsen"""
    
    def __init__(self, results_dir):
        """
        Inicializa el analizador con el directorio de resultados
        
        Args:
            results_dir: Ruta al directorio resultados_YYYY-MM-DD_HHMMSS
        """
        self.results_dir = Path(results_dir)
        self.output_dir = self.results_dir / 'analysis'
        self.output_dir.mkdir(exist_ok=True)
        
        print(f"📂 Directorio de resultados: {self.results_dir}")
        print(f"📊 Análisis se guardará en: {self.output_dir}\n")
        
        # Cargar datos
        self.df = self.load_data()
        self.prepare_data()
        
    def load_data(self):
        """Carga el CSV de métricas agregadas"""
        csv_path = self.results_dir / 'summary_metrics_AGG.csv'
        
        if not csv_path.exists():
            raise FileNotFoundError(f"No se encuentra: {csv_path}")
        
        df = pd.read_csv(csv_path)
        print(f"✓ Cargados {len(df)} escenarios")
        print(f"✓ Variables: {len(df.columns)} columnas\n")
        
        return df
    
    def prepare_data(self):
        """Prepara y limpia los datos para análisis"""
        # Crear variables categóricas ordenadas
        self.df['regimen_cat'] = pd.Categorical(
            self.df['regimen'], 
            categories=['desconectada', 'umbral', 'fuerte'],
            ordered=True
        )
        
        self.df['loc_cat'] = pd.Categorical(
            self.df['loc'],
            categories=['low', 'mid', 'high'],
            ordered=True
        )
        
        # Crear variable de proporción de trolls categórica
        self.df['trolls_pct'] = self.df['fracTrolls'] * 100
        self.df['trolls_cat'] = pd.cut(
            self.df['trolls_pct'],
            bins=[-0.1, 0.1, 15, 25, 35],
            labels=['0%', '10%', '20%', '30%']
        )
        
        # Variables de interés principales
        self.outcome_vars = [
            'mean_rangoFinal', 'mean_stdFinal', 'mean_convTime',
            'mean_NDI', 'mean_P2', 'mean_P4',
            'prop_consenso', 'mean_rho_LW'
        ]
        
        print("✓ Datos preparados\n")

#==============================================================================
# SECCIÓN 2: ANÁLISIS EXPLORATORIO DE DATOS (EDA)
#==============================================================================

    def exploratory_analysis(self):
        """Realiza análisis exploratorio completo"""
        print("="*70)
        print("ANÁLISIS EXPLORATORIO DE DATOS")
        print("="*70 + "\n")
        
        # Estadísticos descriptivos por régimen
        print("📊 Estadísticos Descriptivos por Régimen de Red\n")
        desc_by_regime = self.df.groupby('regimen')[self.outcome_vars].agg(['mean', 'std', 'min', 'max'])
        print(desc_by_regime.round(4))
        print("\n" + "-"*70 + "\n")
        
        # Estadísticos por proporción de trolls
        print("📊 Estadísticos Descriptivos por Proporción de Trolls\n")
        desc_by_trolls = self.df.groupby('trolls_cat')[self.outcome_vars].agg(['mean', 'std'])
        print(desc_by_trolls.round(4))
        print("\n" + "-"*70 + "\n")
        
        # Correlaciones
        print("📊 Matriz de Correlaciones entre Variables Clave\n")
        corr_matrix = self.df[self.outcome_vars].corr()
        print(corr_matrix.round(3))
        print("\n")
        
        # Guardar resultados
        with open(self.output_dir / 'descriptive_stats.txt', 'w') as f:
            f.write("ESTADÍSTICOS DESCRIPTIVOS - Experimentos Friedkin-Johnsen\n")
            f.write("="*70 + "\n\n")
            f.write("Por Régimen de Red:\n")
            f.write(desc_by_regime.to_string())
            f.write("\n\n" + "-"*70 + "\n\n")
            f.write("Por Proporción de Trolls:\n")
            f.write(desc_by_trolls.to_string())
            f.write("\n\n" + "-"*70 + "\n\n")
            f.write("Matriz de Correlaciones:\n")
            f.write(corr_matrix.to_string())
        
        return desc_by_regime, desc_by_trolls, corr_matrix

#==============================================================================
# SECCIÓN 3: TESTS DE HIPÓTESIS
#==============================================================================

    def hypothesis_tests(self):
        """Realiza tests de hipótesis estadísticas"""
        print("="*70)
        print("TESTS DE HIPÓTESIS ESTADÍSTICAS")
        print("="*70 + "\n")
        
        results = {}
        
        # H1: El régimen de conectividad afecta la polarización (NDI)
        print("H1: ¿El régimen de red afecta la polarización (NDI)?")
        print("-" * 60)
        h1_result = self._test_regime_effect('mean_NDI')
        results['H1_regime_polarization'] = h1_result
        print()
        
        # H2: Mayor proporción de trolls aumenta polarización
        print("H2: ¿Mayor proporción de trolls aumenta la polarización?")
        print("-" * 60)
        h2_result = self._test_trolls_correlation('mean_NDI')
        results['H2_trolls_polarization'] = h2_result
        print()
        
        # H3: Lambda (susceptibilidad) afecta convergencia
        print("H3: ¿Lambda afecta el tiempo de convergencia?")
        print("-" * 60)
        h3_result = self._test_lambda_effect('mean_convTime')
        results['H3_lambda_convergence'] = h3_result
        print()
        
        # H4: Centralidad de trolls afecta influencia
        print("H4: ¿La ubicación de trolls (centralidad) afecta el rango final?")
        print("-" * 60)
        h4_result = self._test_location_effect('mean_rangoFinal')
        results['H4_location_influence'] = h4_result
        print()
        
        # H5: Interacción régimen x trolls
        print("H5: ¿Existe interacción entre régimen y proporción de trolls?")
        print("-" * 60)
        h5_result = self._test_interaction('mean_NDI')
        results['H5_interaction'] = h5_result
        print()
        
        # H6: Probabilidad frontera vs polarización
        print("H6: ¿Mayor contacto normales-trolls aumenta polarización?")
        print("-" * 60)
        h6_result = self._test_boundary_effect('mean_P_norm_vs_trolls', 'mean_NDI')
        results['H6_boundary_polarization'] = h6_result
        print()
        
        # Guardar resultados
        self._save_hypothesis_results(results)
        
        return results
    
    def _test_regime_effect(self, outcome_var):
        """ANOVA para efecto del régimen"""
        groups = [self.df[self.df['regimen'] == r][outcome_var].dropna() 
                  for r in ['desconectada', 'umbral', 'fuerte']]
        
        # Eliminar grupos vacíos
        groups = [g for g in groups if len(g) > 0]
        
        if len(groups) < 2:
            print("⚠️  Datos insuficientes para ANOVA")
            return {'test': 'ANOVA', 'F': np.nan, 'p_value': np.nan, 'significant': False}
        
        F_stat, p_value = stats.f_oneway(*groups)
        
        # Medias por grupo
        means = {r: self.df[self.df['regimen'] == r][outcome_var].mean() 
                 for r in ['desconectada', 'umbral', 'fuerte']}
        
        print(f"  ANOVA: F = {F_stat:.4f}, p = {p_value:.4f}")
        print(f"  Medias por régimen: {means}")
        
        if p_value < 0.05:
            print("  ✓ Efecto significativo (p < 0.05)")
            significant = True
        else:
            print("  ✗ No significativo (p ≥ 0.05)")
            significant = False
        
        return {
            'test': 'ANOVA',
            'F': F_stat,
            'p_value': p_value,
            'significant': significant,
            'means': means
        }
    
    def _test_trolls_correlation(self, outcome_var):
        """Correlación de Spearman: trolls vs polarización"""
        data = self.df[['fracTrolls', outcome_var]].dropna()
        
        if len(data) < 3:
            print("⚠️  Datos insuficientes para correlación")
            return {'test': 'Spearman', 'rho': np.nan, 'p_value': np.nan, 'significant': False}
        
        rho, p_value = stats.spearmanr(data['fracTrolls'], data[outcome_var])
        
        print(f"  Spearman ρ = {rho:.4f}, p = {p_value:.4f}")
        
        if p_value < 0.05:
            direction = "positiva" if rho > 0 else "negativa"
            print(f"  ✓ Correlación {direction} significativa (p < 0.05)")
            significant = True
        else:
            print("  ✗ No significativo (p ≥ 0.05)")
            significant = False
        
        return {
            'test': 'Spearman',
            'rho': rho,
            'p_value': p_value,
            'significant': significant
        }
    
    def _test_lambda_effect(self, outcome_var):
        """Correlación: lambda vs tiempo de convergencia"""
        data = self.df[['lam', outcome_var]].dropna()
        
        if len(data) < 3:
            print("⚠️  Datos insuficientes")
            return {'test': 'Pearson', 'r': np.nan, 'p_value': np.nan, 'significant': False}
        
        r, p_value = stats.pearsonr(data['lam'], data[outcome_var])
        
        print(f"  Pearson r = {r:.4f}, p = {p_value:.4f}")
        
        if p_value < 0.05:
            direction = "positiva" if r > 0 else "negativa"
            print(f"  ✓ Correlación {direction} significativa (p < 0.05)")
            significant = True
        else:
            print("  ✗ No significativo (p ≥ 0.05)")
            significant = False
        
        return {
            'test': 'Pearson',
            'r': r,
            'p_value': p_value,
            'significant': significant
        }
    
    def _test_location_effect(self, outcome_var):
        """ANOVA: ubicación de trolls vs influencia"""
        # Filtrar solo escenarios con trolls
        data_with_trolls = self.df[self.df['trolls'] > 0]
        
        if len(data_with_trolls) == 0:
            print("⚠️  No hay escenarios con trolls")
            return {'test': 'ANOVA', 'F': np.nan, 'p_value': np.nan, 'significant': False}
        
        groups = [data_with_trolls[data_with_trolls['loc'] == loc][outcome_var].dropna()
                  for loc in ['low', 'mid', 'high']]
        
        groups = [g for g in groups if len(g) > 0]
        
        if len(groups) < 2:
            print("⚠️  Datos insuficientes")
            return {'test': 'ANOVA', 'F': np.nan, 'p_value': np.nan, 'significant': False}
        
        F_stat, p_value = stats.f_oneway(*groups)
        
        means = {loc: data_with_trolls[data_with_trolls['loc'] == loc][outcome_var].mean()
                 for loc in ['low', 'mid', 'high']}
        
        print(f"  ANOVA: F = {F_stat:.4f}, p = {p_value:.4f}")
        print(f"  Medias por ubicación: {means}")
        
        if p_value < 0.05:
            print("  ✓ Efecto significativo (p < 0.05)")
            significant = True
        else:
            print("  ✗ No significativo (p ≥ 0.05)")
            significant = False
        
        return {
            'test': 'ANOVA',
            'F': F_stat,
            'p_value': p_value,
            'significant': significant,
            'means': means
        }
    
    def _test_interaction(self, outcome_var):
        """Test de interacción: régimen x trolls (ANOVA two-way simplificado)"""
        # Preparar datos
        data = self.df[['regimen', 'trolls_cat', outcome_var]].dropna()
        
        if len(data) < 10:
            print("⚠️  Datos insuficientes para test de interacción")
            return {'test': 'Two-way ANOVA', 'interaction_p': np.nan, 'significant': False}
        
        # Crear grupos para cada combinación
        interaction_means = data.groupby(['regimen', 'trolls_cat'])[outcome_var].mean()
        
        print(f"  Medias por régimen x trolls:")
        print(interaction_means.unstack(fill_value=np.nan).round(4))
        
        # Test simplificado: comparar varianza entre vs dentro grupos
        grand_mean = data[outcome_var].mean()
        
        # Varianza explicada por interacción (aproximada)
        ss_interaction = 0
        for (reg, troll), group_data in data.groupby(['regimen', 'trolls_cat']):
            group_mean = group_data[outcome_var].mean()
            n = len(group_data)
            ss_interaction += n * (group_mean - grand_mean)**2
        
        ss_total = np.sum((data[outcome_var] - grand_mean)**2)
        
        if ss_total > 0:
            eta_squared = ss_interaction / ss_total
            print(f"  η² (efecto) = {eta_squared:.4f}")
            
            # Heurística: interacción fuerte si η² > 0.10
            if eta_squared > 0.10:
                print("  ✓ Interacción aparente fuerte (η² > 0.10)")
                significant = True
            else:
                print("  ✗ Interacción débil (η² ≤ 0.10)")
                significant = False
        else:
            eta_squared = 0
            significant = False
        
        return {
            'test': 'Two-way ANOVA (simplified)',
            'eta_squared': eta_squared,
            'significant': significant,
            'means': interaction_means.to_dict()
        }
    
    def _test_boundary_effect(self, predictor_var, outcome_var):
        """Correlación: probabilidad frontera vs polarización"""
        data = self.df[[predictor_var, outcome_var]].dropna()
        
        if len(data) < 3:
            print("⚠️  Datos insuficientes")
            return {'test': 'Spearman', 'rho': np.nan, 'p_value': np.nan, 'significant': False}
        
        rho, p_value = stats.spearmanr(data[predictor_var], data[outcome_var])
        
        print(f"  Spearman ρ = {rho:.4f}, p = {p_value:.4f}")
        
        if p_value < 0.05:
            direction = "positiva" if rho > 0 else "negativa"
            print(f"  ✓ Correlación {direction} significativa (p < 0.05)")
            significant = True
        else:
            print("  ✗ No significativo (p ≥ 0.05)")
            significant = False
        
        return {
            'test': 'Spearman',
            'rho': rho,
            'p_value': p_value,
            'significant': significant
        }
    
    def _save_hypothesis_results(self, results):
        """Guarda resultados de tests de hipótesis"""
        with open(self.output_dir / 'hypothesis_tests.txt', 'w') as f:
            f.write("TESTS DE HIPÓTESIS - Experimentos Friedkin-Johnsen\n")
            f.write("="*70 + "\n\n")
            
            for hypothesis, result in results.items():
                f.write(f"{hypothesis}:\n")
                f.write(f"  Test: {result['test']}\n")
                
                for key, value in result.items():
                    if key != 'test':
                        f.write(f"  {key}: {value}\n")
                
                f.write("\n" + "-"*70 + "\n\n")
        
        print(f"✓ Resultados guardados en: {self.output_dir / 'hypothesis_tests.txt'}\n")

# #==============================================================================
# # SECCIÓN 4: VISUALIZACIONES
# #==============================================================================

#     def create_all_plots(self):
#         """Genera todas las visualizaciones"""
#         print("="*70)
#         print("GENERACIÓN DE VISUALIZACIONES")
#         print("="*70 + "\n")
        
#         # 1. Efectos principales
#         self.plot_main_effects()
        
#         # 2. Interacciones
#         self.plot_interactions()
        
#         # 3. Polarización
#         self.plot_polarization_metrics()
        
#         # 4. Convergencia
#         self.plot_convergence_analysis()
        
#         # 5. Propiedades de red
#         self.plot_network_properties()
        
#         # 6. Heatmaps
#         self.plot_heatmaps()
        
#         # 7. Correlaciones
#         self.plot_correlation_matrix()
        
#         print("\n✓ Todas las visualizaciones generadas\n")
    
#     def plot_main_effects(self):
#         """Gráficos de efectos principales"""
#         print("📈 Generando gráficos de efectos principales...")
        
#         fig, axes = plt.subplots(2, 2, figsize=(16, 12))
#         fig.suptitle('Efectos Principales sobre Polarización (NDI)', fontsize=16, fontweight='bold')
        
#         # 1. Régimen de red
#         ax = axes[0, 0]
#         sns.boxplot(data=self.df, x='regimen_cat', y='mean_NDI', ax=ax, palette='Set2')
#         ax.set_title('A) Efecto del Régimen de Red')
#         ax.set_xlabel('Régimen de Conectividad')
#         ax.set_ylabel('Polarización (NDI)')
#         ax.grid(True, alpha=0.3)
        
#         # 2. Proporción de trolls
#         ax = axes[0, 1]
#         sns.violinplot(data=self.df, x='trolls_cat', y='mean_NDI', ax=ax, palette='Set3')
#         ax.set_title('B) Efecto de Proporción de Trolls')
#         ax.set_xlabel('Proporción de Trolls (%)')
#         ax.set_ylabel('Polarización (NDI)')
#         ax.grid(True, alpha=0.3)
        
#         # 3. Lambda (susceptibilidad)
#         ax = axes[1, 0]
#         for regime in ['desconectada', 'umbral', 'fuerte']:
#             data = self.df[self.df['regimen'] == regime]
#             ax.scatter(data['lam'], data['mean_NDI'], label=regime, alpha=0.6, s=50)
#         ax.set_title('C) Efecto de Susceptibilidad (λ)')
#         ax.set_xlabel('Lambda (λ)')
#         ax.set_ylabel('Polarización (NDI)')
#         ax.legend(title='Régimen')
#         ax.grid(True, alpha=0.3)
        
#         # 4. Ubicación de trolls
#         ax = axes[1, 1]
#         data_with_trolls = self.df[self.df['trolls'] > 0]
#         if len(data_with_trolls) > 0:
#             sns.boxplot(data=data_with_trolls, x='loc_cat', y='mean_NDI', ax=ax)
#             ax.set_title('D) Efecto de Ubicación de Trolls')
#             ax.set_xlabel('Banda de Centralidad')
#             ax.set_ylabel('Polarización (NDI)')
#             ax.grid(True, alpha=0.3)
        
#         plt.tight_layout()
#         plt.savefig(self.output_dir / 'main_effects.png', dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"  ✓ Guardado: main_effects.png")
    
#     def plot_interactions(self):
#         """Gráficos de interacciones"""
#         print("📈 Generando gráficos de interacciones...")
        
#         fig, axes = plt.subplots(1, 2, figsize=(16, 6))
#         fig.suptitle('Interacciones entre Factores', fontsize=16, fontweight='bold')
        
#         # 1. Régimen x Trolls
#         ax = axes[0]
#         for regime in ['desconectada', 'umbral', 'fuerte']:
#             data = self.df[self.df['regimen'] == regime]
#             means = data.groupby('trolls_pct')['mean_NDI'].mean()
#             ax.plot(means.index, means.values, marker='o', linewidth=2, markersize=8, label=regime)
        
#         ax.set_title('A) Régimen × Proporción de Trolls')
#         ax.set_xlabel('Proporción de Trolls (%)')
#         ax.set_ylabel('Polarización Media (NDI)')
#         ax.legend(title='Régimen')
#         ax.grid(True, alpha=0.3)
        
#         # 2. Lambda x Trolls
#         ax = axes[1]
#         for trolls_cat in ['0%', '10%', '20%', '30%']:
#             data = self.df[self.df['trolls_cat'] == trolls_cat]
#             if len(data) > 0:
#                 means = data.groupby('lam')['mean_NDI'].mean()
#                 ax.plot(means.index, means.values, marker='s', linewidth=2, markersize=8, label=trolls_cat)
        
#         ax.set_title('B) Lambda × Proporción de Trolls')
#         ax.set_xlabel('Lambda (λ)')
#         ax.set_ylabel('Polarización Media (NDI)')
#         ax.legend(title='Trolls')
#         ax.grid(True, alpha=0.3)
        
#         plt.tight_layout()
#         plt.savefig(self.output_dir / 'interactions.png', dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"  ✓ Guardado: interactions.png")
    
#     def plot_polarization_metrics(self):
#         """Comparación de métricas de polarización"""
#         print("📈 Generando gráficos de métricas de polarización...")
        
#         fig, axes = plt.subplots(2, 2, figsize=(16, 12))
#         fig.suptitle('Métricas de Polarización', fontsize=16, fontweight='bold')
        
#         # NDI por régimen y trolls
#         ax = axes[0, 0]
#         pivot_ndi = self.df.pivot_table(values='mean_NDI', index='regimen_cat', columns='trolls_cat')
#         sns.heatmap(pivot_ndi, annot=True, fmt='.3f', cmap='YlOrRd', ax=ax, cbar_kws={'label': 'NDI'})
#         ax.set_title('A) NDI: Network Disagreement Index')
#         ax.set_xlabel('Proporción de Trolls')
#         ax.set_ylabel('Régimen')
        
#         # P2 por régimen y trolls
#         ax = axes[0, 1]
#         pivot_p2 = self.df.pivot_table(values='mean_P2', index='regimen_cat', columns='trolls_cat')
#         sns.heatmap(pivot_p2, annot=True, fmt='.3f', cmap='RdPu', ax=ax, cbar_kws={'label': 'P2'})
#         ax.set_title('B) P2: Energía Media')
#         ax.set_xlabel('Proporción de Trolls')
#         ax.set_ylabel('Régimen')
        
#         # P4 por régimen y trolls
#         ax = axes[1, 0]
#         pivot_p4 = self.df.pivot_table(values='mean_P4', index='regimen_cat', columns='trolls_cat')
#         sns.heatmap(pivot_p4, annot=True, fmt='.1f', cmap='Blues', ax=ax, cbar_kws={'label': 'P4'})
#         ax.set_title('C) P4: Extremismo Absoluto')
#         ax.set_xlabel('Proporción de Trolls')
#         ax.set_ylabel('Régimen')
        
#         # Comparación de métricas normalizadas
#         ax = axes[1, 1]
#         metrics = ['mean_NDI', 'mean_P2', 'mean_P4']
#         df_norm = self.df[metrics].apply(lambda x: (x - x.min()) / (x.max() - x.min()))
#         df_norm.columns = ['NDI (norm)', 'P2 (norm)', 'P4 (norm)']
#         df_norm.boxplot(ax=ax, patch_artist=True)
#         ax.set_title('D) Comparación de Métricas (Normalizadas)')
#         ax.set_ylabel('Valor Normalizado [0,1]')
#         ax.grid(True, alpha=0.3)
        
#         plt.tight_layout()
#         plt.savefig(self.output_dir / 'polarization_metrics.png', dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"  ✓ Guardado: polarization_metrics.png")
    
#     def plot_convergence_analysis(self):
#         """Análisis de convergencia"""
#         print("📈 Generando gráficos de convergencia...")
        
#         fig, axes = plt.subplots(2, 2, figsize=(16, 12))
#         fig.suptitle('Análisis de Convergencia', fontsize=16, fontweight='bold')
        
#         # 1. Tiempo de convergencia vs lambda
#         ax = axes[0, 0]
#         for regime in ['desconectada', 'umbral', 'fuerte']:
#             data = self.df[self.df['regimen'] == regime]
#             ax.scatter(data['lam'], data['mean_convTime'], label=regime, alpha=0.6, s=50)
#         ax.set_title('A) Tiempo de Convergencia vs Lambda')
#         ax.set_xlabel('Lambda (λ)')
#         ax.set_ylabel('Tiempo Medio de Convergencia')
#         ax.legend(title='Régimen')
#         ax.grid(True, alpha=0.3)
        
#         # 2. Radio espectral vs convergencia
#         ax = axes[0, 1]
#         ax.scatter(self.df['mean_rho_LW'], self.df['mean_convTime'], alpha=0.6, s=50, c=self.df['trolls_pct'], cmap='viridis')
#         ax.set_title('B) Radio Espectral vs Convergencia')
#         ax.set_xlabel('Radio Espectral ρ(ΛW)')
#         ax.set_ylabel('Tiempo Medio de Convergencia')
#         cbar = plt.colorbar(ax.collections[0], ax=ax)
#         cbar.set_label('% Trolls')
#         ax.grid(True, alpha=0.3)
        
#         # 3. Proporción de consenso por régimen
#         ax = axes[1, 0]
#         consensus_data = self.df.groupby(['regimen_cat', 'trolls_cat'])['prop_consenso'].mean().unstack()
#         consensus_data.plot(kind='bar', ax=ax, width=0.8, colormap='tab10')
#         ax.set_title('C) Proporción de Consenso')
#         ax.set_xlabel('Régimen')
#         ax.set_ylabel('Proporción de Consenso')
#         ax.legend(title='% Trolls', bbox_to_anchor=(1.05, 1))
#         ax.grid(True, alpha=0.3, axis='y')
#         plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
        
#         # 4. Distribución de tiempo de convergencia
#         ax = axes[1, 1]
#         for regime in ['desconectada', 'umbral', 'fuerte']:
#             data = self.df[self.df['regimen'] == regime]['mean_convTime'].dropna()
#             if len(data) > 0:
#                 ax.hist(data, alpha=0.6, bins=20, label=regime, edgecolor='black')
#         ax.set_title('D) Distribución de Tiempo de Convergencia')
#         ax.set_xlabel('Tiempo de Convergencia')
#         ax.set_ylabel('Frecuencia')
#         ax.legend(title='Régimen')
#         ax.grid(True, alpha=0.3, axis='y')
        
#         plt.tight_layout()
#         plt.savefig(self.output_dir / 'convergence_analysis.png', dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"  ✓ Guardado: convergence_analysis.png")
    
#     def plot_network_properties(self):
#         """Propiedades de las redes"""
#         print("📈 Generando gráficos de propiedades de red...")
        
#         fig, axes = plt.subplots(2, 3, figsize=(18, 12))
#         fig.suptitle('Propiedades de las Redes', fontsize=16, fontweight='bold')
        
#         # 1. Componentes fuertemente conexas (W)
#         ax = axes[0, 0]
#         sns.boxplot(data=self.df, x='regimen_cat', y='mean_nSCC_W', ax=ax, palette='Set2')
#         ax.set_title('A) SCCs en W')
#         ax.set_xlabel('Régimen')
#         ax.set_ylabel('Número de SCCs')
#         ax.grid(True, alpha=0.3, axis='y')
        
#         # 2. Componentes fuertemente conexas (ΛW)
#         ax = axes[0, 1]
#         sns.boxplot(data=self.df, x='regimen_cat', y='mean_nSCC_LW', ax=ax, palette='Set3')
#         ax.set_title('B) SCCs en ΛW')
#         ax.set_xlabel('Régimen')
#         ax.set_ylabel('Número de SCCs')
#         ax.grid(True, alpha=0.3, axis='y')
        
#         # 3. Proporción fuertemente conexa
#         ax = axes[0, 2]
#         strong_data = self.df.groupby('regimen_cat')[['prop_isStrongW', 'prop_isStrongLW']].mean()
#         strong_data.plot(kind='bar', ax=ax, width=0.8)
#         ax.set_title('C) Proporción Fuertemente Conexa')
#         ax.set_xlabel('Régimen')
#         ax.set_ylabel('Proporción')
#         ax.legend(['W', 'ΛW'])
#         ax.grid(True, alpha=0.3, axis='y')
#         plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
        
#         # 4. Radio espectral W vs ΛW
#         ax = axes[1, 0]
#         ax.scatter(self.df['mean_rho_W'], self.df['mean_rho_LW'], 
#                    c=self.df['trolls_pct'], cmap='plasma', alpha=0.6, s=50)
#         ax.plot([0, 1], [0, 1], 'k--', alpha=0.3, label='y=x')
#         ax.set_title('D) Radio Espectral: W vs ΛW')
#         ax.set_xlabel('ρ(W)')
#         ax.set_ylabel('ρ(ΛW)')
#         ax.legend()
#         ax.grid(True, alpha=0.3)
#         cbar = plt.colorbar(ax.collections[0], ax=ax)
#         cbar.set_label('% Trolls')
        
#         # 5. Período por régimen
#         ax = axes[1, 1]
#         period_data = self.df.groupby('regimen_cat')[['mean_periodW', 'mean_periodLW']].mean()
#         period_data.plot(kind='bar', ax=ax, width=0.8)
#         ax.set_title('E) Período Medio de Grafos')
#         ax.set_xlabel('Régimen')
#         ax.set_ylabel('Período')
#         ax.legend(['W', 'ΛW'])
#         ax.grid(True, alpha=0.3, axis='y')
#         plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
        
#         # 6. Proporción primitiva (ΛW)
#         ax = axes[1, 2]
#         prim_data = self.df.groupby(['regimen_cat', 'trolls_cat'])['prop_isPrimLW'].mean().unstack()
#         prim_data.plot(kind='bar', ax=ax, width=0.8, colormap='Spectral')
#         ax.set_title('F) Proporción de Grafos Primitivos (ΛW)')
#         ax.set_xlabel('Régimen')
#         ax.set_ylabel('Proporción Primitiva')
#         ax.legend(title='% Trolls', bbox_to_anchor=(1.05, 1))
#         ax.grid(True, alpha=0.3, axis='y')
#         plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
        
#         plt.tight_layout()
#         plt.savefig(self.output_dir / 'network_properties.png', dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"  ✓ Guardado: network_properties.png")
    
#     def plot_heatmaps(self):
#         """Heatmaps de resultados clave"""
#         print("📈 Generando heatmaps...")
        
#         fig, axes = plt.subplots(2, 2, figsize=(16, 14))
#         fig.suptitle('Heatmaps: Régimen × Trolls', fontsize=16, fontweight='bold')
        
#         # 1. Rango final
#         ax = axes[0, 0]
#         pivot = self.df.pivot_table(values='mean_rangoFinal', 
#                                      index='regimen_cat', 
#                                      columns='trolls_cat')
#         sns.heatmap(pivot, annot=True, fmt='.3f', cmap='RdYlGn_r', ax=ax,
#                     cbar_kws={'label': 'Rango Final'})
#         ax.set_title('A) Rango Final de Opiniones')
#         ax.set_xlabel('Proporción de Trolls')
#         ax.set_ylabel('Régimen')
        
#         # 2. Desviación estándar final
#         ax = axes[0, 1]
#         pivot = self.df.pivot_table(values='mean_stdFinal', 
#                                      index='regimen_cat', 
#                                      columns='trolls_cat')
#         sns.heatmap(pivot, annot=True, fmt='.3f', cmap='YlOrBr', ax=ax,
#                     cbar_kws={'label': 'Std Final'})
#         ax.set_title('B) Desviación Estándar Final')
#         ax.set_xlabel('Proporción de Trolls')
#         ax.set_ylabel('Régimen')
        
#         # 3. Probabilidad frontera normales-trolls
#         ax = axes[1, 0]
#         data_with_trolls = self.df[self.df['trolls'] > 0]
#         if len(data_with_trolls) > 0:
#             pivot = data_with_trolls.pivot_table(values='mean_P_norm_vs_trolls',
#                                                  index='regimen_cat',
#                                                  columns='trolls_cat')
#             sns.heatmap(pivot, annot=True, fmt='.3f', cmap='coolwarm', ax=ax,
#                         cbar_kws={'label': 'P(N→T)'})
#             ax.set_title('C) Probabilidad Frontera Normales→Trolls')
#             ax.set_xlabel('Proporción de Trolls')
#             ax.set_ylabel('Régimen')
        
#         # 4. Radio espectral ΛW
#         ax = axes[1, 1]
#         pivot = self.df.pivot_table(values='mean_rho_LW',
#                                      index='regimen_cat',
#                                      columns='trolls_cat')
#         sns.heatmap(pivot, annot=True, fmt='.3f', cmap='viridis', ax=ax,
#                     cbar_kws={'label': 'ρ(ΛW)'})
#         ax.set_title('D) Radio Espectral de ΛW')
#         ax.set_xlabel('Proporción de Trolls')
#         ax.set_ylabel('Régimen')
        
#         plt.tight_layout()
#         plt.savefig(self.output_dir / 'heatmaps.png', dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"  ✓ Guardado: heatmaps.png")
    
#     def plot_correlation_matrix(self):
#         """Matriz de correlaciones"""
#         print("📈 Generando matriz de correlaciones...")
        
#         # Seleccionar variables numéricas clave
#         vars_of_interest = [
#             'mean_rangoFinal', 'mean_stdFinal', 'mean_NDI', 'mean_P2', 'mean_P4',
#             'mean_convTime', 'mean_rho_LW', 'prop_consenso',
#             'fracTrolls', 'lam', 'p'
#         ]
        
#         # Filtrar variables que existen
#         vars_available = [v for v in vars_of_interest if v in self.df.columns]
        
#         corr_data = self.df[vars_available].corr()
        
#         fig, ax = plt.subplots(figsize=(12, 10))
        
#         # Crear máscara para triángulo superior
#         mask = np.triu(np.ones_like(corr_data, dtype=bool))
        
#         sns.heatmap(corr_data, mask=mask, annot=True, fmt='.2f', 
#                     cmap='coolwarm', center=0, square=True, ax=ax,
#                     linewidths=1, cbar_kws={'label': 'Correlación de Pearson'})
        
#         ax.set_title('Matriz de Correlaciones entre Variables Clave', 
#                      fontsize=14, fontweight='bold', pad=20)
        
#         # Rotar etiquetas para mejor lectura
#         plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
#         plt.setp(ax.get_yticklabels(), rotation=0)
        
#         plt.tight_layout()
#         plt.savefig(self.output_dir / 'correlation_matrix.png', dpi=300, bbox_inches='tight')
#         plt.close()
        
#         print(f"  ✓ Guardado: correlation_matrix.png")

#==============================================================================
# SECCIÓN 5: ANÁLISIS COMPARATIVO POR LAMBDA
#==============================================================================

    def lambda_analysis(self):
        """Análisis detallado del efecto de lambda"""
        print("="*70)
        print("ANÁLISIS DETALLADO DE LAMBDA (SUSCEPTIBILIDAD)")
        print("="*70 + "\n")
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Efecto de Lambda (Susceptibilidad) en las Dinámicas', 
                     fontsize=16, fontweight='bold')
        
        # 1. Lambda vs Polarización por régimen
        ax = axes[0, 0]
        for regime in ['desconectada', 'umbral', 'fuerte']:
            data = self.df[self.df['regimen'] == regime]
            lambda_groups = data.groupby('lam')['mean_NDI'].agg(['mean', 'std'])
            ax.errorbar(lambda_groups.index, lambda_groups['mean'], 
                       yerr=lambda_groups['std'], marker='o', capsize=5, 
                       linewidth=2, markersize=8, label=regime)
        
        ax.set_title('A) Lambda vs Polarización (NDI)')
        ax.set_xlabel('Lambda (λ)')
        ax.set_ylabel('NDI (media ± std)')
        ax.legend(title='Régimen')
        ax.grid(True, alpha=0.3)
        
        # 2. Lambda vs Convergencia
        ax = axes[0, 1]
        for regime in ['desconectada', 'umbral', 'fuerte']:
            data = self.df[self.df['regimen'] == regime]
            lambda_groups = data.groupby('lam')['mean_convTime'].agg(['mean', 'std'])
            ax.errorbar(lambda_groups.index, lambda_groups['mean'],
                       yerr=lambda_groups['std'], marker='s', capsize=5,
                       linewidth=2, markersize=8, label=regime)
        
        ax.set_title('B) Lambda vs Tiempo de Convergencia')
        ax.set_xlabel('Lambda (λ)')
        ax.set_ylabel('Tiempo (media ± std)')
        ax.legend(title='Régimen')
        ax.grid(True, alpha=0.3)
        
        # 3. Lambda vs Consenso
        ax = axes[1, 0]
        for trolls_cat in ['0%', '10%', '20%', '30%']:
            data = self.df[self.df['trolls_cat'] == trolls_cat]
            if len(data) > 0:
                lambda_groups = data.groupby('lam')['prop_consenso'].mean()
                ax.plot(lambda_groups.index, lambda_groups.values,
                       marker='D', linewidth=2, markersize=8, label=trolls_cat)
        
        ax.set_title('C) Lambda vs Proporción de Consenso')
        ax.set_xlabel('Lambda (λ)')
        ax.set_ylabel('Proporción de Consenso')
        ax.legend(title='% Trolls')
        ax.grid(True, alpha=0.3)
        
        # 4. Lambda vs Rango Final
        ax = axes[1, 1]
        scatter = ax.scatter(self.df['lam'], self.df['mean_rangoFinal'],
                            c=self.df['trolls_pct'], cmap='plasma',
                            alpha=0.6, s=50, edgecolor='black', linewidth=0.5)
        ax.set_title('D) Lambda vs Rango Final')
        ax.set_xlabel('Lambda (λ)')
        ax.set_ylabel('Rango Final de Opiniones')
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('% Trolls')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'lambda_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Análisis de lambda completado: lambda_analysis.png\n")

#==============================================================================
# SECCIÓN 6: ANÁLISIS DE UBICACIÓN DE TROLLS
#==============================================================================

    def location_analysis(self):
        """Análisis del efecto de la ubicación de trolls"""
        print("="*70)
        print("ANÁLISIS DE UBICACIÓN DE TROLLS (CENTRALIDAD)")
        print("="*70 + "\n")
        
        data_with_trolls = self.df[self.df['trolls'] > 0]
        
        if len(data_with_trolls) == 0:
            print("⚠️  No hay datos con trolls para analizar ubicación\n")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Efecto de la Ubicación de Trolls (Centralidad)', 
                     fontsize=16, fontweight='bold')
        
        # 1. Ubicación vs Polarización por régimen
        ax = axes[0, 0]
        for regime in ['desconectada', 'umbral', 'fuerte']:
            data = data_with_trolls[data_with_trolls['regimen'] == regime]
            if len(data) > 0:
                loc_groups = data.groupby('loc_cat')['mean_NDI'].mean()
                ax.plot(loc_groups.index.astype(str), loc_groups.values,
                       marker='o', linewidth=2, markersize=10, label=regime)
        
        ax.set_title('A) Ubicación vs Polarización (NDI)')
        ax.set_xlabel('Banda de Centralidad')
        ax.set_ylabel('NDI Medio')
        ax.legend(title='Régimen')
        ax.grid(True, alpha=0.3)
        
        # 2. Ubicación vs Rango Final
        ax = axes[0, 1]
        sns.boxplot(data=data_with_trolls, x='loc_cat', y='mean_rangoFinal', ax=ax, palette='Set2')
        ax.set_title('B) Ubicación vs Rango Final')
        ax.set_xlabel('Banda de Centralidad')
        ax.set_ylabel('Rango Final')
        ax.legend(title='Régimen', bbox_to_anchor=(1.05, 1))
        ax.grid(True, alpha=0.3, axis='y')
        
        # 3. Heatmap: Ubicación x % Trolls
        ax = axes[1, 0]
        pivot = data_with_trolls.pivot_table(values='mean_NDI',
                                             index='loc_cat',
                                             columns='trolls_cat')
        sns.heatmap(pivot, annot=True, fmt='.3f', cmap='YlOrRd', ax=ax,
                    cbar_kws={'label': 'NDI'})
        ax.set_title('C) NDI: Ubicación × % Trolls')
        ax.set_xlabel('Proporción de Trolls')
        ax.set_ylabel('Banda de Centralidad')
        
        # 4. Ubicación vs Consenso
        ax = axes[1, 1]
        loc_consensus = data_with_trolls.groupby(['loc_cat', 'trolls_cat'])['prop_consenso'].mean().unstack()
        loc_consensus.plot(kind='bar', ax=ax, width=0.8, colormap='tab10')
        ax.set_title('D) Ubicación vs Proporción de Consenso')
        ax.set_xlabel('Banda de Centralidad')
        ax.set_ylabel('Proporción de Consenso')
        ax.legend(title='% Trolls', bbox_to_anchor=(1.05, 1))
        ax.grid(True, alpha=0.3, axis='y')
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'location_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Análisis de ubicación completado: location_analysis.png\n")

#==============================================================================
# SECCIÓN 7: INFORME EJECUTIVO
#==============================================================================

    def generate_executive_summary(self, hypothesis_results):
        """Genera informe ejecutivo en texto"""
        print("="*70)
        print("GENERANDO INFORME EJECUTIVO")
        print("="*70 + "\n")
        
        report_path = self.output_dir / 'executive_summary.txt'
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("="*70 + "\n")
            f.write("INFORME EJECUTIVO: EXPERIMENTOS FRIEDKIN-JOHNSEN\n")
            f.write("="*70 + "\n\n")
            
            # 1. Resumen del experimento
            f.write("1. DISEÑO EXPERIMENTAL\n")
            f.write("-" * 70 + "\n")
            f.write(f"   • Total de escenarios: {len(self.df)}\n")
            f.write(f"   • Tamaño de red: n = {self.df['n'].iloc[0]}\n")
            f.write(f"   • Regímenes de red: {', '.join(self.df['regimen'].unique())}\n")
            f.write(f"   • Proporciones de trolls: {sorted(self.df['trolls_pct'].unique())}%\n")
            f.write(f"   • Valores de lambda: {sorted(self.df['lam'].unique())}\n")
            f.write(f"   • Ubicaciones de trolls: {', '.join(self.df['loc'].unique())}\n\n")
            
            # 2. Hallazgos principales
            f.write("2. HALLAZGOS PRINCIPALES\n")
            f.write("-" * 70 + "\n\n")
            
            # Efecto del régimen
            regime_effect = self.df.groupby('regimen')['mean_NDI'].mean()
            f.write(f"   A) EFECTO DEL RÉGIMEN DE CONECTIVIDAD\n")
            f.write(f"      • Desconectada: NDI = {regime_effect['desconectada']:.4f}\n")
            f.write(f"      • Umbral:       NDI = {regime_effect['umbral']:.4f}\n")
            f.write(f"      • Fuerte:       NDI = {regime_effect['fuerte']:.4f}\n")
            
            if hypothesis_results['H1_regime_polarization']['significant']:
                f.write(f"      ✓ Efecto SIGNIFICATIVO (p < 0.05)\n")
            else:
                f.write(f"      ✗ Efecto NO significativo\n")
            f.write("\n")
            
            # Efecto de trolls
            trolls_effect = self.df.groupby('trolls_pct')['mean_NDI'].mean()
            f.write(f"   B) EFECTO DE LA PROPORCIÓN DE TROLLS\n")
            for pct in sorted(trolls_effect.index):
                f.write(f"      • {pct:.0f}% trolls: NDI = {trolls_effect[pct]:.4f}\n")
            
            if hypothesis_results['H2_trolls_polarization']['significant']:
                rho = hypothesis_results['H2_trolls_polarization']['rho']
                f.write(f"      ✓ Correlación SIGNIFICATIVA (ρ = {rho:.3f}, p < 0.05)\n")
            else:
                f.write(f"      ✗ Correlación NO significativa\n")
            f.write("\n")
            
            # Efecto de lambda
            lambda_effect = self.df.groupby('lam')['mean_convTime'].mean()
            f.write(f"   C) EFECTO DE LAMBDA EN CONVERGENCIA\n")
            for lam in sorted(lambda_effect.index):
                f.write(f"      • λ = {lam:.2f}: Tiempo = {lambda_effect[lam]:.2f}\n")
            
            if hypothesis_results['H3_lambda_convergence']['significant']:
                r = hypothesis_results['H3_lambda_convergence']['r']
                f.write(f"      ✓ Correlación SIGNIFICATIVA (r = {r:.3f}, p < 0.05)\n")
            else:
                f.write(f"      ✗ Correlación NO significativa\n")
            f.write("\n")
            
            # Ubicación de trolls
            data_with_trolls = self.df[self.df['trolls'] > 0]
            if len(data_with_trolls) > 0:
                loc_effect = data_with_trolls.groupby('loc')['mean_rangoFinal'].mean()
                f.write(f"   D) EFECTO DE LA UBICACIÓN DE TROLLS\n")
                for loc in ['low', 'mid', 'high']:
                    if loc in loc_effect.index:
                        f.write(f"      • {loc.capitalize()}: Rango = {loc_effect[loc]:.4f}\n")
                
                if hypothesis_results['H4_location_influence']['significant']:
                    f.write(f"      ✓ Efecto SIGNIFICATIVO (p < 0.05)\n")
                else:
                    f.write(f"      ✗ Efecto NO significativo\n")
                f.write("\n")
            
            # 3. Métricas de polarización
            f.write("3. MÉTRICAS DE POLARIZACIÓN\n")
            f.write("-" * 70 + "\n")
            f.write(f"   • NDI medio:  {self.df['mean_NDI'].mean():.4f} ± {self.df['mean_NDI'].std():.4f}\n")
            f.write(f"   • P2 medio:   {self.df['mean_P2'].mean():.4f} ± {self.df['mean_P2'].std():.4f}\n")
            f.write(f"   • P4 medio:   {self.df['mean_P4'].mean():.2f} ± {self.df['mean_P4'].std():.2f}\n")
            f.write(f"   • Rango medio: {self.df['mean_rangoFinal'].mean():.4f} ± {self.df['mean_rangoFinal'].std():.4f}\n\n")
            
            # 4. Convergencia
            f.write("4. CONVERGENCIA\n")
            f.write("-" * 70 + "\n")
            f.write(f"   • Tiempo medio de convergencia: {self.df['mean_convTime'].mean():.2f} iteraciones\n")
            f.write(f"   • Proporción de consenso global: {self.df['prop_consenso'].mean():.2%}\n")
            f.write(f"   • Radio espectral ΛW medio: {self.df['mean_rho_LW'].mean():.4f}\n\n")
            
            # 5. Recomendaciones
            f.write("5. RECOMENDACIONES Y CONCLUSIONES\n")
            f.write("-" * 70 + "\n")
            
            # Análisis automático de qué factores son más importantes
            significant_factors = []
            if hypothesis_results['H1_regime_polarization']['significant']:
                significant_factors.append('régimen de red')
            if hypothesis_results['H2_trolls_polarization']['significant']:
                significant_factors.append('proporción de trolls')
            if hypothesis_results['H3_lambda_convergence']['significant']:
                significant_factors.append('lambda (susceptibilidad)')
            if hypothesis_results['H4_location_influence']['significant']:
                significant_factors.append('ubicación de trolls')
            
            if significant_factors:
                f.write(f"   • Factores significativos identificados: {', '.join(significant_factors)}\n\n")
            else:
                f.write(f"   • No se identificaron factores con efectos significativos fuertes.\n\n")
            
            # Recomendaciones específicas
            f.write("   Recomendaciones para mitigar polarización:\n")
            
            if hypothesis_results['H2_trolls_polarization']['significant']:
                if hypothesis_results['H2_trolls_polarization']['rho'] > 0:
                    f.write("   1. Reducir la proporción de trolls en la red (correlación positiva con polarización)\n")
            
            if hypothesis_results['H4_location_influence']['significant']:
                f.write("   2. Monitorizar trolls en posiciones de alta centralidad (mayor impacto)\n")
            
            if hypothesis_results['H3_lambda_convergence']['significant']:
                f.write("   3. Ajustar susceptibilidad (lambda) para optimizar convergencia\n")
            
            f.write("\n")
            f.write("="*70 + "\n")
            f.write("Fin del informe\n")
        
        print(f"✓ Informe ejecutivo guardado: {report_path}\n")

#==============================================================================
# FUNCIÓN PRINCIPAL
#==============================================================================

def main():
    """Función principal para ejecutar el análisis completo"""
    
    print("\n" + "="*70)
    print(" "*15 + "ANÁLISIS DE EXPERIMENTOS FRIEDKIN-JOHNSEN")
    print("="*70 + "\n")
    
    # Solicitar directorio de resultados
    import sys
    
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = input("Ingrese la ruta del directorio de resultados: ").strip()
    
    if not results_dir:
        print("❌ Error: Debe proporcionar un directorio de resultados")
        return
    
    try:
        # Inicializar analizador
        analyzer = FriedkinAnalyzer(results_dir)
        
        # 1. Análisis exploratorio
        desc_regime, desc_trolls, corr_matrix = analyzer.exploratory_analysis()
        
        # 2. Tests de hipótesis
        hypothesis_results = analyzer.hypothesis_tests()
        
        # 3. Visualizaciones
        #analyzer.create_all_plots()
        
        # 4. Análisis específicos
        analyzer.lambda_analysis()
        analyzer.location_analysis()
        
        # 5. Informe ejecutivo
        analyzer.generate_executive_summary(hypothesis_results)
        
        # Resumen final
        print("="*70)
        print("✓ ANÁLISIS COMPLETADO")
        print("="*70)
        print(f"\n📁 Todos los resultados guardados en: {analyzer.output_dir}")
        print("\nArchivos generados:")
        print("  • descriptive_stats.txt")
        print("  • hypothesis_tests.txt")
        # print("  • executive_summary.txt")
        # print("  • main_effects.png")
        # print("  • interactions.png")
        # print("  • polarization_metrics.png")
        # print("  • convergence_analysis.png")
        # print("  • network_properties.png")
        # print("  • heatmaps.png")
        # print("  • correlation_matrix.png")
        # print("  • lambda_analysis.png")
        # print("  • location_analysis.png")
        print("\n" + "="*70 + "\n")
        
    except FileNotFoundError as e:
        print(f"\n❌ Error: {e}")
        print("Verifique que la ruta del directorio sea correcta.\n")
    except Exception as e:
        print(f"\n❌ Error inesperado: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()