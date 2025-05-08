import pandas as pd
import matplotlib.pyplot as plt

# Charger les données
data = pd.read_csv("webbase-1M_resultats.csv")

# Créer une figure avec deux sous-graphiques
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Graphique 1 : Itérations en fonction de Alpha pour chaque pourcentage de suppression
for suppression in data['Pourcentage_Suppression'].unique():
    subset = data[data['Pourcentage_Suppression'] == suppression]
    ax1.plot(subset['Alpha'], subset['Iterations'], marker='o', label=f'Suppression {suppression*100}%')

ax1.set_xlabel('Alpha')
ax1.set_ylabel('Nombre d\'itérations')
ax1.set_title('Itérations en fonction de Alpha')
ax1.legend()
ax1.grid(True)

# Graphique 2 : Temps en fonction de Alpha pour chaque pourcentage de suppression
for suppression in data['Pourcentage_Suppression'].unique():
    subset = data[data['Pourcentage_Suppression'] == suppression]
    ax2.plot(subset['Alpha'], subset['Temps'], marker='o', label=f'Suppression {suppression*100}%')

ax2.set_xlabel('Alpha')
ax2.set_ylabel('Temps de convergence (s)')
ax2.set_title('Temps en fonction de Alpha')
ax2.legend()
ax2.grid(True)

# Ajuster l'espacement
plt.tight_layout()

# Sauvegarder le graphique
plt.savefig('analyse_webbase-1M.png')

# Afficher le graphique (optionnel, commentez si vous ne voulez pas d'affichage)
# plt.show()