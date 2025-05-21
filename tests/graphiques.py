import matplotlib.pyplot as plt
import csv
import numpy as np

def plot_graphique(df, x_col, y_col, title, xlabel, ylabel):
    plt.figure(figsize=(10, 6))
    plt.plot(df[x_col], df[y_col], marker='o', linestyle='-', color='b')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()
    plt.show()

    return plt

fichier_csv = "2025-05-14T14:21:50.335.csv"

# Listes pour stocker les colonnes
time = []
number_of_photons = []
lifetime = []

# Lecture du fichier CSV
with open(fichier_csv, 'r') as f:
    lecteur_csv = csv.reader(f)
    next(lecteur_csv)  # sauter la première ligne avec l'horodatage
    next(lecteur_csv)  # sauter la ligne d'en-tête
    next(lecteur_csv)  # sauter la ligne d'en-tête

    for ligne in lecteur_csv:
        if len(ligne) < 3:
            continue  # ignorer les lignes incomplètes
        time.append(float(ligne[0]))
        number_of_photons.append(float(ligne[1]))
        lifetime.append(float(ligne[2]))

averaged_time = np.zeros(len(time))

for i in range(len(time) - 5):
    averaged_time[i] = np.mean(lifetime[i:i+5])  # Normaliser le temps

# Tracer Lifetime vs Time
plt.figure(figsize=(10, 5))
# plt.plot(time[1000:], lifetime[1000:], marker=None, linestyle='-', color='black', label='lifetime', alpha=1, linewidth=1)
plt.plot(time[1000:-5], averaged_time[1000:-5], linewidth=1, label='Lifetime')
plt.title('Lifetime vs Time')
plt.xlabel('Time [s]')
plt.ylabel('Lifetime [ns]')
plt.grid(True)
plt.tight_layout()
plt.axvspan(1400, 1770, color='gold', alpha=0.3, label='Signal')
plt.axvspan(1900, 2150, color='gold', alpha=0.3)
plt.legend()
plt.show()
