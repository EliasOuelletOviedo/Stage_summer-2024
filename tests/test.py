import pyabf
import matplotlib.pyplot as plt
import numpy as np

# Replace 'your_file.abf' with the actual path to your ABF file
abf1 = pyabf.ABF("/Users/eliasouellet-oviedo/Desktop/home/2025_08_18_0016.abf")
abf1.setSweep(sweepNumber=0)

abf2 = pyabf.ABF("/Users/eliasouellet-oviedo/Desktop/home/2025_08_18_0017.abf")
abf2.setSweep(sweepNumber=0)

moynnage = 500

y1 = abf1.sweepY
y2 = abf2.sweepY

x1 = abf1.sweepX
x2 = abf2.sweepX

# y1_m = np.zeros(len(y1) - moynnage +1)

# for i in range(len(y1)- moynnage +1):
#     if i < len(y1) - moynnage:
#         y1_m[i] = np.mean(y1[i : i + moynnage])
#     else:
#         y1_m[i] = np.mean(y1[i : len(y1)])

y1_m = []
x1_m = []
y2_m = []
x2_m = []

for i in range(0, len(y1), moynnage):
    # Crée une sous-liste de 5 éléments
    sous_liste = y1[i : i + moynnage]
    sous_liste2 = y2[i : i + moynnage]
    
    # Calcule la moyenne de cette sous-liste
    moyenne_sous_liste = sum(sous_liste) / len(sous_liste)
    moyenne_sous_liste2 = sum(sous_liste2) / len(sous_liste2)
    
    # Ajoute la moyenne à la liste des résultats
    y1_m.append(moyenne_sous_liste)
    x1_m.append(x1[i])
    y2_m.append(moyenne_sous_liste2)
    x2_m.append(x2[i])

fig, ax = plt.subplots(2, 1, figsize=(8, 6))
ax[0].plot(x1_m, y1_m)
ax[0].set_ylim(-75, 150)
ax[1].plot(x2_m, y2_m)
ax[1].set_ylim(-50, 200)

plt.show()