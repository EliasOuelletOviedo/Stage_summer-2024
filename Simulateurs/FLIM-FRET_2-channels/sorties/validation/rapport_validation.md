# Validation interne du simulateur

Généré par `step2_validation.jl`. Modèle : τ_D0=3.689 τ_DA=0.151 τ_A=3.00 ns  x₀=0.429 x_sat=0.049 K_D=30.3 mM n=1.46  Φ_D=0.87 Φ_A=0.60 (effectif 0.60)  η=[0.336 0.0019; 0.164 0.590] ε=0.030 εΦ_dark=0.000  Δt_jaune=0 ps  β_A=0.00 β_D=0.00  M2_lent

| Résultat | Test | Critère | Valeur |
| --- | --- | --- | --- |
| OK | 1a aller-retour sans bruit, modèle complet (7 paramètres libres) | écart relatif max < 1e-8 | 8.9e-16 (8 itérations) |
| OK | 1b estimateur 4 (canal bleu sans fuite d'accepteur) sans bruit rend x | écart sur x < 1e-10 | 2.6e-15 |
| OK | 1b estimateur 5a sans bruit rend x | écart sur x < 1e-10 | 5.6e-16 |
| OK | 1b estimateur 5d sans bruit rend x | écart sur x < 1e-10 | 1.1e-16 |
| OK | 2  conservation : donneur Φ_D[x τ_DA/τ_D0 + 1 − x] | écart relatif < 1e-12 | 2.2e-16 |
| OK | 2  conservation : émission sensibilisée Φ_A·x·E | écart relatif < 1e-12 | 1.1e-16 |
| OK | 3  amplitudes du bleu idéal : x et 1 − x, rien sur τ_A | écart < 1e-14 | 0.0e+00 ; α(τ_A) = 0.0e+00 |
| OK | 4  Cl⁻ saturant : bleu mono-exponentiel à τ_D0 | seule l'amplitude τ_D0 non nulle | oui |
| OK | 4  Cl⁻ saturant : jaune = (η_YD/η_BD)·bleu, nul si η_YD = 0 | écart < 1e-14 | 8.7e-17 ; plus grande amplitude du jaune (η_YD=0) = 0.0e+00 |
| OK | 4  ton ajustement mono sur ce déclin rend τ_D0 | écart sur τ < 0,01 ns | 3.6910 contre 3.6886 ns (écart 0.06 % ≈ rapport des pas 0,048861/0,048828) |
| OK | 5  jaune pur sensibilisé : −A à τ_DA, +A à τ_A, rapport −1 | écart relatif < 1e-8 | τ = (0.151496, 3.000000) ns, A_court/A_long = -1.0000000000 ; écart 4.4e-16 |
| OK | 6  efficacité à 10⁶ photons, 5a conjoint | IC 95 % de CRLB/var recoupe [0,98 ; 1,02] | 1.005 [0.985 ; 1.024], biais sur x -8.4e-06 (20000 réal., 1 s) |
| OK | 6  efficacité à 10⁶ photons, 5d conjoint échelle commune | IC 95 % de CRLB/var recoupe [0,98 ; 1,02] | 0.999 [0.979 ; 1.018], biais sur x -8.7e-07 (20000 réal., 1 s) |
| OK | 6  efficacité à 10⁶ photons, 4 biexp τ figés (bleu sans fuite d'accepteur) | IC 95 % de CRLB/var recoupe [0,98 ; 1,02] | 0.996 [0.977 ; 1.016], biais sur x -3.2e-05 (20000 réal., 0 s) |
| OK | 6  efficacité à 10⁶ photons, 3 biexp τ libres (idem) | IC 95 % de CRLB/var recoupe [0,98 ; 1,02] | 0.995 [0.975 ; 1.014], biais sur x 1.1e-04 (20000 réal., 2 s) |
| OK | 7  dérivées analytiques contre ForwardDiff (6 jeux aléatoires) | écart relatif < 1e-10 | 2.1e-14 |
| OK | 8  phasor : la titration M2 lent est un segment (bleu et jaune) | distance max < 1e-10 (en longueur de segment) | bleu 6.9e-14, jaune 9.9e-15 ; M2 rapide : bleu 1.4e-01, jaune 1.5e-01 |
| OK | 9  noyau contre l'intégrale double exacte (quadrature) | écart relatif < 1e-11 | 8.7e-14 |
| OK | 10 convolution circulaire FFT contre convolution directe | écart relatif < 1e-12 | 6.5e-16 |
| OK | 11 adaptateur : ton τ mono de SCLM 0 mM c1 (raw_data.csv) | identique | 2.0707331135068343 contre 2.0707331135068343 |
| OK | 12 modèle à N états (canaux_etats) ≡ M2 lent à deux états | écart relatif < 1e-12 | 5.1e-16 |
| OK | 13 IRF dilatée : identité à s = 1, aire conservée, largeur monotone | 1e-12 ; croissante | 3.9e-16 ; aires −1 = 2.2e-16 ; largeurs 2.40 < 2.63 < 2.90 canaux |
| OK | 14 retard Δt du jaune : forme décalée et ∂/∂t₀ (différences finies) | 1e-12 ; 1e-5 | 0.0e+00 ; 4.1e-06 |
| OK | 15 blanchiment : β = 0 ≡ M2, β → 0 continu, populations conservées | 1e-12 ; 1e-7 ; 1e-12 | 3.3e-16 ; 3.9e-10 ; 1.1e-16 |
| OK | 16 blanchiment du donneur : moyenne de la survie sur l'acquisition | écart relatif < 1e-12 | 1.1e-16 |
| OK | 17 ingestion : titration synthétique sans bruit, paramètres retrouvés | τ < 1e-6 rel., x < 1e-6, t₀ < 1e-6 ns (aussi t₀ par histogramme) | τ 2.0e-09 ; x 4.6e-10 ; t₀ 2.4e-11 ns ; t₀ par histo 4.1e-10 ns |
| OK | 18 canal jaune : ajustement conjoint d'un banc simulé sans bruit | écart relatif max < 1e-4 | 1.4e-09 |
| OK | 19 dérive de gain : écart-type de log g après τ_g | à 5 % près (4 000 marches) | 0.049 contre 0.050 |
| OK | 20 TAC plus long que la période : convolution linéaire, dérivées | 1e-12 ; 1e-6 | 6.3e-16 ; 1.8e-09 (L = 512) |
