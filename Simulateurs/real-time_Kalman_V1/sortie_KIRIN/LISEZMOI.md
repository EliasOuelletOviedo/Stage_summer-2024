# Simulateur Kalman temps réel, étape 1 — capteur KIRIN

Course du 2026-10-19 : le même simulateur (V1, étape 1 et ses deux passes de corrections), alimenté par la titration
KIRIN du 2026-10-16 à 800 nm au lieu de la titration SCLM.

```bash
julia --project=../FLIM-FRET_2-channels step0_params.jl --titration KIRIN --params params_KIRIN.toml
julia --project=. -t auto run.jl --params params_KIRIN.toml --sorties sortie_KIRIN
```

## Ce qui change dans les paramètres

`params_KIRIN.toml` reprend tous les réglages de `params.toml` (trajectoire 5 → 30 mM, budgets, dérives, filtre de
référence OU, 5a avec diffusion), avec trois différences :

- le bloc « step0 » vient de l'ajustement global de KIRIN (`deux_canaux_2026-10-16_KIRIN`) ;
- la plage calibrée va de 0 à 150 mM (celle de la titration) ;
- les bornes de fraction de 5a passent à [−1,5 ; 1,5] : chez KIRIN, q vaut environ −0,5 (la composante rapide du
  jaune est une montée, sensibilisée), exactement sur la borne de −0,5 d'origine.

| Langmuir (8 points) | y₀ | y∞ | K |
| --- | --- | --- | --- |
| J/B | 1,57 | 3,82 | 3,6 mM |
| p (bleu) | 0,19 | 0,59 | 2,4 mM |
| q (jaune) | −0,53 | −0,48 | 6,0 mM |

J/B **monte** avec le Cl⁻ et sature tôt : la trajectoire 5 → 30 mM se déroule presque entièrement dans la zone
plate (à 25 mM, la pente de J/B n'est plus que 0,8 % par mM). q ne porte presque rien : la règle écrite d'avance
le laisse hors de l'observable de 5a (w = (1, 0)), qui est donc p seul.

## Résultats (8 répétitions, filtre de référence OU)

| N = 1000 | J/B : biais 30 s | 5a : biais 30 s | Kalman : biais 30 s | Kalman : RMSE |
| --- | --- | --- | --- | --- |
| aucune dérive | +0,54 | +4,3 | +0,29 ± 0,53 | 2,69 |
| rampe 10 % | **+9,40** | +4,3 | +0,54 ± 0,54 | 2,57 |
| rampe 20 % | **+35,8** | +4,3 | +0,82 ± 0,56 | 2,56 |
| saut 5 % | +4,76 | +4,3 | +0,11 ± 0,50 | 3,00 |

À N = 3000 : J/B +7,2 mM sous rampe de 10 % et +37 mM sous rampe de 20 % ; le Kalman +0,26 et +0,40 mM.

- **La dérive de gain est catastrophique pour J/B seul chez KIRIN.** Dans la zone saturée, une erreur de gain de
  quelques pour cent se lit comme des dizaines de mM : +9 mM pour une rampe de 10 %, +36 mM pour 20 % (où 24 % des
  lectures sortent même de la plage calibrée). Le biais est positif parce que J/B monte avec le Cl⁻.
- **La fusion tient.** Le Kalman reste sous 1 mM de biais dans toutes les dérives, avec une RMSE de 2,6 mM à
  N = 1000 (1,9 à 3000) contre 11,7 à 41 mM pour J/B seul. 5a, insensible à g, fait son travail d'ancre.
- **Le prix est la précision brute.** σ_C à 20 mM : 9,3 mM pour J/B, 29,5 mM pour 5a à N = 1000 (5,1 et 16,3 à
  3000). Le capteur est peu sensible dans cette plage : le Kalman ne descend pas sous 2 mM de RMSE.
- **Critères de sortie** (`criteres.csv`) : à N = 3000, le filtre OU les passe tous (biais rampe 10 % +0,26 mM,
  P₁₁/err² 0,89 à 1,53, RMSE au plus 0,33 × J/B). À N = 1000, le critère de biais tombe (+0,54 mM pour un seuil de
  0,3) ; les deux autres passent.
- **Fisher contre Monte-Carlo** : variances à 0,2 % (J/B) et 1,2 % (5a) à N = 1000.
- **Amplitudes** : sans bruit, 5a retrouve p et q injectés à 10⁻¹² près, avec et sans diffusion.

## Réserves

- La trajectoire n'a pas été adaptée à KIRIN : elle reste celle de l'étape 1 (repos à 5 mM, sommet à 30 mM), là
  où KIRIN est déjà saturé. Une trajectoire autour de son K (2 à 4 mM) donnerait une autre image.
- `seuils_relance.csv` applique les seuils fixés pour la deuxième passe des corrections, dont les valeurs
  « avant » et les bornes (biais d'inversion dans [2,5 ; 3,1] mM, hausse de σ ≤ 10 %…) avaient été posées pour un
  autre capteur. Ils ne sont pas interprétables ici ; seuls les tests qui ne dépendent pas du capteur (amplitudes
  exactes, accord Fisher) ont un sens.
- Le filtre garde les réglages de Q choisis sur l'autre capteur (OU, σ_w,C = 1,2 mM par pas). Ils n'ont pas été
  revus pour KIRIN.
