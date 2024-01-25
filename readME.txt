/*
    Authors :
    KAOUCH Abdelssamad
    AIT KHELIFA Tanina
*/

Pour compiler : make
Autres commande utile au besoin : make clean

Paramètre d'exécution :
./ver_double d n min max epsilon
./ver_mpfr d n min max epsilon prec

d : si 0, les coefficients tirés au hasard seront entiers
    (stockés dans des doubles/mpfr, ex : 5.0, 1.0, etc...), 0 par défaut
n : taille de la matrice initialisée, 4 par défaut
min : valeur minimale d'un coefficient de la matrice initialisée, 1 par défaut
max : valeur maximale d'un coefficient de la matrice , 10 par défaut
prec : mpfr only, précision des nombres à virgule flottante, 100 par défaut


Exemple de commande pour exécuter le programme :
./ver_double 0 4 1 15 0.001
./ver_mpfr 0 4 1 15 0.001 100