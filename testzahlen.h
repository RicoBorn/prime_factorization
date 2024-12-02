// Testzahlen.h
/*
 Funktionen zur Erzeugung der Testzahlen zum Praktikum Algebra 2022
 */
#ifndef TESTZAHLEN_H // Schutz vor doppeltem Laden
#define TESTZAHLEN_H
#include <gmpxx.h>

/*Funktionen zur Erzeugung von Testzahlen*/
/*Fermat(k) gibt die k-te Fermat-Zahl zurück*/
unsigned long int aux_fast_exp(int n, int p); // Berechnet p^n (Achtung: expoldiert schnell)
mpz_class Fermat(unsigned int k);
/*Cunningham(k) gibt die k-stellige Zahl 111...1111 (cf. Cunningham Project)*/
mpz_class Cunningham(int k);
/* Einfache Testzahlen B1,...,B10 zur Verifikation */
mpz_class TestzahlB(int k);
/* RSA-Zahlen R... */
mpz_class RSAZahl(int k);
#endif /*  TESTZAHLEN_H */
