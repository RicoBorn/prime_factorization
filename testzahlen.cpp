/* Die Datei "testzahlen.cpp" enthält die Funktionen zur Erzeugung der Testzahlen zum Praktikum Algebra 2022 */
#include "testzahlen.h"
#include <gmpxx.h>
#include <string>

/* Hilfsfunktion für schnelle Potenzierung bei der Erzeugung der Fermatzahlen*/
unsigned long int aux_fast_exp(int n, int p)
{
    if(p == 0)
    {
        return 1;
    }
    if(p%2 == 0)
    {
            unsigned long int tmp = aux_fast_exp(n, p/2);
            return tmp*tmp;
    }
   else
    {
            unsigned long int tmp = aux_fast_exp(n, p-1);
            return tmp*n;
    }
}

/*Fermat(k) gibt die k-te Fermat-Zahl zurück. (Hinweis: k sollte nicht größer als 30 gewählt werden) */
mpz_class Fermat(unsigned int k){
    mpz_class f("2");
    mpz_t tmp;
    mpz_init(tmp);
    int exp = aux_fast_exp(2,k);
    mpz_pow_ui(tmp, f.get_mpz_t(),exp);
    mpz_class fneu (tmp);
    return fneu+1; 
}
/*Cunningham(k) gibt die k-stellige Zahl 111...1111 (vgl. Cunningham Project)*/
mpz_class Cunningham(int k){
    
    if(k > 0 and k < 10000)
    {
        std::string str("1");
        for(int j = 1; j < k; j++)
        {
            str.push_back('1');
        }
        mpz_class f(str);
        return f;
    }
    
    return mpz_class("1");
}
/* Einfache Testzahlen B1,...,B10 (siehe Moodle Lernumgebung)*/
mpz_class TestzahlB(int k){
    unsigned int N =6533;
    switch(k) {
      case 2:
            N = 10877;
            break;
      case 3:
            N = 160247;
            break;
      case 4:
            N = 438643;
            break;
      case 5:
            N = 207761;
            break;
        case 6:
            N= 1923903;
            break;
        case 7:
            N=83818463;
            break;
        case 8:
            N=117862669;
            break;
        case 9:
            N=3910882427;
            break;
        case 10:
            N =4054595993;
            break;
    }
    
    return mpz_class(N);
}
/* RSAZahl(k) gibt die RSA-Zahl R_k zurück (siehe Moodle Lernumgebung) */
mpz_class RSAZahl(int k){
    mpz_class N("1");
    switch(k) {
        case 9:
            N = "201160159";
            break;
        case 10:
            N = "3372654409";
            break;
        case 13:
            N = "1971986822969";
            break;
        case 15:
            N = "273776584181489";
            break;
        case 18:
            N = "61510304714090633";
            break;
        case 20:
            N= "4271898648953184077";
            break;
        case 25:
            N= "495559974558953081396273";
            break;
        case 30:
            N="46889352875442505708330935641";
            break;
        case 33:
            N= "74755771941480867141962416808413";
            break;
        case 36:
            N = "70229011028020481271740466987968411";
            break;
        case 40:
            N = "961278864621087840707303638982309488249";
            break;
        case 45:
            N = "55401694465066606485655927468341902477357299";
            break;
        case 50:
            N = "5255000997609473182445518818702325585174177284281";
            break;
        case 55:
            N = "270065411074759065232851708211757533053951921294341483";
            break;
        case 60:
            N = "62455028536135086396841663348071508038778737498841936587127";
            break;
        case 65:
            N = "4403549790854021189162334191093945191136870917974165910220406823";
            break;
        case 70:
            N = "887935823300379785835159151079615938803566496516605572261784750585729";
            break;
        case 100:
            N = "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";
            break;
    }
    
    return N;
}
