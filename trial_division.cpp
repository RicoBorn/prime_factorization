//
// HIER KOMPLETT AUFRÄUMEN !!!!!!!!
//

#include "trial_division.h"


/* Datei "probedivision.cpp" enthält ein Programm zur Faktorisierung von Zahlen mittels Probedivision
Kompilieren: g++ probedivision.cpp testzahlen.cpp -lgmp -lgmpxx
Ausführen mit Parameter:
 -F = Fermatzahl F_k
 -C = Cunninghamzahl (1111...111)
 -B = Testzahl b_k
 -R = RSA-Zahl R_k
 ohne Parameter: freie Eingabe einer Zahl.

 (Beispielprogramm von Steffen Kionke für das Praktikum Algebra 2022)
 */

#include <iostream>
#include <gmpxx.h>
#include <list>
#include <chrono>
#include <iterator>
#include "testzahlen.h"

/*
 Klasse "Factor": Kann genutzt werden um einen (Prim)Faktor mit Exponenten zu speichern.
 */
class Factor {
  public:
    unsigned int exponent;  // Exponent: erlaubt es Faktoren mit Potenzen zu zählen (2^6)
    mpz_class factor; // Faktor (im besten Fall: eine Primzahl)
    int is_prime; // Ist "factor" eine Primzahl? 2 sicher, 1 wahrscheinlich, 0 nein
    Factor()
    {   exponent = 0;
        is_prime = 0;
        factor = mpz_class("1");
    }
    void printpp()
    {
        char open = '(';
        char close = ')';

        if(is_prime == 0) // Faktor ist nicht prim
        {
            open = '[';
            close = ']';
        }


        std::cout << open << factor;

         if(exponent > 1)
                  std::cout << "^" << exponent;
         if(is_prime == 1) // nicht 100% sicher prim
                  std::cout << "_?";

        std::cout << close << " ";
    }
};

/*
 divide_out_maximal_power takes a number T and a (prime) number P returns maximal exponent e such that P^e divides T. In the end T is T/P^e.
 */
unsigned int divide_out_maximal_power(mpz_class *t, mpz_class P)
{
    unsigned int exponent = 0;
    while( (*t)%P == 0)
    {
        *t = (*t)/P;
        exponent++;
    }
    return exponent;
}
/*
 Faktorisiert die Zahl N mittels Probedivision bis zur Schranke B und prüft, ob der Rest eine Primzahl ist. Gibt die Liste der gefundenen Primfaktoren zurück. (Für große Zahlen N muss eine sinnvolle Schranke B gewählt werden, sonst explodiert die Laufzeit!)
*/
std::list<Factor> trial_division_bounded(mpz_class N, mpz_class B) {
    mpz_class tmp = N;
    mpz_class C = sqrt(N)+1;
    if(B > C)
        B = C;

    std::list<Factor> factors;
    Factor ppt;

    // Prüfen: Ist N eine Primzahl (27: 1 Baille-PSW und 2 Rabin-Müller Tests)
    int a = mpz_probab_prime_p(N.get_mpz_t(), 27);
    if( a > 0) // N ist wahrscheinlich prim
    {
        ppt.exponent = 1;
        ppt.factor = N;
        ppt.is_prime = a;
        factors.push_back(ppt);
        return factors;
    }

    // Andernfalls: Faktorisiere N mit Probedivision
    // Erst 2 und 3, dann alle Zahlen kongruent 1 oder 5 modulo 6
    mpz_class P("2");
    ppt.exponent = divide_out_maximal_power(&tmp,P);
    if(ppt.exponent > 0)
    {
        ppt.factor = P;
        ppt.is_prime = 2;
        factors.push_back(ppt);

    }
    P++; // P = 3
    ppt.exponent = divide_out_maximal_power(&tmp,P);
    if(ppt.exponent > 0)
    {
        ppt.factor = P;
        ppt.is_prime = 2;
        factors.push_back(ppt);

    }

    P+=3; // P = 6
    while(tmp > 1 and P <= B )
    {
        ppt.exponent = divide_out_maximal_power(&tmp,P-1);
        if(ppt.exponent > 0)
        {
            ppt.factor = P-1;
            ppt.is_prime = 2;
            factors.push_back(ppt);

        }
        ppt.exponent = divide_out_maximal_power(&tmp,P+1);
        if(ppt.exponent > 0)
        {
            ppt.factor = P+1;
            ppt.is_prime = 2;
            factors.push_back(ppt);
        }
        P += 6;
    }

    // Prüfen, ob der Rest eine Primzahl ist.
    a = mpz_probab_prime_p(tmp.get_mpz_t(), 27);
    if( tmp > 1) // tmp is probably prime
    {
        ppt.exponent = 1;
        ppt.factor = tmp;
        ppt.is_prime=a;
        factors.push_back(ppt);
        return factors;
    }

    return factors;
}
/*
 Faktorisiert die Zahl N vollständig mittels Probedivision und gibt die Liste der Primfaktoren zurück. (Nur für kleine Zahlen (etwa 10^15) funktioniert das sicher in annehmbarer Zeit)
 */
std::list<Factor> trial_division(mpz_class N) {
    return trial_division_bounded(N,N);
}
/*
 Main: Parameter
    -F = Fermatzahlen
    -C Cunninghamzahl
    -B Testzahl b_k
    -R RSA-Zahl R_k
 */
int main(int argc, char *argv[]) {

    int mode = 0; // 0 Eingabe, 1 Fermat, 2 Cunningham, 3 TestzahlB, 4 RSA-Zahl
    std::string frage("Zahl k");
    if(argc > 1)
    {
        std::string arg1(argv[1]);
        if(arg1 == "-F")
        {
            mode = 1;
            frage = "Fermat-Zahl F_k";

        }
        else if(arg1 == "-C")
        {
            mode = 2;
            frage = "Cunningham-Zahl C_k";
        }
        else if(arg1 == "-B")
        {
            mode = 3;
            frage = "Testzahl B_k";
        }
        else if(arg1 == "-R")
        {
            mode = 4;
            frage = "RSA-Zahl R_k";
        }
    }

    std::string str;

    std::cout << "Welche " << frage << " soll faktorisiert werden? \n";

    std::cout  << "Eingabe k:  " ;
    std::cin >> str;
    mpz_class N("0");

    // Hier werden die Testzahlen generiert.
    // Die Funktionen dazu sind in testzahlen.h und testzahlen.cpp zu finden.
    if(mode == 1)
    {
        N = Fermat(std::stoi(str));
    }
    else if(mode == 2)
    {
        N = Cunningham(std::stoi(str));
    }
    else if(mode == 3)
    {
        N = TestzahlB(std::stoi(str));
    }
    else if(mode == 4)
    {
        N = RSAZahl(std::stoi(str));
    }
    else
    {
        N.set_str(str,0);
    }

    // Schranke: Probedivision bis 10^8
    // Kann erhöht werden, die Laufzeit wird aber schnell zu groß
    mpz_class B("100000000");

    // Zeitmessung starten---------
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();

    //Probedivision durchführen
    std::list<Factor> factors = trial_division_bounded(N,B);

    //Zeitmessung stoppen----------
    std::chrono::time_point<std::chrono::high_resolution_clock>
 end = std::chrono::high_resolution_clock::now();
    //Vergangene Zeitspanne bestimmen und ausgeben
    std::chrono::duration<double, std::milli> float_ms = end - start;
    std::cout << "Berechnung in " << float_ms.count() << " Millisekunden abgeschlossen" << std::endl;

    std::cout << "Ermittelte Zerlegung von " << N << ":\n";

    std::list<Factor>::iterator it;
    for (it = factors.begin(); it != factors.end(); it++)
        it->printpp();


  return EXIT_SUCCESS;
}

