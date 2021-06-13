/****************************
Julia Kobusinska 165159
Agnieszka Mazur 165143
****************************/
#include <iostream>
#include <vector>
#include <bitset>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <fstream>
#define PI 3.141592654

using namespace std;

class Stan
{
private:
    bitset<4> wartosc_stanu;

public:
    Stan() {}
    Stan(bitset<4> wartosc_stanu, int hamming)
    {
        this->wartosc_stanu = wartosc_stanu;
        this->hamming = hamming;
    }

    vector<int> sciezka;
    int hamming;
    bool odwiedzony = false;

    string pobierzStan()
    {
        return this->wartosc_stanu.to_string();
    }

    int kolejnyStan(bool wartosc_pobudzenia)
    {
        return (int)(this->wartosc_stanu >> 1).set(3, wartosc_pobudzenia).to_ulong();
    }

    bitset<6> obliczBityWyjsciowe(bool wartosc_pobudzenia)
    {
        bitset<6> w;
        w[5] = (wartosc_pobudzenia + this->wartosc_stanu[3] + this->wartosc_stanu[1] + this->wartosc_stanu[0]) % 2;
        w[4] = (wartosc_pobudzenia + this->wartosc_stanu[2] + this->wartosc_stanu[0]) % 2;
        w[3] = (wartosc_pobudzenia + this->wartosc_stanu[3] + this->wartosc_stanu[2] + this->wartosc_stanu[1] + this->wartosc_stanu[0]) % 2;
        w[2] = w[5];
        w[1] = w[4];
        w[0] = w[3];

        return w;
    }
};

int liczba_bitow_w_ramce, liczba_ramek;
float max_Eb_N0, min_Eb_N0, krok_Eb_N0;
int stan_poczatkowy[4] = {0, 0, 0, 0};
int stan_koncowy[4] = {0, 0, 0, 0};
int bit_wej;
int bity_wyj[6];
int iteracja = 0;
double BER;
double BER_referencyjny;

float gauss(float mean, float sigma);
void kanal(float es_n0, long dl_kan, int *wej, float *wyj);
int obliczOdlHamminga(bitset<6> a, bitset<6> b);

ofstream fileout;

int main()
{
    fileout.open("zapis.txt");
    fileout.close();
    cout << "Podaj liczbe bitow w ramce danych trafiajacych do kodera (100/500/1000): ";
    cin >> liczba_bitow_w_ramce;
    cout << "Podaj liczbe transmitowanych ramek: ";
    cin >> liczba_ramek;
    cout << "Podaj MINIMALNA wartosc 'Eb/N0'[dB]: ";
    cin >> min_Eb_N0;
    cout << "Podaj MAKSYMALNA wartosc 'Eb/N0'[dB]: ";
    cin >> max_Eb_N0;
    cout << "Podaj krok zmiany 'Eb/N0': ";
    cin >> krok_Eb_N0;
    int ciag_informacyjny[liczba_bitow_w_ramce + 4];
    int ciag_kodowy[(liczba_bitow_w_ramce + 4) * 6];
    int ciag_odebrany[(liczba_bitow_w_ramce + 4) * 6];
    for (int i = 0; i < liczba_bitow_w_ramce * 6; i++)
    {
        ciag_kodowy[i] = 0;
    }
    srand(time(NULL));
    float kanal_wyj[(liczba_bitow_w_ramce + 4) * 6];
    for (int krok = min_Eb_N0; krok <= max_Eb_N0; krok = krok + krok_Eb_N0)
    {

        for (int k = 0; k < liczba_ramek; k++)
        {
            //informacja o postepach
            cout << "Pozostalo "<<liczba_ramek - k<<" ramek dla Eb/No = "<<krok<<endl;
            for (int i = 0; i < liczba_bitow_w_ramce; i++)
            {
                ciag_informacyjny[i] = (rand() % 2) + 0;
            }

            for (int i = liczba_bitow_w_ramce; i < liczba_bitow_w_ramce + 4; i++)
            {
                ciag_informacyjny[i] = 0;
            }


            for (int i = 0; i < liczba_bitow_w_ramce + 4; i++)
            {
                bit_wej = ciag_informacyjny[i];
                bity_wyj[0] = (bit_wej + stan_poczatkowy[0] + stan_poczatkowy[2] + stan_poczatkowy[3]) % 2;
                bity_wyj[1] = (bit_wej + stan_poczatkowy[1] + stan_poczatkowy[3]) % 2;
                bity_wyj[2] = (bit_wej + stan_poczatkowy[0] + stan_poczatkowy[1] + stan_poczatkowy[2] + stan_poczatkowy[3]) % 2;
                bity_wyj[3] = bity_wyj[0];
                bity_wyj[4] = bity_wyj[1];
                bity_wyj[5] = bity_wyj[2];
                for (int j = 0; j < 6; j++)
                {
                    ciag_kodowy[iteracja * 6 + j] = bity_wyj[j];
                }
                stan_koncowy[0] = bit_wej;
                for (int j = 1; j < 4; j++)
                {
                    stan_koncowy[j] = stan_poczatkowy[j - 1];
                }
                for (int j = 0; j < 4; j++)
                {
                    stan_poczatkowy[j] = stan_koncowy[j];
                }
                iteracja++;
            }
            for (int j = 0; j < 4; j++)
            {
                stan_poczatkowy[j] = 0;
            }
            iteracja = 0;
            //wejscie_dekodera kanalu
            kanal(krok, 6 * (liczba_bitow_w_ramce + 4), ciag_kodowy, kanal_wyj);

            //detekcja
            for (int i = 0; i < 6 * (liczba_bitow_w_ramce + 4); i++)
            {
                if (kanal_wyj[i] > 0)
                {
                    ciag_odebrany[i] = 1;
                }
                else
                {
                    ciag_odebrany[i] = 0;
                }
            }
            //sprawdzenie czy bity ciagu kodowego za kanalem zgadzaja sie z tymi przed kanalem
            for (int i = 0; i < 6 * (liczba_bitow_w_ramce + 4); i++)
            {
                if (ciag_odebrany[i] != ciag_kodowy[i])
                {
                    BER_referencyjny++;
                }
            }

            //dekoder Viterbiego
            vector<int> wejscie_dekodera(ciag_odebrany, ciag_odebrany + sizeof ciag_odebrany / sizeof ciag_odebrany[0]);

            Stan stany[16];

            for (int i = 0; i < 16; i++)
            {
                stany[i] = Stan(bitset<4>(i), 0);
            }

            stany[0].odwiedzony = true;


            while (wejscie_dekodera.size())
            {
                bitset<6> dane;

                for (int i = 0; i < 6; i++)
                {
                    dane[5 - i] = *wejscie_dekodera.begin();
                    wejscie_dekodera.erase(wejscie_dekodera.begin());
                }

                Stan tablica_pomocnicza[16];
                for (int i = 0; i < 16; i++)
                {
                    tablica_pomocnicza[i] = Stan(bitset<4>(i), 0);
                }

                for (int i = 0; i < 16; i++)
                {
                    if (!stany[i].odwiedzony)
                    {
                        continue;
                    }

                    for (int wartosc_pobudzenia = 0; wartosc_pobudzenia < 2; wartosc_pobudzenia++)
                    {
                        int nowy_stan;
                        bitset<6> bity_wyjsciowe;

                        nowy_stan = stany[i].kolejnyStan(wartosc_pobudzenia);
                        bity_wyjsciowe = stany[i].obliczBityWyjsciowe(wartosc_pobudzenia);

                        int hamming = stany[i].hamming + obliczOdlHamminga(dane, bity_wyjsciowe);

                        if (!tablica_pomocnicza[nowy_stan].odwiedzony || hamming < tablica_pomocnicza[nowy_stan].hamming)
                        {
                            tablica_pomocnicza[nowy_stan].odwiedzony = true;
                            tablica_pomocnicza[nowy_stan].hamming = hamming;
                            tablica_pomocnicza[nowy_stan].sciezka = stany[i].sciezka;
                            tablica_pomocnicza[nowy_stan].sciezka.push_back(wartosc_pobudzenia);
                        }
                    }
                }

                for (int i = 0; i < 16; i++)
                {
                    if (!tablica_pomocnicza[i].odwiedzony)
                    {
                        continue;
                    }

                    stany[i].hamming = tablica_pomocnicza[i].hamming;
                    stany[i].sciezka = tablica_pomocnicza[i].sciezka;
                    stany[i].odwiedzony = tablica_pomocnicza[i].odwiedzony;
                }

                //cout<<"Pozostalo " << wejscie_dekodera.size() << " bitów do zdekodowania w obecnej ramce" << endl;
            }


            //sprawdzenie czy bity ciagu informacyjnego zgadzaja sie z bitami ciagu zdekodowanego
            for (int i = 0; i < (liczba_bitow_w_ramce + 4); i++)
            {
                if (ciag_informacyjny[i] != stany[0].sciezka[i])
                {
                    BER++;
                }
            }

        }
        cout << endl << "Eb/N0=" << krok << " BER referencyjny=" << BER_referencyjny / ((liczba_bitow_w_ramce + 4) * liczba_ramek * 6) << endl;
        cout << endl << "Eb/N0=" << krok << " BER =" << BER / ((liczba_bitow_w_ramce + 4) * liczba_ramek ) << endl;
        fileout.open("zapis.txt", ios_base::app);
        fileout << "Eb/N0=" << krok << " BER referencyjny=" << BER_referencyjny / ((liczba_bitow_w_ramce + 4) * liczba_ramek * 6) << "Eb/N0=" << krok << " BER =" << BER / ((liczba_bitow_w_ramce + 4) * liczba_ramek )<< endl;
        fileout.close();
        BER_referencyjny = 0;
        BER = 0;
    }
    system("pause");
}

//************************************************************************

// Function kanal changes binary values into bipolar ones (-1/+1) and adds noise
// *wej - Input vector of binary values (0/1)
// *wyj - Output vector of real numbers
// es_n0 - Es/N0
// dl_kan - the number of input bits
void kanal(float es_n0, long dl_kan, int *wej, float *wyj)
{
    float mean = 0;
    float es = 1;
    float sygnal;
    float sigma;
    float s_n;
    long y;

    s_n = (float)pow(10, (es_n0 / 10));
    sigma = (float)sqrt(es / (2 * s_n));

    for (y = 0; y < dl_kan; y++)
    {
        sygnal = 2 * *(wej + y) - 1;              // change the binary value (0/1) into symbol (-1/+1)
        *(wyj + y) = sygnal + gauss(mean, sigma); // noise addition
    }
}
//**********************************************************************
float gauss(float mean, float sigma)
{
    double x;
    double z;

    z = (double)rand() / RAND_MAX;
    if (z == 1.0)
        z = 0.9999999;
    x = sigma * sqrt(2.0 * log(1.0 / (1.0 - z)));

    z = (double)rand() / RAND_MAX;
    if (z == 1.0)
        z = 0.9999999;
    return ((float)(mean + x * cos(2 * PI * z)));
}
//*******************************************************************
int obliczOdlHamminga(bitset<6> a, bitset<6> b)
{
    int hamming = 0;
    for (int i = 0; i < 6; i++)
    {
        if (a[i] != b[i])
        {
            hamming++;
        }
    }

    return hamming;
}
//*******************************************************************
