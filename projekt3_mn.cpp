#include <iostream>
#include <math.h>
#include<string>
#include<fstream>
#include<vector>

//#define LICZBA_WEZLOW 5	
#define MAX 512				
//#define LICZBA_PRZEDZIALOW (LICZBA_WEZLOW-1)
#define KROK 1
#define NAZWA "trasa_3"

using namespace std;

void zamiana()
{
	string linia;
	string nazwa = NAZWA;
	ifstream plik(nazwa + ".txt");
	ofstream nowy(nazwa + "_wynik.txt");

	if (plik.good())
	{
		while (!plik.eof())
		{
			getline(plik, linia);
			for (int i = 0; i < linia.size(); i++)
			{
				if (linia[i] == ',')
				{
					linia[i] = ' ';
				}
			}

				nowy << linia<<endl;
		}
	
	
	}
	else
		cout << "Nie mozna wczytaj takiego pliku.\n";
	cout << "skonczylem wczytywac\n";

	plik.close();
	nowy.close();
}

void wczytaj_plik(double *x, double *y, double *punkty, int liczba_przedzialow)
{
	string nazwa = NAZWA;
	ifstream plik(nazwa+".txt");


	if (plik.good())
	{
		double tmp_y[MAX];
		double dystans, wysokosc;
		int i = -1;
		int licznik = 0;
		while (plik >> dystans >> wysokosc)
		{
			i++;
			punkty[i] = dystans;
			tmp_y[i] = wysokosc;
		}


		double skok = MAX / liczba_przedzialow;
		int k = 0;
		for (double j = 0; j < skok*liczba_przedzialow; k++, j += skok)
		{
			x[k] = punkty[(int)j];
			y[k] = tmp_y[(int)j];
		}

		x[k] = punkty[MAX - 1];
		y[k] = tmp_y[MAX - 1];
		
		plik.close();
	}
	else
		cout << "Nie mozna wczytaj takiego pliku.\n";
	cout << "skonczylem wczytywac\n";
}

void zapisz_l(double *x, double *y)
{
	string nazwa_pliku=NAZWA;
	ofstream plik("wyniki_1/rezultat_lg_1.txt");
	for (int i = 0; i < MAX; i++)
	{
		plik << x[i] << " " << y[i] << endl;
	}
	plik.close();
}

void zapisz_s(vector<double> x, vector<double> y)
{
	string nazwa_pliku = NAZWA;
	ofstream plik("wyniki_1/rezultat_wiel_1.txt");
	for (int i = 0; i < x.size(); i++)
	{
		plik << x.at(i) << " " << y.at(i) << endl;
	}
	plik.close();
}

void LU(int size, double **A, double *b, double *v_x)
{
	double **L, **U, **P;
	L = new double*[size];
	U = new double*[size];
	P = new double*[size];

	for (int i = 0; i < size; i++)
	{
		L[i] = new double[size];
		U[i] = new double[size];
		P[i] = new double[size];

		for (int j = 0; j < size; j++)
		{
			if (i == j)
			{
				L[i][j] = 1;
				P[i][j] = 1;
			}
			else
			{
				L[i][j] = 0;
				P[i][j] = 0;
			}

			U[i][j] = A[i][j];
		}
	}

	for (int i = 0; i < size - 1; i++)
	{
		double pivot = 0;
		int p_ind;
		for (int j = i; j < size; j++)
		{
			if (abs(U[j][i]) > pivot)
			{
				pivot = abs(U[j][i]);
				p_ind = j;
			}
		}

		for (int j = 0; j < size; j++)
		{
			if (j >= i)
				swap(U[i][j], U[p_ind][j]);
			else
				swap(L[i][j], L[p_ind][j]);

			swap(P[i][j], P[p_ind][j]);
		}

		for (int j = i + 1; j < size; j++)
		{
			L[j][i] = U[j][i] / U[i][i];
			for (int k = i; k < size; k++)
			{
				U[j][k] -= L[j][i] * U[i][k];
			}
		}
	}

	double *b_temp;
	b_temp = new double[size];
	for (int i = 0; i < size; i++)
	{
		b_temp[i] = 0;
		for (int j = 0; j < size; j++)
			b_temp[i] += b[j] * P[i][j];
	}
	for (int i = 0; i < size; i++)
		b[i] = b_temp[i];

	delete[] b_temp;

	// podstawienie wprzod
	double *d = new double[size];
	for (int i = 0; i < size; i++)
	{
		d[i] = b[i];
		for (int j = 0; j < i; j++)
			d[i] -= L[i][j] * d[j];

		d[i] /= L[i][i];
	}

	//podstawienie wstecz
	for (int i = size - 1; i >= 0; i--)
	{
		v_x[i] = d[i];
		for (int j = i + 1; j < size; j++)
		{
			v_x[i] -= U[i][j] * v_x[j];
		}
		v_x[i] /= U[i][i];
	}

	delete[]d;
	for (int i = 0; i < size; i++)
	{
		delete[]L[i];
		delete[]U[i];
		delete[]P[i];
	}
	delete[]L;
	delete[]U;
	delete[]P;
}

double fi_i(double *dystans, int iteracja, double punkt, int liczba_wezlow)
{
	double wynik = 1;
	double mianownik = 1;
	double licznik = 1;

	for (int i = 0; i < liczba_wezlow; i++)
	{
		if(dystans[iteracja] != dystans[i])
			mianownik *= (dystans[iteracja] - dystans[i]);
		if (iteracja == i)
			continue;
		licznik *= (punkt - dystans[i]);
	}

	wynik = licznik / mianownik;

	return wynik;
}

double interpolcjaLagrange(double *dystans, double *wysokosc, double punkt, int liczba_wezlow)
{
	double *wynik=new double[liczba_wezlow];
	double wielomian;
	double dst_y = 0;


	for (int i = 0; i < liczba_wezlow; i++)
	{
		wielomian = fi_i(dystans, i, punkt,liczba_wezlow);
		dst_y +=  wielomian*wysokosc[i];
	}
	delete[]wynik;
	return dst_y;
}

void wypelnij_macierz(double **macierzA, double *wektorB, double *wysokosc, double *dystans, int liczba_przedzialow)
{
	for (int i = 0; i < liczba_przedzialow * 4; i++)
	{
		wektorB[i] = 0;
		for (int j = 0; j < liczba_przedzialow * 4; j++)
		{
			macierzA[i][j] = 0;
		}
	}

	for (int j = 0; j < liczba_przedzialow; j++)
	{
		double dx = dystans[j + 1] - dystans[j];

		//rownanie Sj(xj)=f(xj)				(1)
		macierzA[j][4 * j] = 1;
		wektorB[j] = wysokosc[j];

		//rownanie Sj(xj+1) = f(xj+1)		(2)
		macierzA[liczba_przedzialow + j][4 * j] = 1;
		macierzA[liczba_przedzialow + j][4 * j + 1] = dx;
		macierzA[liczba_przedzialow + j][4 * j + 2] = dx*dx;
		macierzA[liczba_przedzialow + j][4 * j + 3] = dx*dx*dx;
		wektorB[liczba_przedzialow + j] = wysokosc[j + 1];

		if (j < liczba_przedzialow - 1)
		{
			// rownanie (3)
			macierzA[2 * liczba_przedzialow + j][4 * j + 1] = 1;
			macierzA[2 * liczba_przedzialow + j][4 * j + 2] = 2 * dx;
			macierzA[2 * liczba_przedzialow + j][4 * j + 3] = 3 * dx*dx;
			macierzA[2 * liczba_przedzialow + j][4 * j + 5] = -1;
			wektorB[2 * liczba_przedzialow + j] = 0;

			//rowanie (4)
			macierzA[3 * liczba_przedzialow - 1 + j][4 * j + 2] = 1;
			macierzA[3 * liczba_przedzialow - 1 + j][4 * j + 3] = 3 * dx;
			macierzA[3 * liczba_przedzialow - 1 + j][4 * j + 6] = -1;
			wektorB[3 * liczba_przedzialow - 1 + j] = 0;
		}
		if (j == 0)
		{	//rownanie (5)
			macierzA[4 * liczba_przedzialow - 2][2] = 1;
			wektorB[4 * liczba_przedzialow - 2] = 0;
		}
		else if (j == liczba_przedzialow - 1)
		{	//rownanie (6)
			macierzA[4 * liczba_przedzialow - 1][4 * liczba_przedzialow - 2] = 1;
			macierzA[4 * liczba_przedzialow - 1][4 * liczba_przedzialow - 1] = 3 * dx;
			wektorB[4 * liczba_przedzialow - 1] = 0;
		}
	}

}

void interpolacjaSplajnami(double *dystans, double *wysokosc, vector<double> punkty_s, vector<double> wyniki_s, int liczba_przedzialow)
{
	double** macierzA = new double*[liczba_przedzialow * 4];
	double *wektorX = new double[4 * liczba_przedzialow];
	double *wektorB = new double[4 * liczba_przedzialow];

	for (int i = 0; i < liczba_przedzialow * 4; i++)
	{
		macierzA[i] = new double[liczba_przedzialow * 4];
		wektorX[i] = 0;
	}

	wypelnij_macierz(macierzA, wektorB, wysokosc, dystans, liczba_przedzialow);

	LU(4 * liczba_przedzialow, macierzA, wektorB, wektorX);

	for (int i = 0; i < liczba_przedzialow; i++)
	{
		for (double k = dystans[i]; k < dystans[i + 1]; k += KROK)
		{
			double delta = k - dystans[i];
			double wynik_delty = wektorX[4*i] + wektorX[4*i+1]*delta+wektorX[4*i+2]*delta*delta+wektorX[4*i+3]*delta*delta*delta;
			punkty_s.push_back(delta);
			wyniki_s.push_back(wynik_delty);
		}
	}

	zapisz_s(punkty_s, wyniki_s);
	delete[]wektorB;
	delete[]wektorX;
	for (int j = 0; j < 4 * liczba_przedzialow; j++)
		delete[] macierzA[j];
	delete[]macierzA;

}


int main(int argc, char **argv)
{
	int liczba_wezlow=10;
	if (argc > 1)
		liczba_wezlow = atoi(argv[1]);

	int liczba_przedzialow = liczba_wezlow - 1;

	double *dystans = new double[liczba_wezlow];
	double *wysokosc = new double[liczba_wezlow];
	double punkty[MAX];
	double wyniki[MAX];
	vector<double> wyniki_s;
	vector<double> punkty_s;
	int ostatni;
	
	wczytaj_plik(dystans, wysokosc, punkty,liczba_przedzialow);
	
	for (int i = 0; i < liczba_wezlow; i++)
	{
		cout << dystans[i] << " " << wysokosc[i] << endl;
	}	

	cout.precision(20);

	//for (int i = 0; i < MAX; i++)
		//wyniki[i] = interpolcjaLagrange(dystans, wysokosc, punkty[i], liczba_wezlow);
	//zapisz_l(punkty, wyniki);

	interpolacjaSplajnami(dystans, wysokosc, punkty_s, wyniki_s,liczba_przedzialow);
	delete[]dystans;
	delete[]wysokosc;

	return 0;
}