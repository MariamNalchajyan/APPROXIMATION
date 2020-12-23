#include <iostream>
#include <math.h>
#include <vector>
#include <iomanip>
#define n 4
#define m 10
#define a 0.0
#define b 2.0
#define c 5

using namespace std;

void swap(double &s, double &t)
{
	double tmp = s;
	s = t;
	t = tmp;
}

double func(double x)
{
	return exp(-x) * cos(x);
}

vector<double> xi()
{
	double h = (b - a) / m;
	vector<double> vect_x;
	for (int i = 0; i <= m; i++)
	{
		double x = a + i * h;
		vect_x.push_back(x);
	}
	return vect_x;
}
//-----------------------------------------------------------------------------------
vector<double> yi(vector<double> vect_x)
{
	double z = 0;
	vector<double> vect_y;
	double y;
	for (int i = 0; i <= m; i++)
	{
		if (z < 0.5)
		{
			y = func(vect_x.at(i)) * (1 - z / c);
		}
		else
		{
			y = func(vect_x.at(i)) * (1 + z / c);
		}
		vect_y.push_back(y);
		z += 1.0 / m;
	}
	return vect_y;
}
//--------------------------------------------------------------
vector<double> d(double x)
{
	vector<double> v;
	v.push_back(1);
	double mult = 1; 
		for (int k = 1; k <= 2*n; k++)
		{
			mult =mult * x;
			v.push_back(mult);
		}
		return v;
}


vector<double> Dk(vector<double> x)
{
	vector<double> result;
	for (int i = 0; i <= 2 * n; i++)
	{
		result.push_back(0);
	}
	
	for (int k = 1; k <= m; k++)
	{
		vector<double> v = d(x[k]);
		for (int i = 0; i <= 2 * n; i++)
		{
			result[i] += v[i];
		}
	}
	return result;
}
//-----------------------------------------------------
vector<double> B(double x,double y)
{
	vector<double> v = d(x);
	for (int i = 0; i <=n; i++)
	{
		v[i] *= y;
	}
	return v;
}

vector<double> Bk(vector<double> x, vector<double> y)
{
	vector<double> result;
	for (int i = 0; i <= n; i++)
	{
		result.push_back(0);
	}
	for (int k = 1; k <= m; k++)
	{
		vector<double> v = B(x[k],y[k]);
		for (int i = 0; i <= n; i++)
		{
			result[i] += v[i];
		}
	}
	return result;
}
//----------------------------------------------------------------
vector<vector<double> > matrix(vector<double> Dk, vector<double> Bk)
{
	vector<vector<double> > A;
	for (int i = 0; i <= n; i++)
	{
		vector<double> temp;
		for (int j = 0; j <= n + 1; j++)
		{
			if (j < n + 1)
			{
				temp.push_back(Dk[i + j]);
			}
			else
			{
				temp.push_back(Bk[i]);
			}
		}
		A.push_back(temp);
	}
	return A;
}
//-------------Перестановка строк в матрице---------------------------------
void max0(vector<vector<double> > &matr)
{
	int p[n + 1];
	double max;
	for (int i = 0; i <= n; i++)
	{
		max = 0;
		for (int j = i; j <= n; j++)
		{
			if (fabs(matr[j][i]) > max)
			{
				max = fabs(matr[j][i]);
				p[i] = j;
			}
		}
	}
	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n + 1; j++)
		{
			swap(matr[i][j], matr[p[i]][j]);
		}
	}
}
//-----------------метод Гаусса----------
vector<double> C(vector<vector<double> > matr)
{
	double x[n + 1];
	double kf;
	for (int j = 0; j <= n; j++)
	{
		for (int i = j + 1; i <= n; i++)
		{
			kf = matr[i][j] / matr[j][j];
			for (int k = j; k <= n + 1; k++)
			{
				matr[i][k] -= kf * matr[j][k];
			}
		}
	}
	double sum;
	for (int i = n; i >= 0; i--)
	{
		sum = 0;
		for (int j = i + 1; j <= n; j++)
		{
			sum += x[j] * matr[i][j];
		}
		x[i] = (matr[i][n + 1] - sum) / matr[i][i];
	}
	vector<double> tmp;
	for (int i = 0; i <= n; i++)
	{
		tmp.push_back(x[i]);
	}
	return tmp;
}
//----------Невязки------------
vector<double> AC(vector<vector<double> > matr, vector<double> vect_c)
{
	double sum;
	vector<double> res;
	for (int i = 0; i <= n; i++)
	{
		sum = 0;
		for (int j = 0; j <= n; j++)
		{
			sum += matr[i][j] * vect_c[j];
		}	
		res.push_back(sum);
	}
	return res;
}
vector<double> F(vector<double> vect_ac, vector<double> vect_b)
{
	vector<double> res;
	double f[n + 1];
	for (int i = 0; i <= n; i++)
	{
		f[i] = vect_ac[i] - vect_b[i];
		res.push_back(f[i]);
	}
	return res;
}
//-----------Y--------
double Y1(double x,vector<double> vect_c)
{
	vector<double> v;
	v.push_back(1);
	double mult = 1;
	for (int k = 1; k <= 2 * n; k++)
	{
		mult = mult * x;
		v.push_back(mult);
	}
	vector<double> y;
	for (int i = 0; i <= n; i++)
	{
		mult = v[i] * vect_c[i];
		y.push_back(mult);
	}
	double sum = 0;
	for (int k = 0; k <= n; k++)
	{
		sum += y[k];
	}
	return sum;
}

vector<double> Y(vector<double> vect_x, vector<double> vect_c)
{
	vector<double> result;
	for (int k = 0; k <= m; k++)
	{
		double v = Y1(vect_x[k],vect_c);
		result.push_back(v);
	}
	return result;
}

//--------S--------------------
double S(vector<double> Y, vector<double> fi)
{
	double sum = 0;
	for (int i = 0; i <= m; i++)
	{
		sum += ((Y[i] - fi[i]) * (Y[i] - fi[i]));
	}
	return sum;
}



int main()
{
	vector<double> x = xi();
	vector<double> f = yi(x);	
	cout << "--------------------" << endl;
	cout << left << setw(5) << "x(i)"<< " | " << left << setw(5) << "f(i)" << endl<<endl;
	for (int i = 0; i <= m; i++)
	{
		cout << left << setw(5) << x.at(i);
		cout << " | " << left << setw(5) << f.at(i) << endl;
	}

	cout<< endl << "--------------Vector Dk------------ " << endl;
	vector<double> D = Dk(x);
	for (int i = 0; i <= 2 * n; i++)
	{
		cout << D.at(i) << endl;;
	}

	vector<double> P = Bk(x, f);
	cout<<endl<< "-------------Vector Bk-------------- "<< endl;
	for (int i = 0; i <= n; i++)
	{
		cout<< P.at(i) << endl;
	}

	vector<vector<double> > A = matrix(D, P);
	cout << endl<< "------Augmented Matrix------" << endl;
	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n + 1; j++)
		{
			cout << left << setw(8) << A[i][j] << "  ";
		}
		cout << endl;
	}
	//max0(A);
	/*
	cout << endl << "------New Augmented Matrix------" << endl;
	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n + 1; j++)
		{
			cout << left << setw(8) << A[i][j] << "  ";
		}
		cout << endl;
	}
	*/
	
	vector<double> Xk = C(A);
	cout << endl << "------Vector Ck--------" << endl;
	
	for (int j = 0; j <= n; j++)
	{
		cout << Xk[j] << endl;
	}
	
	vector<double> mult_ac = AC(A,Xk);
	cout << endl << "------Vector A*C--------" << endl;
	for (int j = 0; j <= n; j++)
	{
		cout << mult_ac[j] << endl;
	}
	vector<double> nev = F(mult_ac, P);
	cout << endl << "------Nevyazki---A*C-B-----" << endl;
	for (int j = 0; j <= n; j++)
	{
		cout << nev[j] << endl;
	}
	cout<< endl<< left << setw(5) << "x" << " | " << left << setw(5) << "Y" << endl;
	vector<double> apr = Y(x,Xk);
	for (int i = 0; i <= m; i++)
	{
		cout << left << setw(5) << x[i];
		cout << " | " << left << setw(5) << apr[i] << endl;
	}

	double s = S(apr, f);
	cout << "s = " << s << endl;
	return 0;
}
