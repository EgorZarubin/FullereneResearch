#include <iostream>

#include <math.h>
#include <iomanip>
#include <fstream>
#include <conio.h>
#include <string.h>
#include <string>
#include <vector>

using namespace std;

double const pi = 3.1415926;
double const l = 1.4;
int Amount = 0;
double coeff = 1.14;
double coeff2 = 0;
char str[2] = "C";
double rv;

const bool ROUND_SNAIL = true;
class PolarPoint;

class Point3 // точка в декартовых координатах
{
public:
	double x, y, z;
	Point3() {};
	Point3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};
	~Point3() {};
	double length() const
	{
		return sqrt(x * x + y * y + z * z);
	}
	//PolarPoint TransToPolar (); //функция преобразования в полярные координаты
	//double distance(Point3 A, Point3 B); // функция опредяляющая расстояние между точками в декартовых координатах		
};

typedef Point3 Vector3 ;
Vector3 operator-(const Vector3& p1, const Vector3& p2)
{
	return Vector3(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}
Vector3 operator+(const Vector3& p1, const Vector3& p2)
{
	return Vector3(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}
Point3 operator/(const Vector3& p, double n) {
	return Point3(p.x / n, p.y / n, p.z / n);
}
Point3 operator*(const Vector3& p, double n) {
	return Point3(p.x * n, p.y * n, p.z * n);
}
Point3 operator*(double n, const Vector3& p) {
	return Point3(p.x * n, p.y * n, p.z * n);
}

Point3 round(const Point3& p, int K)
{
	double r_factor = coeff * (asin(p.z / sqrt(p.x * p.x + p.z * p.z)) + pi * K) / (pi * p.length());
	return r_factor * p;
}


class PolarPoint // точка в параметризованных координатах
{
public:
	double r, fi, psi;

	Point3 transToDecart(PolarPoint B); //функция преобразования в декартовы координаты
	double polarDistance(PolarPoint A, PolarPoint B); // функция опредяляющая расстояние между точками в параметрических координатах	
};

Point3 PolarPoint::transToDecart(PolarPoint B)
{
	Point3 A;
	A.x = B.r * sin(B.psi) * cos(B.fi);
	A.y = B.r * cos(B.psi);
	A.z = B.r * sin(B.psi) * sin(B.fi);
	return A;
}


class SnailStructure
{
	FILE* file;
	Point3 A, B, C; // вершины треугольника
	int N; // к-во шестиугольников в наименьшей хорошей стороне
	
	vector<Point3> basePoints;
	int num;

public:
	SnailStructure() {
		fopen_s(&file, "result1.out", "w");
		fprintf(file, "%s %s %s  %s  %s  %s %s \n", "ATOM", "TYPE", "POSITION(Astrm)", "VELOCITY(km/s)", "KE(eV)", "PE(eV)", "TE(eV)");
	}

	SnailStructure(const char * fileName) {
		fopen_s(&file, fileName, "w");
		fprintf(file, "%s %s %s  %s  %s  %s %s \n", "ATOM", "TYPE", "POSITION(Astrm)", "VELOCITY(km/s)", "KE(eV)", "PE(eV)", "TE(eV)");
	}

	~SnailStructure() {
		fprintf(file, "\n \n \n \n");
		fclose(file);
	}

	void writePoint(Point3 p, int K) const
	{
		if (ROUND_SNAIL && K != 0) p = round(p, K);

		Amount++;
		
		fprintf(file, "%d %s %15.4f  %15.4f  %15.4f  %d %d %d %d %d %d %d %d \n", Amount, str, p.x, p.y, p.z, 0, 0, 0, 0, 0, 0, 1, 0);
		//fprintf(file, "%15.4f  %15.4f  %15.4f \n", p.x, p.y, p.z);
	}

	vector<Point3> defineBasePoints_(int k)
	{
		PolarPoint PIcoV[3];            //icosaeder vertex in polar coordinates
   //double x,y,z;
		double tan_value = 1.618;

		if (k % 2 == 0)
		{
			PIcoV[0].fi = pi * k / 2;
			PIcoV[0].psi = atan(tan_value);
			PIcoV[0].r = k * coeff + coeff2;

			PIcoV[1].fi = pi * k / 2;
			PIcoV[1].psi = pi - atan(tan_value);
			PIcoV[1].r = k * coeff + coeff2;

			PIcoV[2].fi = pi * k / 2 + atan(tan_value);
			PIcoV[2].psi = pi / 2;
			PIcoV[2].r = k * coeff + coeff2;
		}
		else if (k != 1)
		{
			PIcoV[0].fi = pi * k / 2;
			PIcoV[0].psi = atan(1 / tan_value);
			PIcoV[0].r = (k - 1) * coeff + coeff2;

			PIcoV[1].fi = pi * k / 2;
			PIcoV[1].psi = pi - atan(1 / tan_value);
			PIcoV[1].r = (k - 1) * coeff + coeff2;

			PIcoV[2].fi = pi * k / 2 + atan(1 / tan_value);
			PIcoV[2].psi = pi / 2;
			PIcoV[2].r = (k + 1) * coeff + coeff2;
		}
		else if (k == 1)
		{
			PIcoV[0].fi = pi * k / 2;
			PIcoV[0].psi = atan(1 / tan_value);
			PIcoV[0].r = 2 * coeff + coeff2;

			PIcoV[1].fi = pi * k / 2;
			PIcoV[1].psi = pi - atan(1 / tan_value);
			PIcoV[1].r = 2 * coeff + coeff2;

			PIcoV[2].fi = pi * k / 2 + atan(1 / tan_value);
			PIcoV[2].psi = pi / 2;
			PIcoV[2].r = 2 * coeff + coeff2;
		}

		vector<Point3> IcoV;
		for (int i = 0; i < 3; i++)
		{
			double x = PIcoV[i].r * sin(PIcoV[i].psi) * cos(PIcoV[i].fi);
			double y = PIcoV[i].r * cos(PIcoV[i].psi);
			double z = PIcoV[i].r * sin(PIcoV[i].psi) * sin(PIcoV[i].fi);
			IcoV.push_back(Point3(x, y, z));
			//IcoV[i]=TransToDecart(PIcoV[i]);

		}
		return IcoV;
	}

	void defineBasePoints(int k)
	{
		num = k;
		basePoints.resize(3 * (num + 5));
		for (int i = 1; i < num + 6; i++)                 // вычисление координат вершин, и записи их в массив начиная с угла 2пи
		{
			vector<Point3> IcoV = defineBasePoints_(i);
			for (int j = 0; j < 3; j++)
			{
				basePoints[3 * (i - 1) + j] = IcoV[j];
			}
		}
	}

	

	void triangle1(Point3 A, Point3 B, Point3 C, int N, int K = 0);
	void triangle2(Point3 A, Point3 B, Point3 C, int N, int K = 0);
	void triangle3(Point3 A, Point3 B, Point3 C, int N, int K = 0);
	void triangle2_1(Point3 A, Point3 B, Point3 C, int N, int K = 0);
	void triangle2_2(Point3 A, Point3 B, Point3 C, int N, int K = 0);

	void build();

};

void SnailStructure::triangle1(Point3 A, Point3 B, Point3 C, int N, int K)   //функция замощения молекулами правильного треугольника со стороной состоящей из N шестиугольников
{
	Point3 I = (B - A) / (N + 1);
	Point3 J = (C - A) / (N + 1);              //basis vectors in triangle plane
	
	for (int i = 0; i < N + 1; i++)
	{
		Point3 p = A + (I + J) / 3 + i * I;

		if (ROUND_SNAIL && K != 0)
		{
			p = round(p, K);
		}
		
		Amount++;
		fprintf(file, "%d %s %15.4f  %15.4f  %15.4f  %d %d %d %d %d %d %d %d \n", Amount, str, p.x, p.y, p.z, 0, 0, 0, 0, 0, 0, 1, 0);
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N - i; j++)
		{
			Point3 p = A + (I + J) * 2 / 3 + j * I + i * J;			
			writePoint(p, K);
			
			p = p + (2 * J - I) / 3;
			writePoint(p, K);
		}
}

void SnailStructure::triangle2(Point3 A, Point3 B, Point3 C, int N, int K)   //функция замощения молекулами треугольника типа 2 со "хорошей" стороной состоящей из N шестиугольников
{
	Point3 C1, P;
	C1 = C - (C - B) / (N + 2);
	triangle1(A, B, C1, N, K);                                    //замощаем правильный треугольник, далее добавим недостающие точки
	
	Point3 I = (B - A) / (N + 1);
	Point3 J = (C1 - A) / (N + 1);              //basis vectors in triangle plane

	if (N == 0)
	{
		Point3 p = C1 + J - (I - J) * (2 / 3);
		writePoint(p, K);
	}
	else
	{
		for (int i = 0; i < N + 1; i++)
		{
			int i1 = i / 2;
			Point3 p = C1 + J - (I - J) * (2 / 3) - (i1 + (i % 2)) * (2 * J - I) / 3 - i1 * (I + J) / 3;
			writePoint(p, K);
		}
	}
}

void SnailStructure::triangle3(Point3 A, Point3 B, Point3 C, int N, int K)   //замощение треугольника типа 3
{
	Point3 B1 = B - (B - A) / (N + 2);
	Point3 C1 = C - (C - B1) / (N + 2);
 	

	triangle2(A, B1, C, N, K);
	Point3 I = (A - C1) / (N + 1);
	Point3 J = (B1 -C1) / (N + 1);
	
	for (int i = 0; i < N + 2; i++)
	{
		int i1 = i / 2;
		Point3 p = B1 + J - (I + J) * 2 / 3 - (i1 + (i % 2)) * (2 * J - I) / 3 - i1 * (I + J) / 3;
		writePoint(p, K);		
	}
}

void SnailStructure::triangle2_1(Point3 A, Point3 B, Point3 C, int N, int K)   //функция замощения молекулами треугольника типа 2 со "хорошей" стороной состоящей из N шестиугольников
{
	Point3 C1 = C - (C - B) / (N+2);
	
	C1.x = C.x - (C.x - B.x) / (N + 2);
	C1.y = C.y - (C.y - B.y) / (N + 2);
	C1.z = C.z - (C.z - B.z) / (N + 2);

	Point3 I = (B -A) / (N + 1);
	Point3 J = (C1 - A) / (N + 1);      

	Point3 p = A + (I + J) / 3;
	writePoint(p, K);
	
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N - i; j++)
		{
			Point3 p = A + (I + J) * 2 / 3 + i * I + j * J;
			if (j >= i) //мостим только половину треугольника
			{     
				writePoint(p, K);
			}
			
			p = p + (2 * J - I) / 3;			
			if ((j + 1) >= i) {
				writePoint(p, K);
			}
		}


	for (int i = 0; i < N + 1; i++)
	{
		int i1 = i / 2;
		Point3 p = C1 + J - (I + J) * 2 / 3 - (i1 + (i % 2)) * (2 * J - I) / 3 - i1 * (I + J) / 3;
		writePoint(p, K);		
	}
}

void SnailStructure::triangle2_2(Point3 A, Point3 B, Point3 C, int N, int K)   //функция замощения молекулами треугольника типа 2 со "хорошей" стороной состоящей из N шестиугольников
{
	Point3 C1 = C - (C - B) / (N + 2);

	Point3 I = (B - A) / (N + 1);
	Point3 J = (C1 - A) / (N + 1);;                                                   //вводим вектора в локальной плоскости треугольника

	Point3 p = A + (I + J) / 3;
	writePoint(p, K);
		
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N - i; j++)
		{
			p = A + (I + J) * 2 / 3 + i * I + j * J;

			//мостим только половину треугольника
			if (j >= i) writePoint(p, K);
			else writePoint(p, K + 1);

			p = p + (2 * J - I) / 3;
			if ((j + 1) >= i) writePoint(p, K);
			else writePoint(p, K + 1);
		}
	}

	for (int i = 0; i < N + 1; i++)
	{
		int i1 = i / 2;
		p = C1 + J - (I + J) * 2 / 3 - (i1 + (i % 2)) * (2 * J - I) / 3 - i1 * (I + J) / 3;	
		writePoint(p, K);
	}
}

void SnailStructure::build()
{
	Point3 A, B, C, B1, C1;
	A = basePoints[0];
	B = basePoints[3];
	C = basePoints[6];
	triangle1(A, B, C, 0, 1);
	A = basePoints[1];
	B = basePoints[4];
	C = basePoints[7];
	triangle1(A, B, C, 0, 1);
	A = basePoints[2];
	B = basePoints[3];
	C = basePoints[4];
	triangle1(A, B, C, 0, 1);

	for (int j = 1; j < num - 1; j = j + 2)
	{
		int i;
		i = (j - 1) / 2;
		A = basePoints[3 * j];
		B = basePoints[3 * j + 1];
		C = basePoints[3 * j + 2];
		triangle1(A, B, C, i, i + 2);
		A = basePoints[3 * j];
		B = basePoints[3 * j + 3];
		C = basePoints[3 * j + 2];
		triangle1(A, B, C, i, i + 2);
		A = basePoints[3 * j + 4];
		B = basePoints[3 * j + 1];
		C = basePoints[3 * j + 2];
		triangle1(A, B, C, i, i + 2);
		A = basePoints[3 * j + 3];
		B = basePoints[3 * j + 2];
		C = basePoints[3 * j + 5];
		triangle2(A, B, C, i, i + 2);
		A = basePoints[3 * j + 4];
		B = basePoints[3 * j + 2];
		C = basePoints[3 * j + 5];
		triangle2(A, B, C, i, i + 2);
		A = basePoints[3 * j + 5];
		B = basePoints[3 * j + 6];
		C = basePoints[3 * j + 3];
		triangle3(A, B, C, i, i + 2);
		A = basePoints[3 * j + 5];
		B = basePoints[3 * j + 7];
		C = basePoints[3 * j + 4];
		triangle3(A, B, C, i, i + 2);
		A = basePoints[3 * j + 5];
		B = basePoints[3 * j + 6];
		C = basePoints[3 * j + 7];
		triangle1(A, B, C, i + 1, i + 2);
		A = basePoints[3 * j + 6];
		B = basePoints[3 * j + 9];
		C = basePoints[3 * j + 3];
		triangle2_2(A, B, C, i + 1, i + 2);
		A = basePoints[3 * j + 7];
		B = basePoints[3 * j + 10];
		C = basePoints[3 * j + 4];
		triangle2_2(A, B, C, i + 1, i + 2);
	}
	int j = num - 1;
	int i = (j - 1) / 2;
	A = basePoints[3 * j];
	B = basePoints[3 * j + 1];
	C = basePoints[3 * j + 2];
	triangle1(A, B, C, i, i + 2);
	A = basePoints[3 * j];
	B = basePoints[3 * j + 3];
	C = basePoints[3 * j + 2];
	triangle1(A, B, C, i, i + 2);
	A = basePoints[3 * j + 4];
	B = basePoints[3 * j + 1];
	C = basePoints[3 * j + 2];
	triangle1(A, B, C, i, i + 2);
	A = basePoints[3 * j + 3];
	B = basePoints[3 * j + 2];
	C = basePoints[3 * j + 5];
	triangle2(A, B, C, i, i + 2);
	A = basePoints[3 * j + 4];
	B = basePoints[3 * j + 2];
	C = basePoints[3 * j + 5];
	triangle2(A, B, C, i, i + 2);
	A = basePoints[3 * j + 5];
	B = basePoints[3 * j + 6];
	C = basePoints[3 * j + 3];
	triangle3(A, B, C, i, i + 2);
	A = basePoints[3 * j + 5];
	B = basePoints[3 * j + 7];
	C = basePoints[3 * j + 4];
	triangle3(A, B, C, i, i + 2);
	A = basePoints[3 * j + 5];
	B = basePoints[3 * j + 6];
	C = basePoints[3 * j + 7];
	triangle1(A, B, C, i + 1, i + 2);
	A = basePoints[3 * j + 6];
	B = basePoints[3 * j + 9];
	C = basePoints[3 * j + 3];
	triangle2_1(A, B, C, i + 1, i + 2);
	A = basePoints[3 * j + 7];
	B = basePoints[3 * j + 10];
	C = basePoints[3 * j + 4];
	triangle2_1(A, B, C, i + 1, i + 2);
}

/*Polarpoint::Polarpoint(double R,double Fi,double Psi) //конструктор
{r=R;
 fi=Fi;
 psi=Psi;}

Polarpoint::~Polarpoint()                               //деструктор
{}*/






/*double point::distance(point A, point B)                 // функция опредяляющая расстояние между точками в декартовых координатах
{ double dist;
	dist=sqrt((A.x-B.x)^2 + (A.y-B.y)^2 +(A.z-B.z)^2);
  return dist;}

double Polarpoint::Polardistance(Polarpoint A, Polarpoint B) // функция опредяляющая расстояние между точками в параметрических координатах
{ double dist;
dist=sqrt(A.r^2+B.r^2-2*A.r*B.r*cos(A.psi-B.psi)*sin(A.fi-B.fi));
  return dist;}*/



int main()
{
	int k;   //количество кусокочков спирали с углом пи/2
	cout << "enter k(number of 1/2 spitar turns)";
	cin >> k;
	k = k * 2;
	
	SnailStructure ss = SnailStructure();                 //замощение соответственных треугольников молекулами на соответствующих витках
	ss.defineBasePoints(k);
	ss.build();
	/*
	Point3 A, B, C, B1, C1;
	A = basePoints[0];
	B = basePoints[3];
	C = basePoints[6];
	triangle1(A, B, C, 0, 1);
	A = basePoints[1];
	B = basePoints[4];
	C = basePoints[7];
	triangle1(A, B, C, 0, 1);
	A = basePoints[2];
	B = basePoints[3];
	C = basePoints[4];
	triangle1(A, B, C, 0, 1);


	for (int j = 1; j < k - 1; j = j + 2)
	{
		int i;
		i = (j - 1) / 2;
		A = basePoints[3 * j];
		B = basePoints[3 * j + 1];
		C = basePoints[3 * j + 2];
		ss->triangle1(A, B, C, i, i + 2);
		A = basePoints[3 * j];
		B = basePoints[3 * j + 3];
		C = basePoints[3 * j + 2];
		ss->triangle1(A, B, C, i, i + 2);
		A = basePoints[3 * j + 4];
		B = basePoints[3 * j + 1];
		C = basePoints[3 * j + 2];
		ss->triangle1(A, B, C, i, i + 2);
		A = basePoints[3 * j + 3];
		B = basePoints[3 * j + 2];
		C = basePoints[3 * j + 5];
		ss->triangle2(A, B, C, i, i + 2);
		A = basePoints[3 * j + 4];
		B = basePoints[3 * j + 2];
		C = basePoints[3 * j + 5];
		ss->triangle2(A, B, C, i, i + 2);
		A = basePoints[3 * j + 5];
		B = basePoints[3 * j + 6];
		C = basePoints[3 * j + 3];
		ss->triangle3(A, B, C, i, i + 2);
		A = basePoints[3 * j + 5];
		B = basePoints[3 * j + 7];
		C = basePoints[3 * j + 4];
		ss->triangle3(A, B, C, i, i + 2);
		A = basePoints[3 * j + 5];
		B = basePoints[3 * j + 6];
		C = basePoints[3 * j + 7];
		ss->triangle1(A, B, C, i + 1, i + 2);
		A = basePoints[3 * j + 6];
		B = basePoints[3 * j + 9];
		C = basePoints[3 * j + 3];
		ss->triangle2_2(A, B, C, i + 1, i + 2);
		A = basePoints[3 * j + 7];
		B = basePoints[3 * j + 10];
		C = basePoints[3 * j + 4];
		ss->triangle2_2(A, B, C, i + 1, i + 2);
	}
	int j = k - 1;
	int i = (j - 1) / 2;
	A = basePoints[3 * j];
	B = basePoints[3 * j + 1];
	C = basePoints[3 * j + 2];
	ss->triangle1(A, B, C, i, i + 2);
	A = basePoints[3 * j];
	B = basePoints[3 * j + 3];
	C = basePoints[3 * j + 2];
	ss->triangle1(A, B, C, i, i + 2);
	A = basePoints[3 * j + 4];
	B = basePoints[3 * j + 1];
	C = basePoints[3 * j + 2];
	ss->triangle1(A, B, C, i, i + 2);
	A = basePoints[3 * j + 3];
	B = basePoints[3 * j + 2];
	C = basePoints[3 * j + 5];
	ss->triangle2(A, B, C, i, i + 2);
	A = basePoints[3 * j + 4];
	B = basePoints[3 * j + 2];
	C = basePoints[3 * j + 5];
	ss->triangle2(A, B, C, i, i + 2);
	A = basePoints[3 * j + 5];
	B = basePoints[3 * j + 6];
	C = basePoints[3 * j + 3];
	ss->triangle3(A, B, C, i, i + 2);
	A = basePoints[3 * j + 5];
	B = basePoints[3 * j + 7];
	C = basePoints[3 * j + 4];
	ss->triangle3(A, B, C, i, i + 2);
	A = basePoints[3 * j + 5];
	B = basePoints[3 * j + 6];
	C = basePoints[3 * j + 7];
	ss->triangle1(A, B, C, i + 1, i + 2);
	A = basePoints[3 * j + 6];
	B = basePoints[3 * j + 9];
	C = basePoints[3 * j + 3];
	ss->triangle2_1(A, B, C, i + 1, i + 2);
	A = basePoints[3 * j + 7];
	B = basePoints[3 * j + 10];
	C = basePoints[3 * j + 4];
	ss->triangle2_1(A, B, C, i + 1, i + 2);
	*/
	/*
		 int i=3;
		 int j=7;
		   A=AllCover[3*j];
		   B=AllCover[3*j+1];
		   C=AllCover[3*j+2];

		 ss->triangle1(A,B,C,i);*/
		 /* int i=1;
		  int j=1;
			A=AllCover[3*j+6];
			C1=AllCover[3*j+3];
			B1.x=AllCover[3*j+9].x-(AllCover[3*j+9].x-AllCover[3*j+3].x)/(i+2);
			B.x=AllCover[3*j+9].x+(B1.x+C1.x-2*A.x)/(3*(i+2));
			C.x=AllCover[3*j+3].x+(B1.x+C1.x-2*A.x)/(3*(i+2));

			B1.y=AllCover[3*j+9].y-(AllCover[3*j+9].y-AllCover[3*j+3].y)/(i+2);
			B.y=AllCover[3*j+9].y+(B1.y+C1.y-2*A.y)/(3*(i+2));
			C.y=AllCover[3*j+3].y+(B1.y+C1.y-2*A.y)/(3*(i+2));

			B1.z=AllCover[3*j+9].z-(AllCover[3*j+9].z-AllCover[3*j+3].z)/(i+2);
			B.z=AllCover[3*j+9].z+(B1.z+C1.z-2*A.z)/(3*(i+2));
			C.z=AllCover[3*j+3].z+(B1.z+C1.z-2*A.z)/(3*(i+2));
		  ss->triangle2(A,B,C,i); */
		  /*else if (k%4==1)
		  { for (int j=1; j<k-1; j=j+2)
			   {   int i;
				   if ((j-1)%4==0)
				   { i=(j-1)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+2];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+2];
					 B=AllCover[3*j+4];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+3];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+7];
					 C=AllCover[3*j+4];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+7];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+6];
					 B=AllCover[3*j+9];
					 C=AllCover[3*j+3];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+7];
					 B=AllCover[3*j+10];
					 C=AllCover[3*j+4];
				   ss->triangle1(A,B,C,i);
				   }
				   else
				   { i=(j-3)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+3];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+3];
				   ss->triangle3(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+7];
					 C=AllCover[3*j+4];
				   ss->triangle3(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+7];
				   ss->triangle1(A,B,C,i+1);
					 A=AllCover[3*j+6];
					 B=AllCover[3*j+9];
					 C=AllCover[3*j+3];
				   ss->triangle2(A,B,C,i+1);
					 A=AllCover[3*j+7];
					 B=AllCover[3*j+10];
					 C=AllCover[3*j+4];
				   ss->triangle2(A,B,C,i+1);
			   }
			   }
			   int j=k-1;
			   int i=j/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+3];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2_2(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2_2(A,B,C,i);
			   }
		  else if (k%4==2)
		  { for (int j=1; j<k-1; j=j+2)
			   {   int i;
				   if ((j-1)%4==0)
				   { i=(j-1)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+2];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+2];
					 B=AllCover[3*j+4];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+3];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+7];
					 C=AllCover[3*j+4];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+7];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+6];
					 B=AllCover[3*j+9];
					 C=AllCover[3*j+3];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+7];
					 B=AllCover[3*j+10];
					 C=AllCover[3*j+4];
				   ss->triangle1(A,B,C,i);
				   }
				   else
				   { i=(j-3)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+3];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+3];
				   ss->triangle3(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+7];
					 C=AllCover[3*j+4];
				   ss->triangle3(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+7];
				   ss->triangle1(A,B,C,i+1);
					 A=AllCover[3*j+6];
					 B=AllCover[3*j+9];
					 C=AllCover[3*j+3];
				   ss->triangle2(A,B,C,i+1);
					 A=AllCover[3*j+7];
					 B=AllCover[3*j+10];
					 C=AllCover[3*j+4];
				   ss->triangle2(A,B,C,i+1);
			   }
			   }
			   int j=k-1;
			   int i=(j-3)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+3];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+3];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+7];
					 C=AllCover[3*j+4];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+7];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+6];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+9];
				   ss->triangle2_1(A,B,C,i);
					 A=AllCover[3*j+7];
					 B=AllCover[3*j+4];
					 C=AllCover[3*j+10];
				   ss->triangle2_1(A,B,C,i);  }

		  else //if (k%4==3)
		  {for (int j=1; j<k-1; j=j+2)
			   {   int i;
				   if ((j-1)%4==0)
				   { i=(j-1)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+2];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+2];
					 B=AllCover[3*j+4];
					 C=AllCover[3*j+5];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+3];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+7];
					 C=AllCover[3*j+4];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+7];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+6];
					 B=AllCover[3*j+9];
					 C=AllCover[3*j+3];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+7];
					 B=AllCover[3*j+10];
					 C=AllCover[3*j+4];
				   ss->triangle1(A,B,C,i);
				   }
				   else
				   { i=(j-3)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+3];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+3];
				   ss->triangle3(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+7];
					 C=AllCover[3*j+4];
				   ss->triangle3(A,B,C,i);
					 A=AllCover[3*j+5];
					 B=AllCover[3*j+6];
					 C=AllCover[3*j+7];
				   ss->triangle1(A,B,C,i+1);
					 A=AllCover[3*j+6];
					 B=AllCover[3*j+9];
					 C=AllCover[3*j+3];
				   ss->triangle2(A,B,C,i+1);
					 A=AllCover[3*j+7];
					 B=AllCover[3*j+10];
					 C=AllCover[3*j+4];
				   ss->triangle2(A,B,C,i+1);
			   }
			   }
			   int j=k-1;
			   int i=(j-2)/4;
					 A=AllCover[3*j];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j];
					 B=AllCover[3*j+3];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+1];
					 C=AllCover[3*j+2];
				   ss->triangle1(A,B,C,i);
					 A=AllCover[3*j+3];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2_2(A,B,C,i);
					 A=AllCover[3*j+4];
					 B=AllCover[3*j+2];
					 C=AllCover[3*j+5];
				   ss->triangle2_2(A,B,C,i);
			   }


			   A.x=0;
			   A.y=0;
			   A.z=0;
			   B.x=0;
			   B.y=5;
			   B.z=0;
			   C.x=10/1.75;
			   C.y=2.5;
			   C.z=0;
			   C1.x=10;
			   C1.y=0;
			   C1.z=0;
			   B1.x=-(10/1.75+0.5);
			   B1.y=-3;
			   B1.z=0;*/


			   //ss->triangle1(A,B,C,0);    //проверка по треугольникам отдельно

				//ss->triangle2(A,B,C1,0);

				//ss->triangle3(A,B1,C1,0);
				//ss->triangle2_1(A,B,C,1);

	
	
	cout << "Amount of carbon atoms = " << Amount << "\n";

	/*  FILE*   file;                               // записываем в файл вершины икосаедрической спирали
 fopen_s(&file, "result.txt", "w");
	 for (int j=0; j<3*(k+3); j++)
	 {fprintf( file, "%15.4f  %15.4f  %15.4f  \n", AllCover[j].x,AllCover[j].y,AllCover[j].z);}
fclose(file);*/

	cout << "tops of icosaedrical spiral from angle 2pi to 2pi*(k+1), are placed in file result.txt \n Coordinates of atoms are placed in file result1.out .\n \n Delete file result1.out before recalculation \n";
	system("pause");
	return 0;
}
