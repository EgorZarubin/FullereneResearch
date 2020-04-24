

#include <iostream>
using namespace std;
#include <math.h>
#include <iomanip>
#include <fstream>
#include <conio.h>
#include <string.h>
#include <string>

 double const pi=3.14;
 double const l=1.4; 
 int Amount=0;
 double coeff=1.647278;//1.14;
 double coeff2=0;
 char str[2]="C";


class point // точка в декартовых координатах
{ public:
   double x,y,z;
   point(){};
   ~point(){};
  //void TransToPolar (x,y,z); //функция преобразования в полярные координаты
   double distance(point A, point B); // функция опредяляющая расстояние между точками в декартовых координатах
};


class Polarpoint // точка в параметризованных координатах
{public:
	double r,fi,psi;
	int k;
	
	point IcoV[3]; // массив для записи трех точек на одном четверть-витке

	Polarpoint(){};
	~Polarpoint(){};

  point TransToDecart(Polarpoint B); //функция преобразования в декартовы координаты
  double Polardistance(Polarpoint A, Polarpoint B); // функция опредяляющая расстояние между точками в параметрических координатах
  void build(int k);
};

 
class structure
{ point A,B,C; // вершины треугольника
  int N; // к-во шестиугольников в наименьшей хорошей стороне
public:
		void triangle1(point A, point B, point C, int N);
	    void triangle2(point A, point B, point C, int N);
		void triangle3(point A, point B, point C, int N);
		void triangle2_1(point A, point B, point C, int N);
		void triangle2_2(point A, point B, point C, int N);
	
};

void structure::triangle1(point A, point B, point C, int N)   //функция замощения молекулами правильного треугольника со стороной состоящей из N шестиугольников
{point I,J;                                                   //вводим вектора в локальной плоскости треугольника
 double x,y,z;
	
	I.x=(-A.x+B.x)/(N+1);
	I.y=(-A.y+B.y)/(N+1);
	I.z=(-A.z+B.z)/(N+1);


    J.x=(-A.x+C.x)/(N+1);
	J.y=(-A.y+C.y)/(N+1);
	J.z=(-A.z+C.z)/(N+1);

FILE*   file;                                                 //координаты молекул не хранятся в массиве, а сразу записываются в файл
 fopen_s(&file, "result1.txt", "a");
	for (int i=0; i<N+1; i++)
	      {x=A.x+(I.x+J.x)/3+i*I.x;
	       y=A.y+(I.y+J.y)/3+i*I.y;
		   z=A.z+(I.z+J.z)/3+i*I.z;
		   Amount++;
		   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}
	for (int i=0; i<N; i++)
		for (int j=0; j<N-i; j++)
		  {x=A.x+(I.x+J.x)*2/3+j*I.x+i*J.x;
           y=A.y+(I.y+J.y)*2/3+j*I.y+i*J.y;
		   z=A.z+(I.z+J.z)*2/3+j*I.z+i*J.z;
		   Amount++;
		   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);
	       x=x+(2*J.x-I.x)/3;
		   y=y+(2*J.y-I.y)/3;
		   z=z+(2*J.z-I.z)/3;
		   Amount++;
		   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}
fclose(file);}

void structure::triangle2(point A, point B, point C, int N)   //функция замощения молекулами треугольника типа 2 со "хорошей" стороной состоящей из N шестиугольников
{point C1;
 C1.x=C.x-(C.x-B.x)/(N+2);
 C1.y=C.y-(C.y-B.y)/(N+2);
 C1.z=C.z-(C.z-B.z)/(N+2);
 triangle1(A,B,C1,N);                                    //замощаем правильный треугольник, далее добавим недостающие точки
 point I,J;
 double x,y,z;
	
	I.x=(-A.x+B.x)/(N+1);
	I.y=(-A.y+B.y)/(N+1);
	I.z=(-A.z+B.z)/(N+1);


    J.x=(-A.x+C1.x)/(N+1);
	J.y=(-A.y+C1.y)/(N+1);
	J.z=(-A.z+C1.z)/(N+1);

	
	FILE*   file;                                       //дописываем в файл недостающие точки
 fopen_s(&file, "result1.txt", "a");
    if (N==0)
	{   x=C1.x+J.x-(I.x+J.x)*2/3;
		y=C1.y+J.y-(I.y+J.y)*2/3;
		z=C1.z+J.z-(I.z+J.z)*2/3;
		Amount++;
	    fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}
	else
	{
	for (int i=0; i<N+1; i++)
	{int i1=i/2;
	 x=C1.x+J.x-(I.x+J.x)*2/3-(i1+(i%2))*(2*J.x-I.x)/3-i1*(I.x+J.x)/3;
     y=C1.y+J.y-(I.y+J.y)*2/3-(i1+(i%2))*(2*J.y-I.y)/3-i1*(I.y+J.y)/3;
     z=C1.z+J.z-(I.z+J.z)*2/3-(i1+(i%2))*(2*J.z-I.z)/3-i1*(I.z+J.z)/3;
     Amount++;
	 fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}}
fclose(file);

}
void structure::triangle3(point A, point B, point C, int N)   //замощение треугольника типа 3
{point C1,B1;

 C1.x=C.x-(C.x-B.x)/(N+2);
 C1.y=C.y-(C.y-B.y)/(N+2);
 C1.z=C.z-(C.z-B.z)/(N+2);

 B1.x=B.x-(B.x-A.x)/(N+2);
 B1.y=B.y-(B.y-A.y)/(N+2);
 B1.z=B.z-(B.z-A.z)/(N+2);
 triangle2(A,B1,C,N);
 point I,J;
 double x,y,z;
	
	I.x=(-C1.x+A.x)/(N+1);
	I.y=(-C1.y+A.y)/(N+1);
	I.z=(-C1.z+A.z)/(N+1);


    J.x=(-C1.x+B1.x)/(N+1);
	J.y=(-C1.y+B1.y)/(N+1);
	J.z=(-C1.z+B1.z)/(N+1);

	
	FILE*   file; 
 fopen_s(&file, "result1.txt", "a");

 for (int i=0; i<N+2; i++)
 {int i1=i/2;
	 x=B1.x+J.x-(I.x+J.x)*2/3-(i1+(i%2))*(2*J.x-I.x)/3-i1*(I.x+J.x)/3;
     y=B1.y+J.y-(I.y+J.y)*2/3-(i1+(i%2))*(2*J.y-I.y)/3-i1*(I.y+J.y)/3;
     z=B1.z+J.z-(I.z+J.z)*2/3-(i1+(i%2))*(2*J.z-I.z)/3-i1*(I.z+J.z)/3;
     Amount++;
	 fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}
fclose(file);

}


void structure::triangle2_1(point A, point B, point C, int N)   //функция замощения молекулами треугольника типа 2 со "хорошей" стороной состоящей из N шестиугольников
{point C1;
 C1.x=C.x-(C.x-B.x)/(N+2);
 C1.y=C.y-(C.y-B.y)/(N+2);
 C1.z=C.z-(C.z-B.z)/(N+2);

 point I,J;                                                   //вводим вектора в локальной плоскости треугольника
 double x,y,z;
	
	I.x=(-A.x+B.x)/(N+1);
	I.y=(-A.y+B.y)/(N+1);
	I.z=(-A.z+B.z)/(N+1);


    J.x=(-A.x+C1.x)/(N+1);
	J.y=(-A.y+C1.y)/(N+1);
	J.z=(-A.z+C1.z)/(N+1);

FILE*   file;                                                 //координаты молекул не хранятся в массиве, а сразу записываются в файл
 fopen_s(&file, "result1.txt", "a");
	
	       x=A.x+(I.x+J.x)/3;
	       y=A.y+(I.y+J.y)/3;
		   z=A.z+(I.z+J.z)/3;
		   Amount++;
		   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);
	for (int i=0; i<N; i++)
		for (int j=0; j<N-i; j++)
		  {x=A.x+(I.x+J.x)*2/3+i*I.x+j*J.x;
           y=A.y+(I.y+J.y)*2/3+i*I.y+j*J.y;
		   z=A.z+(I.z+J.z)*2/3+i*I.z+j*J.z;
		  
		   if(j>=i) {                                                     //мостим только половину треугольника
			   Amount++;
			   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}
	       x=x+(2*J.x-I.x)/3;
		   y=y+(2*J.y-I.y)/3;
		   z=z+(2*J.z-I.z)/3;
		   if((j+1)>=i) {
			   Amount++;
			   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}}
                                 

 for (int i=0; i<N+1; i++)
 {int i1=i/2;
	 x=C1.x+J.x-(I.x+J.x)*2/3-(i1+(i%2))*(2*J.x-I.x)/3-i1*(I.x+J.x)/3;
     y=C1.y+J.y-(I.y+J.y)*2/3-(i1+(i%2))*(2*J.y-I.y)/3-i1*(I.y+J.y)/3;
     z=C1.z+J.z-(I.z+J.z)*2/3-(i1+(i%2))*(2*J.z-I.z)/3-i1*(I.z+J.z)/3;
     Amount++;
	 fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}
fclose(file);

}

void structure::triangle2_2(point A, point B, point C, int N)   //функция замощения молекулами треугольника типа 2 со "хорошей" стороной состоящей из N шестиугольников
{point C1;
 C1.x=C.x-(C.x-B.x)/(N+2);
 C1.y=C.y-(C.y-B.y)/(N+2);
 C1.z=C.z-(C.z-B.z)/(N+2);

 point I,J;                                                   //вводим вектора в локальной плоскости треугольника
 double x,y,z;
	
	I.x=(-A.x+B.x)/(N+1);
	I.y=(-A.y+B.y)/(N+1);
	I.z=(-A.z+B.z)/(N+1);


    J.x=(-A.x+C1.x)/(N+1);
	J.y=(-A.y+C1.y)/(N+1);
	J.z=(-A.z+C1.z)/(N+1);

FILE*   file;                                                 //координаты молекул не хранятся в массиве, а сразу записываются в файл
 fopen_s(&file, "result1.txt", "a");
	
	       x=A.x+(I.x+J.x)/3;
	       y=A.y+(I.y+J.y)/3;
		   z=A.z+(I.z+J.z)/3;
		   Amount++;
		   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);
	for (int i=0; i<N; i++)
		for (int j=0; j<N-i; j++)
		  {x=A.x+(I.x+J.x)*2/3+j*I.x+i*J.x;
           y=A.y+(I.y+J.y)*2/3+j*I.y+i*J.y;
		   z=A.z+(I.z+J.z)*2/3+j*I.z+i*J.z;
		  
		   if(i>=j) {                                                     //мостим только половину треугольника
			   Amount++;
			   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}
	       x=x+(2*J.x-I.x)/3;
		   y=y+(2*J.y-I.y)/3;
		   z=z+(2*J.z-I.z)/3;
		   if((i+1)>=j) {
			   Amount++;
			   fprintf( file, "%15.4f  %15.4f  %15.4f \n", x,y,z);}}
                                 

}
/*Polarpoint::Polarpoint(double R,double Fi,double Psi) //конструктор
{r=R;
 fi=Fi;
 psi=Psi;}

Polarpoint::~Polarpoint()                               //деструктор
{}*/

point Polarpoint::TransToDecart(Polarpoint B)
{ point A;
  A.x=B.r*sin(B.psi)*cos(B.fi);
  A.y=B.r*cos(B.psi);
  A.z=B.r*sin(B.psi)*sin(B.fi);
  return A; 
}



void Polarpoint::build(int k)     //функция, позволяющая найти вершивы икосаедрической спирали, ее вершины являются вершинами треугольников типа 1,2,3
{ Polarpoint PIcoV[3];            //вершины задаются аналогично стандартным вершинам икосаедра
  //double x,y,z;
	
if (k%2==0) 
	{PIcoV[0].fi=pi*k/2;
	PIcoV[0].psi=atan(1.6);
	PIcoV[0].r=k*coeff+coeff2;
	
	PIcoV[1].fi=pi*k/2;
	PIcoV[1].psi=pi-atan(1.618);
	PIcoV[1].r=k*coeff+coeff2;
	
	PIcoV[2].fi=pi*k/2+atan(1.618);
	PIcoV[2].psi=pi/2;
	PIcoV[2].r=k*coeff+coeff2;}
else if (k!=1)
   {PIcoV[0].fi=pi*k/2;
	PIcoV[0].psi=atan(1/1.618);
	PIcoV[0].r=(k-1)*coeff+coeff2;
	
	PIcoV[1].fi=pi*k/2;
	PIcoV[1].psi=pi-atan(1/1.618);
	PIcoV[1].r=(k-1)*coeff+coeff2;
	
	PIcoV[2].fi=pi*k/2+atan(1/1.618);
	PIcoV[2].psi=pi/2;
	PIcoV[2].r=(k+1)*coeff+coeff2;}
else if (k==1)
    {PIcoV[0].fi=pi*k/2;
	PIcoV[0].psi=atan(1/1.618);
	PIcoV[0].r=2*coeff+coeff2;
	
	PIcoV[1].fi=pi*k/2;
	PIcoV[1].psi=pi-atan(1/1.618);
	PIcoV[1].r=2*coeff+coeff2;
	
	PIcoV[2].fi=pi*k/2+atan(1/1.618);
	PIcoV[2].psi=pi/2;
	PIcoV[2].r=2*coeff+coeff2;}

for (int i=0; i<3; i++)
{   
	   IcoV[i].x=PIcoV[i].r*sin(PIcoV[i].psi)*cos(PIcoV[i].fi);
       IcoV[i].y=PIcoV[i].r*cos(PIcoV[i].psi);
       IcoV[i].z=PIcoV[i].r*sin(PIcoV[i].psi)*sin(PIcoV[i].fi);
	//IcoV[i]=TransToDecart(PIcoV[i]);

}	
	}





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
    cout<<"enter k(number of 1/2 spitar turns)";
    cin>>k;
	k=k*2;
	point *AllCover=new point[3*(k+5)];            //массив вершин спирального икосаедра
  
	 for (int i=1; i<k+6; i++)                 // вычисление координат вершин, и записи их в массив начиная с угла 2пи
        {Polarpoint* kk = new Polarpoint();
	     kk->build(i);
		 for (int j=0; j<3; j++)
		 {AllCover[3*(i-1)+j]=kk->IcoV[j];}
	     delete kk;
		   }
  FILE*   file;                                                 //координаты молекул не хранятся в массиве, а сразу записываются в файл
 fopen_s(&file, "result1.txt", "w");   
 fprintf( file, "%s \n", "POSITION(Astrm)");	
 fclose(file);
	 structure* ss = new structure();                 //замощение соответственных треугольников молекулами на соответствующих витках
	 point A,B,C,B1,C1;
	       A=AllCover[0];
		   B=AllCover[3];
		   C=AllCover[6];
	 ss->triangle1(A,B,C,0);
	       A=AllCover[1];
		   B=AllCover[4];
		   C=AllCover[7];
	 ss->triangle1(A,B,C,0);
	       A=AllCover[2];
		   B=AllCover[3];
		   C=AllCover[4];
	 ss->triangle1(A,B,C,0);
	 
	 
 for (int j=1; j<k-1; j=j+2)
	 {   int i;
		   i=(j-1)/2;
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
     int j=k-1;
	 int i=(j-1)/2;
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
		 ss->triangle2_1(A,B,C,i+1);
		   A=AllCover[3*j+7];
		   B=AllCover[3*j+10];
		   C=AllCover[3*j+4];
		 ss->triangle2_1(A,B,C,i+1);    
	
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
	 }*/


     A.x=0;
	 A.y=0;
	 A.z=0;
	 B.x=0;
	 B.y=5;
	 B.z=0;
	 C.x=-5*1.73/2;
	 C.y=2.5;
	 C.z=0;
	 C1.x=10;
	 C1.y=5;
	 C1.z=0;
	 B1.x=0;
	 B1.y=-10;
	 B1.z=0;


	//ss->triangle1(A,B,C,0);    //проверка по треугольникам отдельно
	
	//ss->triangle2(B,A,C1,0);
	 
	//ss->triangle3(A,C1,B1,0);
	 //ss->triangle2_1(A,B,C,1);
   
	delete ss;
  cout<<"Amount of carbon atoms = "<<Amount<<"\n";
	  
	/*  FILE*   file;                               // записываем в файл вершины икосаедрической спирали
 fopen_s(&file, "result.txt", "w");
	 for (int j=0; j<3*(k+3); j++)
	 {fprintf( file, "%15.4f  %15.4f  %15.4f  \n", AllCover[j].x,AllCover[j].y,AllCover[j].z);}
fclose(file);*/
		 
delete [] AllCover;
cout<<"tops of icosaedrical spiral from angle 2pi to 2pi*(k+1), are placed in file result.txt \n Coordinates of atoms are placed in file result1.txt .\n \n Delete file result1.txt before recalculation \n";
  system("pause");
  return 0;
}
	