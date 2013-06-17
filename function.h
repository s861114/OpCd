#ifndef FUNCTION_H
#define FUNCTION_H

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <cmath>


using namespace std;

template<class T>
class umat
{
	public:
		T* Tmatrix1(int i);
		T** Tmatrix2(int i,int j);
		T*** Tmatrix3(int i,int j,int k);
		T**** Tmatrix4(int i,int j,int k,int l);
		void Tfree1(T* a);
		void Tfree2(T** a);
		void Tfree3(T*** a);
		void Tfree4(T**** a);
};

template<class T>
T* umat<T>::Tmatrix1(int i)
{
	return (T*)malloc(sizeof(T)*i);
}
template<class T>
T** umat<T>::Tmatrix2(int i,int j)
{
	T** a = (T**)malloc(sizeof(T*)*i);
	a[0] = (T*)malloc(sizeof(T)*i*j);
	for(int p=1;p<i; p++)
		a[p] = a[p-1]+j;
	return a;
}
template<class T>
T*** umat<T>::Tmatrix3(int i,int j,int k)
{
	T*** a = (T***)malloc(sizeof(T**)*i);
	if(a==NULL)
	{
		printf("nnnnn0\n");
		exit(1);
	}
	a[0] = (T**)malloc(sizeof(T*)*i*j);
	if(a[0]==NULL)
	{
		printf("nnnnn1\n");
		exit(1);
	}
	a[0][0] = (T*)malloc(sizeof(T)*i*j*k);
	if(a[0][0]==NULL)
	{
		printf("nnnnn2\n");
		exit(1);
	}
	for(int p=0;p<i; p++)
	{
		a[p] = a[0]+p*j;
		for(int q=0; q<j; q++)
			a[p][q] = a[0][0]+p*j*k+q*k;
	}
	return a;
}
template<class T>
T**** umat<T>::Tmatrix4(int i,int j,int k,int l)
{
	T**** a = (T****)malloc(sizeof(T***)*i);
	a[0] = (T***)malloc(sizeof(T**)*i*j);
	a[0][0] = (T**)malloc(sizeof(T*)*i*j*k);
	a[0][0][0]= (T*)malloc(sizeof(T)*i*j*k*l);
	for(int p=0;p<i; p++)
	{
		a[p] = a[0]+p*j;
		for(int q=0; q<j; q++)
		{
			a[p][q] = a[0][0]+p*j*k+q*k;
			for(int r=0; r<k; r++)
			{
				a[p][q][r] = a[0][0][0] + p*j*k*l + q*k*l + r*l;
			}
		}
	}
	return a;
}
template<class T>
void umat<T>::Tfree1(T* a)
{
	free(a);
}
template<class T>
void umat<T>::Tfree2(T** a)
{
	free(a[0]);
	free(a);
}
template<class T>
void umat<T>::Tfree3(T*** a)
{
	free(a[0][0]);
	free(a[0]);
	free(a);
}
template<class T>
void umat<T>::Tfree4(T**** a)
{
	free(a[0][0][0]);
	free(a[0][0]);
	free(a[0]);
	free(a);
}

#endif
