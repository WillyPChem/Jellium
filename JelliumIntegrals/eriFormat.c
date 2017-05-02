//   JPIS.c       Physical Chemistry II Project
#include<stdio.h>
#include<cmath>
#include<math.h>
#include<stdlib.h>
#include<malloc.h>
#include<complex.h>
#include<time.h>
#include<string.h>
#include <iostream>
#include <fstream>



int main() {

FILE *EEfp, *stream;
int i,j,k,l;
double val;
string *list;
int n;


int *a, *b, *c, *d;
double *value;
int index = 0;

	string list[1230];

	ifstream inFile;
	inFile.open("ERI.dat");

	while (!inFile.eof()) {

		inFile >> list[n];
		n++;

	}

/*

	a = (int *)malloc(26*26*sizeof(int)); // a corresponding to phi
    b = (int *)malloc(26*26*sizeof(int)); // b corresponding to phi
    c = (int *)malloc(26*26*sizeof(int)); // c corresponding to phi
    d = (int *)malloc(26*26*sizeof(int)); // d corresponding to phi
    value = (double *)malloc(26*26*sizeof(double)); // value energy corresponding

	stream = fopen("ERI.dat", "r");


	std::ifstream file("ERI.dat");
	std::string str;


	for(int index=0; index<10; index++) {
		
			fscanf(stream, "%d\n", &i);
			fscanf(stream, "%d\n", &j);
			fscanf(stream, "%d\n", &k);
			fscanf(stream, "%d\n", &l);
			fscanf(stream, "%lf\n", &val);

			a[index] = i;
			b[index] = j;
			c[index] = k;
			d[index] = l;
			value[index] = val;


	}

*/
	for(int r=0; r<5; r++)
	{

	//printf(" %d %d %d %d %lf\n",a[index],b[index],c[index],d[index],value[index]);

		printf(" %s\n ",list[r]);

	}


}