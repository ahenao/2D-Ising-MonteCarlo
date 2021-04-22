/*
This Program perform a single spin flip Glauber Dynamics simulation of the Ising
Magnetic system, by the technique of Monte Carlo

Andres Henao Aristizabal
Physics Engineer- Universidad Nacional de Colombia
Studying Master in Computational Physics- Universitat Politecnica Catalunya, Unversitat Barcelona
ahenaoa@unal.edu.co 
*/




#include <iostream>
#include <fstream> //Contiene la función para escribir ficheros .txt
#include <cmath>
#include <ctime>                      // define time()

#include "randomc.h"    //Library with the Random Number Generators

/*To compile put in the same directory of the .cpp the library randomc.h, mersenne.o, userintf.o
and write: 
g++ 2DIsing_principal.cpp -o 2DIsing  mersenne.o usernintf.o

The codes mersenne.cpp and userintf.cpp are constructed as .o to avoid the warnings on the differences between C and C++*/

using namespace std;

void usage ( void )
{
  cout<<"\n--------------------------------------------------------------------------------"<<endl;
  cout<<"\t\t\tMonte Carlo Simulation of\n\t\t\t     the Ising Model\n";
  cout<<"\n\t\t\t Written by Andres Henao";
  cout<<"\n\t       Master in Computational and Applied Physics\n\t\t\t    UPC-UB June 2011\n\n";
  cout<<"./2DIsing usage:\n";
  cout<<">> 2DIsing [options]\n\n";
  cout<<"Options:\n";
  cout<<"\t-d\t[integer]\tDimension\n";
  cout<<"\t-L\t[integer]\tLinear Size\n";
  cout<<"\t-T\t[real]\t\tTemperature\n";
  cout<<"\t-nmcs\t[integer]\tNumber of MonteCarlo steps\n";
  cout<<"\t-nmeas\t[integer]\tNumber of MC steps between successive measures\n";
  cout<<"\t-seed\t[integer]\tSeed of the random number generator\n";
  cout<<"\t-h  \t\t\tPrint this info\n";
  cout<<"--------------------------------------------------------------------------------\n"<<endl;
}

/********************Initialization of Variables****************/
int qq,kk;//The random spin to change
float r,p,J=1,h=0;//Random number, probability, exchange constant, external field
double Hd,delta;//To measure delta of energy
int d=2;
double L=20;
float T=2.0;
double nmcs=1e5;
int nmeas=1;
int seed=0;
double N;
int n;
double energy=0,magnet=0;

ofstream a1("Energy.txt");  //Put here the name of the energy output
ofstream a2("Magnet.txt");  //Put here the name of the Magnet output
/****************************************************************/


/******Here we parse the command line arguments; If you add an option, document it in the usage() function!****/
void read(int argc, char *argv[])  
{
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-d")) d = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-L")) L = atof(argv[++i]);
    else if (!strcmp(argv[i],"-T")) T = atof(argv[++i]);
    else if (!strcmp(argv[i],"-nmcs")) nmcs = atof(argv[++i]);
    else if (!strcmp(argv[i],"-nmeas")) nmeas = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-seed")) seed = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-h")) 
      {
      usage(); 
      exit(0);
      }
    else {
      fprintf(stderr,"Error: Command-line argument '%s' not recognized.\n",
	      argv[i]);
      exit(-1);
    }
  }
}  
/****************************************************************/
CRandomMersenne RanGen(seed);       // make instance of random number generator

/****************Initialize the Lattice***************************/
void initialize_square(int *spinij[],int n)
{
  for(int i=1;i<=n;i++)        //Initial Configuration
    {   
      for(int j=1;j<=n;j++)
	{
	  r= RanGen.Random();
	  if(r<0.5)
	    {
	      spinij[i][j]=1;
	      magnet+=1;
	    }
	  else
	    {
	      spinij[i][j]=-1;
	      magnet-=1;
	    }
	  if(i==1){spinij[n+1][j]=spinij[1][j];}//Boundary neighbours
	  if(i==n){spinij[0][j]=spinij[n][j];}
	  if(j==1){spinij[i][n+1]=spinij[i][1];}
	  if(j==n){spinij[i][0]=spinij[i][n];}
	}
    }
  for(int i = 1; i <= n; i++)       //Measuring Initial energy
    {
      for(int j = 1; j <= n; j++)
	{
	  energy=energy-h*spinij[i][j];
	  if(spinij[i][j]==spinij[i-1][j])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}
	  if(spinij[i][j]==spinij[i+1][j])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}
	  if(spinij[i][j]==spinij[i][j-1])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}
	  if(spinij[i][j]==spinij[i][j+1])
	    {energy=energy-J;}
	  else
	    {energy=energy+J;}	                     
	}
    }
  energy*=0.5;
}
/********************************************************************/


/********************Monte Carlo Update******************************/
void Monte_Carlo(int *spinij[],int n,double N)
{     
  for(int m = 1; m <= n; m++)
    {
      for(int o = 1; o <=n; o++)
	{
  //------------------Chosing a random spin to flip--------------------------//     
  //qq=RanGen.IRandom(1,n);
  //kk=RanGen.IRandom(1,n);//Vamos a cambiar el spin qq,kk
	  qq=m;
	  kk=o;
  
  //--------------------Finding the new Hamiltonian--------------------------//
  Hd=0;delta=0;
  Hd=spinij[qq-1][kk]+spinij[qq+1][kk]+spinij[qq][kk-1]+spinij[qq][kk+1];
  delta=2*Hd*spinij[qq][kk]*J + 2*h*spinij[qq][kk]; 
  
  r= RanGen.Random();
  //p=(exp(-(delta)/T));//Metropolis
  p=1.0/(1+exp((delta)/T));//Glauber Dynamics
  if(r<=p)
    {
      spinij[qq][kk]=-1*spinij[qq][kk];
      magnet+=spinij[qq][kk];
      energy+=delta;
      if(qq==1&&kk!=1&&kk!=n){spinij[n+1][kk]=-spinij[n+1][kk];}
      if(qq==n&&kk!=1&&kk!=n){spinij[0][kk]=-spinij[0][kk];}
      if(qq==1&&kk==1){spinij[n+1][n+1]=-spinij[n+1][n+1];}
      if(qq==n&&kk==n){spinij[0][0]=-spinij[0][0];}
      if(qq==1&&kk==n){spinij[n+1][0]=-spinij[n+1][0];}
      if(qq==n&&kk==1){spinij[0][n+1]=-spinij[0][n+1];}
      if(kk==1&&qq!=1&&qq!=n){spinij[qq][n+1]=-spinij[qq][n+1];}
      if(kk==n&&qq!=1&&qq!=n){spinij[qq][0]=-spinij[qq][0];}
    } 
	}
    } 
}
/*******************************************************************/


/************************Output Files***************************/
void write(int k, double NN)
{
  a1<<k<<"\t"<<"\t"<<energy<<"\t"<<energy/NN<<endl;
  a2<<k<<"\t"<<magnet<<"\t"<<magnet/NN<<endl;
}
/**************************************************************/







int main( int argc, char * argv[] )
{
  read(argc,&argv[0]);//Reads the Input
  N=pow(L,d);         
  n=int(L);

  int **spinij;	
  spinij=new int*[n+2];
  for(int i = 0; i <= n+2; i++)
    {
      spinij[i]=new int[n+2];
    }
    
  initialize_square(spinij,n); //Initialize the Lattice and measures the Initial Conditions  

  for(int k = 1; k <= nmcs; k++)
    {
      Monte_Carlo(spinij,n,N);
      if(k%nmeas==0){
	write(k,N);}
    }                       
  a1.close();
  a2.close();
  
  //*****Release the memory*********//
  for(int i = 0; i <= n; i++){
    delete [] spinij[i];}
  delete [] spinij;
  //********************************//
  
  return 0;
}


