#include<cmath>
#include<iomanip>
#include <iostream>
#include <fstream>
#include"Gauss_Jordan_pivoting.h"

using namespace std;
void Gauss_Jordan_pivoting(int* X_echelon_label, int** stochio_matrix, double** stochio_matrix_echelon, int& number_vertices, int& nbre_reactions){


// pour la matrice stochiometrique, les lignes correspondent au nombre d'espèces, number_vertices, les colonnes aux nombre de reactions
int pivot;
int norow =number_vertices;
int nocolumn =  nbre_reactions;

pivot = 1;
for (int r = 1;r<=norow; r++ )
{
 if (nocolumn</*=*/ pivot){cout << "outp1"<< endl;cout << "row  = "<< r << " pivot-column = "<< pivot <<endl ;return;}
 int i=r;
 while(abs(stochio_matrix_echelon[i][pivot]) == 0)
 {
   ++i;//i = i+1;
   if (norow /*==*/< i)
   {
     i   =r;
	++pivot;// pivot = pivot+1;
	 if (nocolumn /*==*/< pivot){/*pivot = pivot -1; break;*/cout << "outp2"<< endl;cout << "row  = "<< r << " pivot-column = "<< pivot <<endl ; return;}
   } // end if (norow == i)
 } // end while(stochio_matrix[i][pivot] == 0)

   
   int exchang_label = X_echelon_label[r]  ;
   X_echelon_label[r] =  X_echelon_label[i];
   X_echelon_label[i] =exchang_label  ;
   // cout << "r label = "<<r << " va en i label = " << i<<endl; 
	 //for (int j=1;j<= norow;++j){cout <<  X_echelon_label[j]<< endl;}

 for (int j = 1; j<=nocolumn;j++)
 {
  double temp = stochio_matrix_echelon[r][j];
  stochio_matrix_echelon[r][j] =  stochio_matrix_echelon[i][j];
  stochio_matrix_echelon[i][j] = temp;
 
   
  double temp2 = stochio_matrix[r][j];
  stochio_matrix[r][j] =  stochio_matrix[i][j];
  stochio_matrix[i][j] = temp2;
 

 }
 double div = stochio_matrix_echelon[r][pivot];
 /*if (div != 0)*/ for (int j = 1; j<=nocolumn;j++){stochio_matrix_echelon[r][j] /=(div*1.0);} 
 for (int j = 1; j<=norow;j++)
 {
    if (j!=r){double sub = stochio_matrix_echelon[j][pivot];for(int k=1;k<=nocolumn; k++){stochio_matrix_echelon[j][k] -=(sub* stochio_matrix_echelon[r][k]);}}
 }
pivot = pivot+1;//if (nocolumn<= pivot){return;}else{ pivot = pivot+1;}
}// end for (int r = 1,r<=norow; r++ ) 


/*
for (int i=1;i<= number_vertices;++i){
  for (int j=1;j<= nbre_reactions;++j){
   cout << stochio_matrix[i][j]<<"  ";
  }
  cout <<""<< endl;
}
*/



return;}



















































