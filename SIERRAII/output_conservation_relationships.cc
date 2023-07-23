 /* contient les éléments cin et cout entre autre*/
#include <iostream>
/* fstream: stream  class to both read and write  from/to  files */
#include <fstream>
/* utilisation des chaines de caractères*/
#include <string>
/* utilisation des manipulateurs*/
#include<iomanip>
/*conversion nombre string et tring nombre (stringstream,  etc)*/
#include<sstream>
#include"output_conservation_relationships.h"
using namespace std;
void output_conservation_relationships(double** matrix_coef_L, int* X_echelon_label, int& dim_Nr, int& number_vertices){
 int*    nbre_species_conservation;
 nbre_species_conservation= new int[number_vertices-dim_Nr+1];
 
string OUTPUT_FILE_NAME ="output_data/";
string  filename1 =OUTPUT_FILE_NAME+"CONSERVATION_RELATIONSHIPS.csv";
string  filename2 =OUTPUT_FILE_NAME+"matrix_coef_L.csv";
string  filename3 =OUTPUT_FILE_NAME+" X_echelon_label.csv";
string  filename4 =OUTPUT_FILE_NAME+" species_involved_in_conservation_relations.csv";
ofstream output_file1(filename1.c_str());
ofstream output_file4(filename4.c_str());
if (output_file1.is_open()){
  for ( int i=1; i<= (number_vertices-dim_Nr) ;++i){
  nbre_species_conservation[i] = 0;
     for (int j=1;j<=dim_Nr;++j)
    {
       if (matrix_coef_L[i][j] != 0){ nbre_species_conservation[i]++; if (matrix_coef_L[i][j]> 0.0){output_file1 <<" +"<<  matrix_coef_L[i][j] << "*X[" <<X_echelon_label[j] <<"] ";}
	   else{output_file1 <<  matrix_coef_L[i][j] << "*X[" <<X_echelon_label[j] <<"] ";}
	   }
     }// end  for (int j=1;j<=dim_Nr ;++j)
	 nbre_species_conservation[i]++;
 output_file1 << " -X[" <<  X_echelon_label[dim_Nr+i]<<"] = Cte"<< endl;
  }// end for ( int i=1; i<= (number_vertices-dim_Nr) ;++i)
  output_file1.close();
}else{cout << "file can not be opened"  << endl;}

if (output_file4.is_open()){
  for ( int i=1; i<= (number_vertices-dim_Nr) ;++i){
     output_file4<< nbre_species_conservation[i]<<",";
     for (int j=1;j<=dim_Nr;++j)
    {
       if (matrix_coef_L[i][j] != 0){  if (matrix_coef_L[i][j]> 0.0){output_file4  <<X_echelon_label[j]<<",";}
	   else{output_file4<<X_echelon_label[j]<<",";}
	   }
     }// end  for (int j=1;j<=dim_Nr ;++j)
 output_file4 <<  X_echelon_label[dim_Nr+i]<< endl;
  }// end for ( int i=1; i<= (number_vertices-dim_Nr) ;++i)
  output_file4.close();
}else{cout << "file can not be opened"  << endl;}

ofstream output_file2(filename2.c_str());
if (output_file2.is_open()){
  for ( int i=1; i<= (number_vertices-dim_Nr) ;++i){
     for (int j=1;j<=dim_Nr ;++j)
     {
      if (j !=dim_Nr){ output_file2 <<  matrix_coef_L[i][j] <<",";}else{ output_file2 <<  matrix_coef_L[i][j] <<endl;}
     }// end  for (int j=1;j<=nbre_reactions ;++j)
  }// end for ( int i=1; i<= nbre_reactions ;++i)
  output_file2.close();
}else{cout << "file can not be opened"  << endl;}
ofstream output_file3(filename3.c_str());
if (output_file3.is_open()){
  for ( int i=1; i<= number_vertices ;++i){
     output_file3 <<X_echelon_label[i] <<endl;
  }// end for ( int i=1; i<= number_vertices ;++i)
 output_file3.close();
}else{cout << "file can not be opened"  << endl;}
delete [] nbre_species_conservation;
return;}





