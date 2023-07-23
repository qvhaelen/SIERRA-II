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
#include"output_matrices_fluxes.h"
using namespace std;
void  output_matrices_fluxes(double** matrix_GAMMA,  double** matrix_C, int& nbre_reactions,int& number_vertices, int& time_point,int* X_echelon_label){
string OUTPUT_FILE_NAME ="output_data/screenshot/";
std::ostringstream tp_shot;
tp_shot << time_point;
std::string label_name = tp_shot.str();
string  filename1 =OUTPUT_FILE_NAME+label_name+ "-MATRIX-C.csv";
string  filename2 =OUTPUT_FILE_NAME+label_name+ "-MATRIX-GAMMA.csv";
ofstream output_file1(filename1.c_str());
if (output_file1.is_open()){
  for ( int i=1; i<= nbre_reactions ;++i){
     for (int j=1;j<=nbre_reactions ;++j)
     {
      if (j != nbre_reactions){ output_file1 <<  matrix_C[i][j] <<",";}else{ output_file1 <<  matrix_C[i][j] <<endl;}
     }// end  for (int j=1;j<=nbre_reactions ;++j)
 }// end for ( int i=1; i<= nbre_reactions ;++i)
  output_file1.close();
}else{cout << "file can not be opened"  << endl;}
ofstream output_file2(filename2.c_str());
if (output_file2.is_open()){
  for ( int i=1; i<= number_vertices ;++i){
     for (int j=1;j<=nbre_reactions ;++j)
     {
      if (j != nbre_reactions){ output_file2 <<  matrix_GAMMA[X_echelon_label[i]][j] <<",";}else{ output_file2 << matrix_GAMMA[X_echelon_label[i]][j] <<endl;}
     }// end  for (int j=1;j<=nbre_reactions ;++j)
  }// end for ( int i=1; i<= number_vertices ;++i)
  output_file2.close();
}else{cout << "file can not be opened"  << endl;}
return;}





