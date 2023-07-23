//**************************************************************************************************
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include<iomanip>
#include<sstream>
#include<cmath>
#include<numeric>
#include<algorithm>
#include<stack>
#include<queue>
#include<set>
#include<assert.h>

#include"linalg.h"
#include"ap.h"
#include"stdafx.h"

//**************************************************************************************************
#include "Initial_set_up.h"
#include"conversion_reaction_index.h"
#include"Stochiometric_matrix_set_up.h"
#include"Fluxes_rates_analytical.h"
#include "Fluxes_rates_analytical_init.h"
#include"Initialization_kinetic_constant.h"
#include"output_matrices_fluxes.h"
#include"output_conservation_relationships.h"
#include"output_statistics_sierraII.h"
#include"Gauss_Jordan_pivoting.h"
using namespace std;


string INPUT_FILE_NAME = "input_data/";
string OUTPUT_FILE_NAME = "output_data/";
string TIME_SERIE_DATA = INPUT_FILE_NAME + "0800-LIFR-GP130-LOW-R.txt";
string SCREENSHOT_DATA = INPUT_FILE_NAME + "LIST-SCREENSHOT-ID.txt";


int main(){
clock_t start, end=0; 
assert( (start = clock()) != -1);
cout <<""<< endl; 
cout <<"                            ****************************************************                "<< endl;
cout <<"                            *                                                  *           "<< endl;
cout <<"                            *                  CODE SIERRA II                  *        "<< endl;
cout <<"                            *                                                  *     "<< endl;
cout <<"                            ****************************************************           "<< endl;
cout <<"                  "<< endl;
cout <<"========================================================================================================================="<< endl;
cout <<"                                  RESUME DES PRINCIPALES CARACTERISTIQUES:             "<< endl;
cout <<"========================================================================================================================="<< endl;
cout << "                   MCA COMPUTATION OF CONSERVATION RELATIONSHIPS AND FLUXES" << endl;
cout <<"========================================================================================================================="<< endl;
cout << ""<< endl; 

int  number_vertices; 
int  nbre_colonne;
int nbre_reactions; 
int  number_time_points;
double two_pi = 2.0*3.14159265358979323846;
int nbre_kinetic_cte; 
Initial_set_up(number_vertices, nbre_colonne,number_time_points, nbre_reactions,nbre_kinetic_cte);
bool MCA = true;

if (MCA){
/////////////////////////////////////////////////////////////////////////////////////////////////////
// MCA COMPUTATION OF CONSERVATION RELATIONSHIPS AND FLUXES
/////////////////////////////////////////////////////////////////////////////////////////////////////
   int dim_Nr = 0;
   int dim_Nr_SVD = 0;
   double*    X;
   X= new double[number_vertices+1];
   for (int i=1; i<=number_vertices;++i){X[i] = 0;}
   
   int*    X_echelon_label;
   X_echelon_label= new int[number_vertices+1];
   for (int i=1; i<=number_vertices;++i){X_echelon_label[i] = i;}
   
   int*    X_echelon_label_test;
   X_echelon_label_test= new int[number_vertices+1];
   for (int i=1; i<=number_vertices;++i){X_echelon_label_test[i] = i;}
   
   double*    base_row;
   base_row= new double[number_vertices+1];
   for (int i=1; i<=number_vertices;++i){base_row[i] = 0;}
   
   double*    k;
   k= new double[nbre_kinetic_cte+1];

   double** DvDx_matrix=NULL;
   DvDx_matrix = new double*[nbre_reactions+1];
   for (int i=1;i<= nbre_reactions; ++i){DvDx_matrix[i] = new double[number_vertices+1];}
 
     for (int i=1;i<= nbre_reactions;++i){
      for (int j=1;j<= number_vertices;++j){
       DvDx_matrix[i][j] =0;
      }
     }
   
   
   int** stochio_matrix=NULL;
   stochio_matrix = new int*[number_vertices+1];
   for (int i=1;i<= number_vertices; ++i){stochio_matrix[i] = new int[nbre_reactions+1];}
   
   for (int i=1;i<= number_vertices;++i){
    for (int j=1;j<= nbre_reactions;++j){
     stochio_matrix[i][j] =0;
    }
   }
   double** stochio_matrix_echelon=NULL;
   stochio_matrix_echelon = new double*[number_vertices+1];
   for (int i=1;i<= number_vertices; ++i){stochio_matrix_echelon[i] = new double[nbre_reactions+1];}



Initialization_kinetic_constant(k);
Fluxes_rates_analytical_init(DvDx_matrix,X,k);
Stochiometric_matrix_set_up(stochio_matrix);

int value_row; 
for (int i=1;i<= number_vertices;++i){
value_row = 0;
  for (int j=1;j<= nbre_reactions;++j){
   value_row =value_row +  abs(stochio_matrix[i][j]);
  
  }
  
  if (value_row == 0){cout << "ligne singuliere, nrow  =" << i<< endl;}
}


for (int i=1;i<= number_vertices;++i){
  for (int j=1;j<= nbre_reactions;++j){
   stochio_matrix_echelon[i][j] =  stochio_matrix[i][j];
  }
} 


 Gauss_Jordan_pivoting(X_echelon_label, stochio_matrix, stochio_matrix_echelon, number_vertices,  nbre_reactions);
 


for (int i=1;i<= number_vertices;++i){
  for (int j=1;j<= nbre_reactions;++j){
	base_row[i] = base_row[i] +abs(stochio_matrix_echelon[i][j]);
  }
  if (base_row[i] !=0){dim_Nr++;}
}
cout << ""<< endl;
cout << " Nbre de  vecteurs independants suivant le resultat de Gauss pivoting =  "<<dim_Nr << endl;
cout << ""<< endl;

 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 alglib::real_2d_array B_test;

alglib::real_2d_array V_test;
alglib::real_2d_array U_test;
alglib::real_1d_array Singular_test;
B_test.setlength(number_vertices,nbre_reactions);

  for (int i=0;i<number_vertices;++i){
    for (int j=0; j< nbre_reactions; ++j){
      B_test[i][j] =stochio_matrix[i+1][j+1];
    }
  }
alglib::rmatrixsvd(B_test, number_vertices,  nbre_reactions, 2, 2, 2, Singular_test, U_test, V_test);
double *a_row_test = Singular_test.getcontent();
for (int i=0;  i<number_vertices;  ++i){if (a_row_test[i] != 0.0){dim_Nr_SVD++;}}
cout << " Nbre de  vecteurs independants suivant le resultat de la SVD =  "<< dim_Nr_SVD << endl;	


if (dim_Nr != dim_Nr_SVD){cout<< "PROBLEM WITH GAUSS AND RELATED SVD RESULT.... PROCEDURE FAILURE EMERGENCY STOP "<< endl;return 0;}

double** stochio_matrix_Nr=NULL;
stochio_matrix_Nr = new double*[dim_Nr+1];
for (int i=1;i<= dim_Nr; ++i){stochio_matrix_Nr[i] = new double[nbre_reactions+1];}

double** stochio_matrix_N_dependant=NULL;
stochio_matrix_N_dependant = new double*[number_vertices-dim_Nr+1];
for (int i=1;i<= (number_vertices-dim_Nr); ++i){stochio_matrix_N_dependant[i] = new double[nbre_reactions+1];}

double** matrix_coef_L=NULL;
matrix_coef_L = new double*[number_vertices-dim_Nr+1];
for (int i=1;i<= (number_vertices-dim_Nr); ++i){matrix_coef_L[i] = new double[dim_Nr+1];}

double** matrix_coef_L_complet=NULL;
matrix_coef_L_complet = new double*[number_vertices+1];
for (int i=1;i<= number_vertices; ++i){matrix_coef_L_complet[i] = new double[dim_Nr+1];}

double** Bit=NULL;
Bit = new double*[nbre_reactions];
for (int i=0;i< nbre_reactions; ++i){Bit[i] = new double[dim_Nr];}

double** Nr_pseudoinverse=NULL;
Nr_pseudoinverse = new double*[nbre_reactions];
for (int i=0;i< nbre_reactions; ++i){Nr_pseudoinverse[i] = new double[dim_Nr];}

int new_l = 0;
int new_l_d = 0;
for (int i=1;i<= number_vertices;++i){
  if(i<= dim_Nr){
  new_l++;
    for (int j=1;j<= nbre_reactions;++j){
      stochio_matrix_Nr[new_l][j] =  stochio_matrix[i][j];
    }
 }else{
   new_l_d++;
    for (int j=1;j<= nbre_reactions;++j){
      stochio_matrix_N_dependant[new_l_d][j] =  stochio_matrix[i][j];
    }
} 
}
///////////////////CALCUL SVD ET PSEUDO INVERSE POUR  stochio_matrix_Nr ////////////////////////////////////

alglib::real_2d_array B_Nr;

alglib::real_2d_array V;
alglib::real_2d_array U;
alglib::real_1d_array Singular;
B_Nr.setlength(dim_Nr,nbre_reactions);

  for (int i=0;i<dim_Nr;++i){
    for (int j=0; j< nbre_reactions; ++j){
      B_Nr[i][j] =stochio_matrix_Nr[i+1][j+1];
    }
  }
alglib::rmatrixsvd(B_Nr, dim_Nr,  nbre_reactions, 2, 2, 2, Singular, U, V);
double *a_row = Singular.getcontent();


for (int i=0;  i<nbre_reactions;  ++i){
  for (int j=0; j< dim_Nr; ++j){
    Bit[i][j] = 0;
	
   if(a_row[i] !=0 && i< dim_Nr ){ Bit[i][j] = (1/a_row[i])*U[j][i];}else{ Bit[i][j] = 0;}
  }
}

for (int i=0;  i<nbre_reactions;  ++i){
  for (int j=0; j<dim_Nr; ++j){
  Nr_pseudoinverse[i][j] = 0;
   for (int k=0; k<nbre_reactions; ++k){
     Nr_pseudoinverse[i][j] += V[k][i]*Bit[k][j];
   }  
  }
}

/////////////////////////////////////////////// CALCUL COEFFICIENTS L ///////////////////

for (int i=0;  i<(number_vertices-dim_Nr);  ++i){
  for (int j=0; j<dim_Nr; ++j){
  matrix_coef_L[i+1][j+1] = 0;
   for (int k=0; k<nbre_reactions; ++k){
    matrix_coef_L[i+1][j+1] +=stochio_matrix_N_dependant[i+1][k+1]*Nr_pseudoinverse[k][j];
   }  
  }
}
 for (int i=1; i<=dim_Nr; ++i){
  for (int j=1;j<=dim_Nr;++j){
     if(i == j){   matrix_coef_L_complet[i][j]  = 1;}else{  matrix_coef_L_complet[i][j]  = 0;}
  }
 }

  for (int i=1; i<=(number_vertices-dim_Nr); ++i){
  for (int j=1;j<=dim_Nr;++j){
   if (  abs(matrix_coef_L[i][j]) <= 10e-10){matrix_coef_L[i][j]  = 0;}
   matrix_coef_L_complet[i+dim_Nr][j] = matrix_coef_L[i][j];
  }
}

output_conservation_relationships(matrix_coef_L, X_echelon_label, dim_Nr, number_vertices);  
//////////////////////////////////////////////////////////////////////////////////////////  
/////////////////  DEBUT DE LA BOUCLE TEMPORELLE POUR CALCUL DES COEFFICIENTS DE FLUX/////
//////////////////////////////////////////////////////////////////////////////////////////


//////////////////////  CHARGEMENT DES FICHIERS DE TIME SERIES ET LES SCREENSHOTS ////////

  
   string line_timeserie, tempnum_timeserie;
   int end_line_timeserie = 0;
   vector<double>  input_vector_timeserie;
   vector<double>::iterator  input_vector_timeserie_it;
   double brol_timeserie;
   ifstream name_file(  TIME_SERIE_DATA.c_str(), ios::in);
 if (name_file.is_open()){
    while(!name_file.eof()){
       getline(name_file, line_timeserie);
     // tab  pas un espace. il sagit de la matrice traitée en Excell  donc txt tab delimited
     while((end_line_timeserie=line_timeserie.find('	',0))!=string::npos){
       tempnum_timeserie = line_timeserie.substr(0,end_line_timeserie);
       stringstream(tempnum_timeserie) >> brol_timeserie;
       input_vector_timeserie.push_back(brol_timeserie);
       line_timeserie.erase(0,end_line_timeserie+1);
     }
     tempnum_timeserie = line_timeserie.substr(0,line_timeserie.length());
     stringstream(tempnum_timeserie) >>brol_timeserie;
     input_vector_timeserie.push_back(brol_timeserie);
  }
 name_file.close();
 }else{cout << "file can not be opened"  << endl;}


 string line_screenshot, tempnum_screenshot;
   int end_line_screenshot = 0;
   vector<double>  input_vector_screenshot;
   vector<double>::iterator  input_vector_screenshot_it;
   double brol_screenshot;
   ifstream name_file3(  SCREENSHOT_DATA.c_str(), ios::in);
 if (name_file3.is_open()){
    while(!name_file3.eof()){
       getline(name_file3, line_screenshot);
     // tab  pas un espace. il sagit de la matrice traitée en Excell  donc txt tab delimited
     while((end_line_screenshot=line_screenshot.find('	',0))!=string::npos){
       tempnum_screenshot = line_screenshot.substr(0,end_line_screenshot);
       stringstream(tempnum_screenshot) >> brol_screenshot;
       input_vector_screenshot.push_back(brol_screenshot);
       line_screenshot.erase(0,end_line_screenshot+1);
     }
     tempnum_screenshot = line_screenshot.substr(0,line_screenshot.length());
     stringstream(tempnum_screenshot) >>brol_screenshot;
     input_vector_screenshot.push_back(brol_screenshot);
  }
 name_file3.close();
 }else{cout << "file can not be opened"  << endl;}


/////MATRICE DE CONTROLE DES ETATS GAMMA////////////////////////////////////////////////////////
double** stochio_matrix_Nr_DvDx=NULL;
stochio_matrix_Nr_DvDx = new double*[dim_Nr+1];
for (int i=1;i<= dim_Nr; ++i){stochio_matrix_Nr_DvDx[i] = new double[number_vertices+1];}

alglib::real_2d_array stochio_matrix_Nr_DvDx_L_complet;
stochio_matrix_Nr_DvDx_L_complet.setlength(dim_Nr,dim_Nr);
alglib::ae_int_t info_inv;
alglib::matinvreport rep;

double** L_complet_stochio_matrix_Nr_DvDx_L_complet=NULL;
L_complet_stochio_matrix_Nr_DvDx_L_complet = new double*[number_vertices+1];
for (int i=1;i<= number_vertices; ++i){L_complet_stochio_matrix_Nr_DvDx_L_complet[i] = new double[dim_Nr+1];}

double** matrix_GAMMA=NULL;
matrix_GAMMA = new double*[number_vertices+1];
for (int i=1;i<= number_vertices; ++i){matrix_GAMMA[i] = new double[nbre_reactions+1];}

double** matrix_C=NULL;
matrix_C = new double*[nbre_reactions+1];
for (int i=1;i<= nbre_reactions; ++i){matrix_C[i] = new double[nbre_reactions+1];}

double** Bit2=NULL;
Bit2 = new double*[dim_Nr];
for (int i=0;i< dim_Nr; ++i){Bit2[i] = new double[dim_Nr];}


int number_selected_tp = floor(input_vector_screenshot.size()/3.0);
int time_point;
int label_X;

////////////////DEBUT DE LA BOUCLE SUR LES TIMESTEPS SELECTIONNES -  CF FICHIER SCREENSHOT ID///////////////////////////
for (int main_loop = 1; main_loop<=number_selected_tp;++main_loop ){

input_vector_screenshot_it =input_vector_screenshot.begin()+ (main_loop-1)*3; 
time_point = *input_vector_screenshot_it;

label_X = 0;
for (input_vector_timeserie_it = (input_vector_timeserie.begin()+(time_point - 1)*(number_vertices+1));input_vector_timeserie_it!= (input_vector_timeserie.begin()+time_point*(number_vertices+1));++input_vector_timeserie_it)
 {
   label_X++;
   if (label_X <=number_vertices ){X[label_X] =*input_vector_timeserie_it ;}
 }

Fluxes_rates_analytical(DvDx_matrix,X,k);


///1) calcul  Nr*DvDx
for (int i=1;  i<=dim_Nr;  ++i){
  for (int j=1; j<=number_vertices; ++j){
  stochio_matrix_Nr_DvDx[i][j] = 0;
   for (int k=1; k<=nbre_reactions; ++k){
     stochio_matrix_Nr_DvDx[i][j] += stochio_matrix_Nr[i][k]*DvDx_matrix[k][j];
	 
   }  
  }
}

///2)   Nr*DvDx est multipliee a droite par L_complet
for (int i=1;  i<=dim_Nr;  ++i){
  for (int j=1; j<=dim_Nr; ++j){
  stochio_matrix_Nr_DvDx_L_complet[i-1][j-1] = 0;
   for (int k=1; k<=number_vertices; ++k){
     stochio_matrix_Nr_DvDx_L_complet[i-1][j-1] += stochio_matrix_Nr_DvDx[i][k]*matrix_coef_L_complet[k][j];
	 
   }  
  }
}

///3) inversion de la matrice  stochio_matrix_Nr_DvDx_L_complet, comme elle peut presenter des singularites, on utilise SVD

alglib::rmatrixsvd(stochio_matrix_Nr_DvDx_L_complet, dim_Nr, dim_Nr, 2, 2, 2, Singular, U, V);
a_row = Singular.getcontent();

for (int i=0;  i<dim_Nr;  ++i){
  for (int j=0; j< dim_Nr; ++j){
    Bit2[i][j] = 0;
   if(a_row[i] !=0  ){ Bit2[i][j] = (1/a_row[i])*U[j][i];}else{ Bit2[i][j] = 0;}
  }
}

for (int i=0;  i<dim_Nr;  ++i){
  for (int j=0; j<dim_Nr; ++j){
  stochio_matrix_Nr_DvDx_L_complet[i][j] = 0;
   for (int k=0; k<dim_Nr; ++k){
    stochio_matrix_Nr_DvDx_L_complet[i][j] += V[k][i]*Bit[k][j];
   }  
  }
}


//4) la matrice inverse Nr*DvDx*L est multipliee a gauche par L_complet
for (int i=1;  i<=number_vertices;  ++i){
  for (int j=1; j<=dim_Nr; ++j){
  L_complet_stochio_matrix_Nr_DvDx_L_complet[i][j] = 0;
   for (int k=1; k<=dim_Nr; ++k){
     L_complet_stochio_matrix_Nr_DvDx_L_complet[i][j] += matrix_coef_L_complet[i][k]*stochio_matrix_Nr_DvDx_L_complet[k-1][j-1];
   }  
  }
}

//5) la matrice L_complet*(Nr*DvDx*L)^-1 est multipliee a droite par Nr , le tout est multiplie par -1 ce qui donne GAMMA
for (int i=1;  i<=number_vertices;  ++i){
  for (int j=1; j<=nbre_reactions; ++j){
  matrix_GAMMA[i][j] = 0;
   for (int k=1; k<=dim_Nr; ++k){
     matrix_GAMMA[i][j] += L_complet_stochio_matrix_Nr_DvDx_L_complet[i][k]*stochio_matrix_Nr[k][j];
   }  
  matrix_GAMMA[i][j] = -1.0*matrix_GAMMA[i][j];
  }
}

//6) calcul dela matrice des coefficients de flux C

 for (int i=1;  i<=nbre_reactions;  ++i){
  for (int j=1; j<=nbre_reactions; ++j){
  matrix_C[i][j] = 0;
   for (int k=1; k<=number_vertices; ++k){
      matrix_C[i][j] +=DvDx_matrix[i][k]*matrix_GAMMA[k][j];
   }  
  if (i == j ){matrix_C[i][j] += 1.0;} 
  }
}
 
 output_matrices_fluxes(matrix_GAMMA,  matrix_C, nbre_reactions, number_vertices, time_point, X_echelon_label);
 }// fin boucle sur main_loop
 
/////////////////////////////////////// FIN BOUCLE SUR LES TIMESTEPS SELECTIONNES ///////////////

///////////////////////////////////////  ELIMINATION VARIABLES ETC...///////////////////////////
for (int i=1;i<= nbre_reactions;++i){delete[] DvDx_matrix[i];}
delete []  DvDx_matrix;
for (int i=1;i<= number_vertices;++i){delete[] stochio_matrix[i];}
delete [] stochio_matrix;
for (int i=1;i<= number_vertices;++i){delete[] stochio_matrix_echelon[i];}
delete [] stochio_matrix_echelon;
for (int i=1;i<= dim_Nr;++i){delete[] stochio_matrix_Nr[i];}
delete [] stochio_matrix_Nr;
for (int i=1;i<= (number_vertices-dim_Nr);++i){delete[] stochio_matrix_N_dependant[i];}
delete [] stochio_matrix_N_dependant;
for (int i=1;i<= (number_vertices-dim_Nr);++i){delete[] matrix_coef_L[i];}
delete [] matrix_coef_L;
for (int i=1;i<= number_vertices;++i){delete[] matrix_coef_L_complet[i];}
delete [] matrix_coef_L_complet;


for (int i=1;i<= dim_Nr;++i){delete[] stochio_matrix_Nr_DvDx[i];}
delete [] stochio_matrix_Nr_DvDx;
for (int i=1;i<= number_vertices;++i){delete[]  L_complet_stochio_matrix_Nr_DvDx_L_complet[i];}
delete []  L_complet_stochio_matrix_Nr_DvDx_L_complet;
for (int i=1;i<= number_vertices;++i){delete[]  matrix_GAMMA[i];}
delete [] matrix_GAMMA;
for (int i=1;i<= nbre_reactions;++i){delete[]  matrix_C[i];}
delete [] matrix_C;
for (int i=0;i<= dim_Nr;++i){delete[]  Bit2[i];}
delete [] Bit2;
for (int i=0;i< nbre_reactions; ++i){delete [] Bit[i];}
delete [] Bit;
for (int i=0;i< nbre_reactions; ++i){delete [] Nr_pseudoinverse[i];}
delete [] Nr_pseudoinverse;

delete []  X_echelon_label;
delete [] X_echelon_label_test;
delete []  X;
delete []  k; 
delete []  base_row;

input_vector_timeserie.clear();
input_vector_screenshot.clear();

} // end if MCA

end = clock();
double cpu_duration = (end-start)/CLOCKS_PER_SEC;
output_statistics_sierraII(cpu_duration, MCA);

cout << ""<< endl;
cout <<"==================================================================================================================="<< endl;
cout << "               ECRITURE DES RESULTATS  TERMINEE, FIN DU PROGRAMME.  DUREE DE L'OPERATION: "<< cpu_duration<<" SECONDES"<< endl;
cout <<"==================================================================================================================="<< endl;
cout << ""<< endl;



 
return 0;
} // end main function

