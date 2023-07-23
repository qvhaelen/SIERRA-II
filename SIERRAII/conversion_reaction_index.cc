#include<cmath>
#include<iomanip>
#include"conversion_reaction_index.h"
void  conversion_reaction_index(double** Matrice_Adt_transition){
int* vector_conversion;
vector_conversion  = new int[37];
double** Matrice_Adt_transition_bis = NULL;
Matrice_Adt_transition_bis = new double*[37];
for (int i=1;i<=36;++i){Matrice_Adt_transition_bis[i] = new double [37];}
vector_conversion[1] = 25;
vector_conversion[2] = 26;
vector_conversion[3] = 24;
vector_conversion[4] = 27;
vector_conversion[5] = 19;
vector_conversion[6] = 31;
vector_conversion[7] = 36;
vector_conversion[8] = 8;
vector_conversion[9] = 33;
vector_conversion[10] = 30;
vector_conversion[11] = 35;
vector_conversion[12] = 20;
vector_conversion[13] = 7;
vector_conversion[14] = 34;
vector_conversion[15] = 21;
vector_conversion[16] = 32;
vector_conversion[17] = 6;
vector_conversion[18] = 10;
vector_conversion[19] = 11;
vector_conversion[20] = 5;
vector_conversion[21] = 3;
vector_conversion[22] = 4;
vector_conversion[23] = 12;
vector_conversion[24] = 13;
vector_conversion[25] = 1;
vector_conversion[26] = 2;
vector_conversion[27] = 23;
vector_conversion[28] = 29;
vector_conversion[29] = 28;
vector_conversion[30] = 22;
vector_conversion[31] = 17;
vector_conversion[32] = 16;
vector_conversion[33] = 14;
vector_conversion[34] = 15;
vector_conversion[35] = 9;
vector_conversion[36] = 18;
for (int i = 1; i<=36;++i){
   for (int j = 1; j<=36;++j){
        Matrice_Adt_transition_bis[vector_conversion[i]][j] = Matrice_Adt_transition[i][j];
   }
}

for (int i = 1; i<=36;++i){
   for (int j = 1; j<=36;++j){
        Matrice_Adt_transition[i][vector_conversion[j]] = Matrice_Adt_transition_bis[i][j];
   }
}

 for (int i=1;i<=36;++i){delete [] Matrice_Adt_transition_bis[i];}
 delete [] Matrice_Adt_transition_bis;
 delete [] vector_conversion;
return;}


