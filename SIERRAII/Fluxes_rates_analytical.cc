#include<cmath>
#include "Fluxes_rates_analytical.h"

void Fluxes_rates_analytical(double** DvDx_matrix, double* X, double* k){
DvDx_matrix[9][2]=k[9]*X[3];
DvDx_matrix[9][3]=k[9]*X[2];
DvDx_matrix[14][5]=k[14]*X[8];
DvDx_matrix[14][8]=k[14]*X[5];
DvDx_matrix[15][5]=k[15]*X[16];
DvDx_matrix[15][16]=k[15]*X[5];
DvDx_matrix[16][9]=k[16]*X[17];
DvDx_matrix[16][17]=k[16]*X[9];
DvDx_matrix[17][11]=k[17]*X[18];
DvDx_matrix[17][18]=k[17]*X[11];
DvDx_matrix[18][1]=k[18]*X[20];
DvDx_matrix[18][20]=k[18]*X[1];
return;}


