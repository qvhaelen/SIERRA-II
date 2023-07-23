 /* fstream: stream  class to both read and write  from/to  files */
#include <fstream>
/* utilisation des chaines de caractères*/
/* contient les éléments cin et cout entre autre*/
#include <iostream>
#include <string>
/* utilisation des manipulateurs*/
#include<iomanip>
/*conversion nombre string et tring nombre (stringstream,  etc)*/
#include<sstream>
#include<ctime>
#include"output_statistics_sierraII.h"
using namespace std;
void output_statistics_sierraII(double& cpu_duration,bool& MCA ){
string OUTPUT_FILE_NAME ="output_data/";
string  filename =OUTPUT_FILE_NAME+ "OUTPUT-SIERRAII-STATISTIQUES.txt";
time_t rawtime;
struct tm *timeinfo;
char format[32];
time(&rawtime);
timeinfo = localtime(&rawtime);
ofstream output_file(filename.c_str());
if (output_file.is_open()){
  output_file<<"*************************************************************************************************"<< endl;
  output_file<< "Date de la simulation (écriture des résultats et fin du programme): "<< endl; 
  output_file<< asctime(timeinfo) <<endl;
 output_file<<"*************************************************************************************************"<< endl;
  output_file<<"Elapsed CPU time for the  completion of the job is: " <<cpu_duration << " seconds"<< " = "<<  cpu_duration/60.0 <<" minutes " << endl;
  if(MCA){output_file<< "Computation of MCA analysis: yes"<<endl;}else{output_file<< "Computation of MCA analysis:  no"<<endl;}
  output_file<<"*********************************** FIN DES STATISTIQUES *****************************************" << endl;
  output_file.close();
}else{cout << "file can not be opened"  << endl;}
return;}


























