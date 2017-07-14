//Wrapper for the biomechanical model function

#include "worm_trial.h"
#include <iostream>
#include <fstream>
#include <omp.h>
#include "string"
#include <vector>

using namespace std;

int main();
void run_trials_parallel();

int main()
{
    //Everything to set up the parallel processing section
    int id;
    double wtime;

    printf ( "\n" );
    printf ( "Worm simulations with OpenMP\n" );
    printf ( "  C/OpenMP version\n" );
    printf ( "\n" );
    printf ( "  Number of processors available = %d\n", omp_get_num_procs ( ) );
    printf ( "  Number of threads =              %d\n", omp_get_max_threads ( ) );

    wtime = omp_get_wtime ( );

    //Run the parallelized trials
    run_trials_parallel();

    wtime = omp_get_wtime ( ) - wtime;
    printf ( "\n" );
    printf ( "  Elapsed wall clock time = %f\n", wtime );

    return(0);
}

void run_trials_parallel(){
    //Runs the trials in a parallelized way
    //Uses openMP
    
# pragma omp parallel
{
    int id = omp_get_thread_num ( );
    printf ( "Process %d successfully started up.\n", id ) ;

# pragma omp for collapse(2)
    //LOTS of nested loops, for each body variable
    for(int iPerc=-5;iPerc<6;++iPerc){
        for(int iVar=0;iVar<7;++iVar){
            //Loops through percent differences and variable names
            //NOT USING THE ENVIRONMENTAL VARIABLES

            vector<char*> vNames = {(char*)"L_SEG",(char*)"K_PE",(char*)"D_PE",(char*)"AE_PE_ratio",(char*)"K_DE",(char*)"D_DE",(char*)"T_MUSCLE"};

            char* thisVName;
            thisVName = vNames.at(iVar);

            char pramsFileName[50];
            //For now, just do one instance of the suppression wave
            //  These are NOT looped over
            char DURATION[5] = "10.0";
            char SR_B_TSTART[5] = "1.6";
            char SR_B_TEND[5] = "5.1";
            char SR_B_SLOPE[5] = "2.0";
            char NMJ_B_TSTART[5] = "1.6";
            char NMJ_B_TEND[5] = "5.1";
            char NMJ_B_SLOPE[5] = "2.0";
            char NMJ_WEIGHT_B_FACTOR[5] = "0.5";
            //These will be looped over
            char vNum[5] = "1.0";

            int flag = 1;
        
            //Create the strings to be written to the file
            //  These should be percent changes: 90-110%
            sprintf(vNum,"%.2f",1.0+float(2*iPerc)/100.0); 

            //Let's give the file a good name
            sprintf(pramsFileName,"Jac_%s%s",
                thisVName,vNum);
            
            flag = write_prams_file(
                    20,
                    "FILENAME",pramsFileName,
                    "DURATION",DURATION,
                    "SR_B_TSTART",SR_B_TSTART,
                    "SR_B_TEND",SR_B_TEND,
                    "SR_B_SLOPE",SR_B_SLOPE,
                    "NMJ_B_TSTART",NMJ_B_TSTART,
                    "NMJ_B_TEND",NMJ_B_TEND,
                    "NMJ_B_SLOPE",NMJ_B_SLOPE,
                    "NMJ_WEIGHT_B_FACTOR",NMJ_WEIGHT_B_FACTOR,
                    thisVName,vNum);
            
            if(flag == 0 || flag == 2){
                worm_trial(pramsFileName);
            }
                
            }//All 7 for loops
        }
    }
}

