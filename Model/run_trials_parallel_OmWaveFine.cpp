//Wrapper for the biomechanical model function

#include "worm_trial.h"
#include <iostream>
#include <fstream>
#include <omp.h>

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
    
# pragma omp parallel default(none)
{
    int id = omp_get_thread_num ( );
    printf ( "Process %d successfully started up.\n", id ) ;

# pragma omp for collapse(2)
    //Need to have two perfectly nested loops in order to use "collapse"
    for(int iDur=0;iDur<10;++iDur){ //Cycle through wave durations
       for(int iWeight=0;iWeight<5;++iWeight){ //Cycle through NMJ weights

        char pramsFileName[50];
        char DURATION[5] = "10.0";
        char SR_B_TSTART[5] = "0.0";
        char SR_B_TEND[5] = "5.0";
        char SR_B_SLOPE[5] = "2.0";
        char NMJ_B_TSTART[5] = "0.0";
        char NMJ_B_TEND[5] = "5.0";
        char NMJ_B_SLOPE[5] = "2.0";
        char NMJ_WEIGHT_B_FACTOR[5] = "0.0";
        int flag = 1;
        float waveEnd, waveDuration, waveStart;

            for(int iStart=0;iStart<20;++iStart){ //Cycle through start times  
                waveStart = 2.0+float(iStart)/10.0; //Start times from 3.0 to 5.0
                waveDuration = 3.0+float(iDur)/2.0; //From 3.0 to 8.0
                waveEnd = waveStart + waveDuration;

                //Create the strings to be written to the file
                sprintf(SR_B_TSTART,"%.1f",waveStart);
                sprintf(NMJ_B_TSTART,"%.1f",waveStart);
                sprintf(SR_B_TEND,"%.1f",waveEnd);
                sprintf(NMJ_B_TEND,"%.1f",waveEnd);
                sprintf(DURATION,"%.1f",5.0+waveEnd);
                sprintf(NMJ_WEIGHT_B_FACTOR,"%.1f",float(iWeight*2)/10); //From 0.0 to 1.0

                //Let's give the file a good (unique) name
                sprintf(pramsFileName,"om_tst_%s_tend%s_NMJw%s_tyr",SR_B_TSTART,SR_B_TEND,NMJ_WEIGHT_B_FACTOR);
                
                flag = write_prams_file(
                        18, //IMPORTANT: this sets how many variables we're writing to the parameter file
                        "FILENAME",pramsFileName,
                        "DURATION",DURATION,
                        "SR_B_TSTART",SR_B_TSTART,
                        "SR_B_TEND",SR_B_TEND,
                        "SR_B_SLOPE",SR_B_SLOPE,
                        "NMJ_B_TSTART",NMJ_B_TSTART,
                        "NMJ_B_TEND",NMJ_B_TEND,
                        "NMJ_B_SLOPE",NMJ_B_SLOPE,
                        "NMJ_WEIGHT_B_FACTOR",NMJ_WEIGHT_B_FACTOR);
                
                if(flag == 0 || flag == 2){
                    worm_trial(pramsFileName);
                }
                
                }
            }
        }
    }
}
