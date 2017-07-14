//Wrapper for the biomechanical model function

#include "worm_trial.h"
#include <iostream>
#include <fstream>
#include <omp.h>
#include "string"

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

# pragma omp for
    for(int i=1;i<1;++i){

        char pramsFileName[50];
        //For now, just do one instance of the suppression wave
        //  These are NOT looped over
        char DURATION[5] = "10.0";
        char SR_B_TSTART[5] = "1.6";
        char SR_B_TEND[5] = "5.1";
        char SR_B_SLOPE[5] = "2.0";

        int flag = 1;

        //Let's give the file a good name
        sprintf(pramsFileName,"test_single_run");
        
        flag = write_prams_file(
                10, //This is the number of parameters to write in the file, including the name!
                "FILENAME",pramsFileName,
                "DURATION",DURATION,
                "SR_B_TSTART",SR_B_TSTART,
                "SR_B_TEND",SR_B_TEND,
                "SR_B_SLOPE",SR_B_SLOPE);
        
        if(flag == 0 || flag == 2){
            worm_trial(pramsFileName);
        }

        }
    }
}
