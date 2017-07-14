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

# pragma omp for collapse(7)
    //Loop over each body variable, including the medium, to get sensitivities
    // Main (last) loop is what we want to understand: the model with different lengths of stretch receptors
    for(int iMedium=1;iMedium<2;++iMedium){
        for(int iL_SEG=-1;iL_SEG<2;++iL_SEG){
            for(int ik_PE=-1;ik_PE<2;++ik_PE){
                for(int iD_PE=-1;iD_PE<2;++iD_PE){
                    for(int iAE_PE_ratio=-1;iAE_PE_ratio<2;++iAE_PE_ratio){
                        for(int iT_muscle=-1;iT_muscle<2;++iT_muscle){
                            for(int iSR=0;iSR<11;++iSR){

                                char pramsFileName[50];
                                char DURATION[5] = "10.0";
                                //For now, just do one instance of the suppression wave
                                //  These are NOT looped over
                                char N_SR[5] = "6";
                                char MEDIUM[5] = "1.0";
                                //These will be looped over
                                char L_SEG[5] = "1.0";
                                char k_PE[5] = "1.0";
                                char D_PE[5] = "1.0";
                                char AE_PE_ratio[5] = "1.0";
                                char k_DE[5] = "1.0";
                                char D_DE[5] = "1.0";
                                char T_muscle[5] = "1.0";

                                float MediumNum;

                                int flag = 1;
                                
                                //Create the strings to be written to the file
                                //  These should be percent changes: 98%-102%
                                sprintf(L_SEG,"%.1f",1.0+float(iL_SEG)/100.0);
                                sprintf(k_PE,"%.1f",1.0+float(ik_PE)/100.0);
                                sprintf(D_PE,"%.1f",1.0+float(iD_PE)/100.0);
                                sprintf(AE_PE_ratio,"%.1f",1.0+float(iAE_PE_ratio)/100.0);
                                sprintf(T_muscle,"%.1f",1.0+float(iT_muscle)/100.0);

                                //We want to cycle through the possible number of stretch receptors
                                sprintf(N_SR,"%d",iSR+1);
                                MediumNum = 1.0 - float(iMedium);
                                sprintf(MEDIUM,"%.1f",MediumNum);//Only 2 different media: 0 and 1

                                //Let's give the file a good name (or at least unique)
                                sprintf(pramsFileName,"env_%s_%s_%s_%s_%s_M%s_SR%s",
                                    L_SEG,k_PE,D_PE,AE_PE_ratio,T_muscle,MEDIUM,N_SR);
                                
                                //Actually write the file
                                flag = write_prams_file(
                                        18, //IMPORTANT: this sets how many variables we're writing to the parameter file
                                        "FILENAME",pramsFileName,
                                        "DURATION",DURATION,
                                        "L_SEG",L_SEG,
                                        "K_PE",k_PE,
                                        "D_PE",D_PE,
                                        "AE_PE_RATIO",AE_PE_ratio,
                                        "T_MUSCLE",T_muscle,
                                        "N_SR",N_SR,
                                        "MEDIUM",MEDIUM);
                                
                                if(flag == 0 || flag == 2){
                                    worm_trial(pramsFileName);
                                }
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

