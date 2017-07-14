/* 
This software is freely available for use in research, teaching and other
non-commercial purposes.  Users have the right to modify, alter, improve,
or enhance the software without limitation, under the condition that you
do not remove or alter the copyright information in this section.

If you publish results which were obtained using the software, or its
source code, please cite the software with a reference to the associated 
publication:

J.H. Boyle, S. Berri and N. Cohen (2012), Gait modulation in C. elegans:
an integrated neuromechanical model, Front. Comput. Neurosci, 6:10 
doi: 10.3389/fncom.2012.00010

This code has been significantly modified by Charles Fieseler (University 
of Washington, Physics department), and if this version of the code is used
please also cite:

.............

Licence information for Sundials IDA can be found in "sundials-2.6.0/LICENCE"
*/

// Required includes
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <stdarg.h>     /* va_list, va_start, va_arg, va_end */
#include <omp.h>


// My own additions (to get it to work with sundials 2.6.0)
#include "worm_trial.h"
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

using namespace std;

// Simulation paramaters
#define OBJECTS 0				//set number of Objects (>= 0)
#define LAYOUT 0				//change between 0 (square) 1 (hex) 2 (random)
#define N_UNITS 12				//Number of neural segments

// Simulator constants
#define NSEG 48					//Number of physical segments (4 per neural segment)
#define NBAR NSEG+1				//Number of bars, which are between each physical segment
#define NSEG_MINUS_1 NSEG-1
#define NEQ   3*(NBAR)			//Number of equations to solve
#define DELTAT 0.001			//Timestep
#define HALFPI M_PI/2.0

struct UserEnviron {

// Simulation paramaters
realtype DURATION = 10.0;
realtype MEDIUM = 1.0;

//initializer boolean and counts for downsampling the output data files
bool neurons_initialized = false;
int tCount = 0;
int tCount_meta = 0;

// USES DEFINE STATEMENT ABOVE
int State_A[N_UNITS][2];
int State_B[N_UNITS][2]; 

//Number of stretch receptors
int N_SR = 6;

//	Used in: update_neurons
float SRvelocity = 1.0;
float Hyst = 0.5;
float Gap_coupling = 0.0;

/*****   WAVE OF DEPRESSION VARIABLES   ******/
//	 on the stretch receptors, both A and B sides
bool isInh_I_SR = true; //If the stretch receptor current can be inhibitory

// Stretch receptor suppression variables for the B-class (forward motion) neurons
float SR_B_vel = 0.25; //Speed of the wave, in weird units
float SR_B_tStart = 0.0; //Start of the wave
float SR_B_tEnd = 0.0; //End of the wave
float SR_B_slope = 0.0; //If the slope is 0.0, then the doubleSigmoid will exit early and return 0.0 (i.e. no wave)
float SR_weight_factor_B = 2.0; //If the stretch receptors can't be inhibitory, this should be increased

//Same as above but for A-class (backwards motion)
//	Note that backwards omega turns are never observed in normal worms
float SR_A_vel = 0.25;
float SR_A_tStart = 0.0;
float SR_A_tEnd = 0.0;
float SR_A_slope = 0.0;
float SR_weight_factor_A = 2.0;

//Wave of excitation or depression on the I_AVB current
float I_B_sign = 1.0;
float I_B_tStart = 0.0;
float I_B_tEnd = 0.0;
float I_B_slope = 0.0;

//Wave of excitation or depression on the I_bias term
//	This term contributes to the asymmetry of the model
float I_bias_D_inh = 0.0;
float I_bias_sign = 1.0;
float I_bias_tStart = 0.0;
float I_bias_tEnd = 0.0;
float I_bias_slope = 0.0;

//Wave of excitation on the NMJs, both A and B
float NMJ_B_vel = 0.0625; //We want this to be 1/4 of the SR_vel by default because there are 4x the segments here
float NMJ_B_tStart = 0.0;
float NMJ_B_tEnd = 0.0;
float NMJ_B_slope = 0.0;
float NMJ_weight_B_factor = 1.0;

float NMJ_A_vel = 0.0625;
float NMJ_A_tStart = 100.0;
float NMJ_A_tEnd = 111.0;
float NMJ_A_slope = 0.0;
float NMJ_weight_A_factor = 1.0;

//tyramine (can prevent VD activation)
//	Note that this hasn't been fully tested and doesn't really give good/realistic results
float tyr_sign = 1.0;
float tyr_tStart = 100.0;
float tyr_tEnd = 111.0;
float tyr_slope = 0.0;


/*****   BODY CONSTANTS   ******/
//	Note that if these are changed in the parameter file, they are changed by that percent
//		Note also that some of them, like the body radius aren't directly set here
// General body constants
realtype D = 80e-6; //Max body diameter
realtype R[NBAR]; //Body radius for each bar
realtype L_seg = 1e-3/NSEG; //Length of a segment

// Horizontal element constants
realtype k_PE = (NSEG/24.0)*RCONST(10.0e-3);  	//PASSIVE ELEMENT (PE)
realtype D_PE = RCONST(0.025)*k_PE;				//DAMPING (D_**)
realtype AE_PE_ratio = RCONST(20.0);    		//Ratio of active element to passive element
realtype k_AE = AE_PE_ratio*k_PE;				//ACTIVE ELEMENT (AE)
realtype D_AE = RCONST(5.0)*AE_PE_ratio*D_PE;	

// Diagonal element constants
realtype k_DE = RCONST(350.0)*k_PE;				//DIAGONAL ELEMENT (DE)
realtype D_DE = RCONST(0.01)*k_DE;	

// Length constant holders
realtype L0_P[NSEG]; 							//Default spring length
realtype L_min[NSEG];							//Minimum length
realtype L0_P_minus_L_min[NSEG];	
realtype L0_D[NSEG];							//Default spring length

// Stretch receptor constant holders
realtype SR_shape_compensation[NSEG]; 			//Stretch receptor sensitivity is increased to compensate for the
												// decreasing number of sensory inputs towards the tail

// Muscle time constant
realtype T_muscle = RCONST(0.1);	


/*****   ENVIRONMENTAL CONSTANTS   ******/
//	Note that if these are changed in the parameter file, they are changed by that percent
// UserEnviron constants
realtype CL_water = RCONST(3.3e-6/(2.0*NBAR));
realtype CN_water = RCONST(5.2e-6/(2.0*NBAR));
realtype CL_agar = RCONST(3.2e-3/(2.0*NBAR));
realtype CN_agar = RCONST(128e-3/(2.0*NBAR));

// UserEnviron variables
realtype K_agar = CN_agar/CL_agar;
realtype CN[NBAR];
realtype CL[NBAR];
realtype ContactForce;
static const int N_objects = OBJECTS;
int root_N = round(sqrt(N_objects));
float Objects[N_objects][3];
realtype k_Object = k_PE*5;


/*****   COMMUNICATION VARIABLES   ******/
// Communication variables  
realtype L_SR[NSEG][2]; 						//Length of the stretch receptors
realtype I_SR[NSEG][2]; 						//Current from the stretch receptors

// Neuron and muscle state variables
realtype V_muscle[NSEG][2];
realtype V_neuron[NSEG][2];

// Declare output file name (these are default)
char * outfileName = (char *)"simdata_test.csv";
char * outfileNameMeta = (char *)"metadata_test.csv";

};

//Overload the output operator so we can see some actual values
ostream& operator<<(ostream& os, const UserEnviron env){
	return os << "DURATION " << env.DURATION << endl
			  << "MEDUIM " << env.MEDIUM << endl
			  << "SR_B_TSTART " << env.SR_B_tStart << endl
			  << "SR_B_TEND " << env.SR_B_tEnd << endl
			  << "SR_B_SLOPE " << env.SR_B_slope << endl
			  << "NMJ_B_TSTART " << env.NMJ_B_tStart << endl
			  << "NMJ_B_TEND " << env.NMJ_B_tEnd << endl
			  << "NMJ_B_SLOPE " << env.NMJ_B_slope << endl
			  << "NMJ_WEIGHT_B_FACTOR " << env.NMJ_weight_B_factor << endl;
}



// Prototypes of functions called by IDA (Copied from Sundials examples)
void update_external(realtype timenow, UserEnviron* envi);
	//Updates body model, i.e. calculates the rod-spring forces and viscosity interactions
int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector resval, void *rdata);
	//Residual function of the form F-ma=0
	//To get real motion, the IDA solver tries to minimize the return of this function
void update_neurons(realtype timenow, ofstream& outfile, UserEnviron* envi);
	//Updates bistable neuron model
void update_muscles(realtype timenow, UserEnviron* envi);
	//Updates muscles using a linear sigmoidal approximation
void update_SR(realtype timenow, UserEnviron* envi);
	//Updates the stretch receptors using the body segment lengths

int write_prams_file(int n, ...);
	//Writes parameter files with "n" user-passed variables

// Prototypes of private functions (Copied from Sundials examples)
static void SaveOutput(void *mem, realtype t, N_Vector y, ofstream& outfile, UserEnviron* envi);
static void PrintFinalStats(void *mem);
static int check_flag(void *flagvalue, char *funcname, int opt);
double randn(double mu, double sigma);

float doubleSigmoid(float timenow, float tstart, float tend, float slope);
	//Implements the double sigmoidal wave function

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */


int worm_trial(const char* pramsFileName)
{
	//This function initializes all the variables, sets up the IDA solver, and loops until the simulation is finished
	//The only input is the file name to write parameters to

  // IDA variables (Copied from Sundials examples)
  void *mem;
  N_Vector yy, yp, avtol;
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, tout, tret;
  int iout, retval, retvalr;  
  FILE *failFile;
  failFile=fopen("fail_log.txt","a");

  //Should be initialized to the correct values
  UserEnviron envObj;

  ofstream outfile;
  ofstream outfile_meta;

  yy = yp = avtol = NULL;
  yval = ypval = atval = NULL;

  // Allocate N-vectors (Copied from Sundials examples)
  yy = N_VNew_Serial(NEQ);
  if(check_flag((void *)yy, (char *)"N_VNew_Serial", 0)) return(1);
  yp = N_VNew_Serial(NEQ);
  if(check_flag((void *)yp, (char *)"N_VNew_Serial", 0)) return(1);
  avtol = N_VNew_Serial(NEQ);
  if(check_flag((void *)avtol, (char *)"N_VNew_Serial", 0)) return(1);

  //Import parameters from text file with two columns: variable name (hardcoded in here);	value (when written, will have been a string)
  if(pramsFileName[0]){

	//Create the correct full parameter file name
	char* extStr = (char *)".csv";
	char* preStr = (char *)"prams_";
	char * pramsFullFileName = (char *) malloc(1 + strlen(preStr) + strlen(pramsFileName)+ strlen(extStr) );
	strcpy(pramsFullFileName, preStr);
	strcat(pramsFullFileName, pramsFileName);
	strcat(pramsFullFileName, extStr);
	ifstream pFile;
	pFile.open(pramsFullFileName);

	//Create the metadata file name, which has the SAME inner name as the prams file
	preStr = (char *)"mdata_";
	extStr = (char *)".csv";
	envObj.outfileNameMeta = (char *) malloc(1+ strlen(preStr) + strlen(pramsFileName)+ strlen(extStr) );
	strcpy(envObj.outfileNameMeta, preStr);
	strcat(envObj.outfileNameMeta, pramsFileName);
	strcat(envObj.outfileNameMeta, extStr);
	
	//Variables for line parsing
	string thisVarName, restOfLine;

	//Variable for changing an environmental constant by a percentage
	float percentChange = 0.0;

	// Read in variables from parameter file, to be saved in the struct UserEnviron
	while(pFile >> thisVarName)
	{
		cout << "Read variable " << thisVarName;
		if(!thisVarName.compare("DURATION")){
			pFile >> envObj.DURATION;}
		else if(!thisVarName.compare("FILENAME")){
			string tmp;
			pFile >> tmp;
			const char* tmpOutfileName = tmp.c_str();
			const char* extStr = ".csv";
			const char* preStr = "simdata_";
			envObj.outfileName = (char *) malloc(1 + strlen(preStr) + strlen(tmpOutfileName)+ strlen(extStr) );
			strcpy(envObj.outfileName, preStr);
			strcat(envObj.outfileName, tmpOutfileName);
			strcat(envObj.outfileName, extStr);
		}
		else if(!thisVarName.compare("MEDIUM")){
			pFile >> envObj.MEDIUM;}
		/*****   DEPRESSION WAVE VARIABLES   ******/
		else if(!thisVarName.compare("SRVELOCITY")){
			pFile >> envObj.SRvelocity;}
		else if(!thisVarName.compare("HYST")){
			pFile >> envObj.Hyst;}
		else if(!thisVarName.compare("GAP_COUPLING")){
			pFile >> envObj.Gap_coupling;}
		else if(!thisVarName.compare("ISINH_I_SR")){
			pFile >> envObj.isInh_I_SR;}
		else if(!thisVarName.compare("SR_B_VEL")){
			pFile >> envObj.SR_B_vel;}
		else if(!thisVarName.compare("SR_B_TSTART")){
			pFile >> envObj.SR_B_tStart;}
		else if(!thisVarName.compare("SR_B_TEND")){
			pFile >> envObj.SR_B_tEnd;}
		else if(!thisVarName.compare("SR_B_SLOPE")){
			pFile >> envObj.SR_B_slope;}
		else if(!thisVarName.compare("SR_WEIGHT_FACTOR_B")){
			pFile >> envObj.SR_weight_factor_B;}
		else if(!thisVarName.compare("SR_A_VEL")){
			pFile >> envObj.SR_A_vel;}
		else if(!thisVarName.compare("SR_A_TSTART")){
			pFile >> envObj.SR_A_tStart;}
		else if(!thisVarName.compare("SR_A_TEND")){
			pFile >> envObj.SR_A_tEnd;}
		else if(!thisVarName.compare("SR_A_SLOPE")){
			pFile >> envObj.SR_A_slope;}
		else if(!thisVarName.compare("SR_WEIGHT_FACTOR_A")){
			pFile >> envObj.SR_weight_factor_A;}
		else if(!thisVarName.compare("I_B_SIGN")){
			pFile >> envObj.I_B_sign;}
		else if(!thisVarName.compare("I_B_TSTART")){
			pFile >> envObj.I_B_tStart;}
		else if(!thisVarName.compare("I_B_TEND")){
			pFile >> envObj.I_B_tEnd;}
		else if(!thisVarName.compare("I_B_SLOPE")){
			pFile >> envObj.I_B_slope;}
		else if(!thisVarName.compare("I_BIAS_D_INH")){
			pFile >> envObj.I_bias_D_inh;}
		else if(!thisVarName.compare("I_BIAS_SIGN")){
			pFile >> envObj.I_bias_sign;}
		else if(!thisVarName.compare("I_BIAS_TSTART")){
			pFile >> envObj.I_bias_tStart;}
		else if(!thisVarName.compare("I_BIAS_TEND")){
			pFile >> envObj.I_bias_tEnd;}
		else if(!thisVarName.compare("I_BIAS_SLOPE")){
			pFile >> envObj.I_bias_slope;}
		else if(!thisVarName.compare("NMJ_B_VEL")){
			pFile >> envObj.NMJ_B_vel;}
		else if(!thisVarName.compare("NMJ_B_TSTART")){
			pFile >> envObj.NMJ_B_tStart;}
		else if(!thisVarName.compare("NMJ_B_TEND")){
			pFile >> envObj.NMJ_B_tEnd;}
		else if(!thisVarName.compare("NMJ_B_SLOPE")){
			pFile >> envObj.NMJ_B_slope;}
		else if(!thisVarName.compare("NMJ_WEIGHT_B_FACTOR")){
			pFile >> envObj.NMJ_weight_B_factor;}
		else if(!thisVarName.compare("NMJ_A_VEL")){
			pFile >> envObj.NMJ_A_vel;}
		else if(!thisVarName.compare("NMJ_A_TSTART")){
			pFile >> envObj.NMJ_A_tStart;}
		else if(!thisVarName.compare("NMJ_A_TEND")){
			pFile >> envObj.NMJ_A_tEnd;}
		else if(!thisVarName.compare("NMJ_A_SLOPE")){
			pFile >> envObj.NMJ_A_slope;}
		else if(!thisVarName.compare("NMJ_WEIGHT_A_FACTOR")){
			pFile >> envObj.NMJ_weight_A_factor;}
		else if(!thisVarName.compare("N_SR")){
			pFile >> envObj.N_SR;}
		else if(!thisVarName.compare("TYR_SIGN")){
			pFile >> envObj.tyr_sign;}
		else if(!thisVarName.compare("TYR_TSTART")){
			pFile >> envObj.tyr_tStart;}
		else if(!thisVarName.compare("TYR_TEND")){
			pFile >> envObj.tyr_tEnd;}
		else if(!thisVarName.compare("TYR_SLOPE")){
			pFile >> envObj.tyr_slope;}
		/*****   BODY CONSTANTS   ******/
		//	All of these will be changed by percentage, not just set by the user
		else if(!thisVarName.compare("D")){
			pFile >> percentChange;
			envObj.D *= percentChange;}
		else if(!thisVarName.compare("L_SEG")){
			pFile >> percentChange;
			envObj.L_seg *= percentChange;}
		else if(!thisVarName.compare("K_PE")){
			pFile >> percentChange;
			envObj.k_PE *= percentChange;}
		else if(!thisVarName.compare("D_PE")){
			pFile >> percentChange;
			envObj.D_PE *= percentChange;}
		else if(!thisVarName.compare("AE_PE_RATIO")){
			pFile >> percentChange;
			envObj.AE_PE_ratio *= percentChange;}
		else if(!thisVarName.compare("K_DE")){
			pFile >> percentChange;
			envObj.k_DE *= percentChange;}
		else if(!thisVarName.compare("D_DE")){
			pFile >> percentChange;
			envObj.D_DE *= percentChange;}
		else if(!thisVarName.compare("T_MUSCLE")){
			pFile >> percentChange;
			envObj.T_muscle *= percentChange;}
		/*****   ENVIRONMENTAL CONSTANTS   ******/
		//	All of these will be changed by percentage, not just set by the user
		else if(!thisVarName.compare("CL_WATER")){
			pFile >> percentChange;
			envObj.CL_water *= percentChange;}
		else if(!thisVarName.compare("CN_WATER")){
			pFile >> percentChange;
			envObj.CN_water *= percentChange;}
		else if(!thisVarName.compare("CL_AGAR")){
			pFile >> percentChange;
			envObj.CL_agar *= percentChange;}
		else if(!thisVarName.compare("CN_AGAR")){
			pFile >> percentChange;
			envObj.CN_agar *= percentChange;}
		else{
			cout << "... Variable " << thisVarName << " was unused.";
			getline(pFile,restOfLine);
		}
		cout << endl;
	}
	pFile.close();
  }
  else{
	  cout << "No input file found; using default values." << endl;
  }

  // Create and initialize  y, y', and absolute tolerance vectors (Copied from Sundials examples)
  yval  = NV_DATA_S(yy);
  ypval = NV_DATA_S(yp);  
  rtol = (envObj.MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-12);
  atval = NV_DATA_S(avtol);

  
  for(int i = 0; i < NBAR; ++i){
	// Initialize body in straight line
	yval[i*3] = i*(envObj.L_seg);
	yval[i*3+1] = RCONST(0.0);
	yval[i*3+2] = M_PI/RCONST(2.0);

	// Initialize derivative values (Copied from Sundials examples)
 	ypval[i*3] = RCONST(0.0);
  	ypval[i*3+1] = RCONST(0.0);
  	ypval[i*3+2] = RCONST(0.0);

	// Set absolute tolerances for solver (Copied from Sundials examples)
	// Tolerance must be set lower when simulating in water, due to lower drag coefficients
	atval[i*3] = (envObj.MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-9);
	atval[i*3+1] = (envObj.MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-9);
	atval[i*3+2] = (envObj.MEDIUM < 0.015 ? 0.1 : 1)*RCONST(1.0e-5);
  }
  
  // Initialize model variables
  for(int i = 0; i < NSEG; ++i){
	envObj.V_muscle[i][0] = 0.0;
	envObj.V_muscle[i][1] = 0.0;
	envObj.V_neuron[i][0] = 0.0;
	envObj.V_neuron[i][1] = 0.0;
	envObj.I_SR[i][0] = 0.0;
	envObj.I_SR[i][1] = 0.0;
  }  
  
  // Set local body radius values based on elliptical approximation
  for(int i = 0; i < NBAR; ++i){	
	envObj.R[i] = envObj.D/2.0*fabs(sin(acos((i-NSEG/2.0)/(NSEG/2.0 + 0.2))));		
  } 

  // Set stretch receptor weightings that compensate for the elliptical shape,
  // giving approximately the same SR response to segment mending angle
  for(int i = 0; i < NSEG; ++i){
	envObj.SR_shape_compensation[i] = envObj.D/(envObj.R[i] + envObj.R[i+1]);  
  } 
 
  // Set muscle constants (rest length, minimum length etc) accounting 
  // for length differences due to the elliptical body shape
  for(int i = 0; i < NSEG; ++i){
	float scale = 0.65*((envObj.R[i] + envObj.R[i+1])/envObj.D);
	envObj.L0_P[i] = sqrt(pow(envObj.L_seg,2) + pow((envObj.R[i] - envObj.R[i+1]),2));
	envObj.L_min[i] = (1.0-scale)*(envObj.L0_P)[i];
	envObj.L0_P_minus_L_min[i] = envObj.L0_P[i] - envObj.L_min[i];
	envObj.L0_D[i] = sqrt(pow(envObj.L_seg,2) + pow((envObj.R[i] + envObj.R[i+1]),2));
  }

  // Set drag constants according to env->MEDIUM
  for(int i = 0; i < NBAR; ++i){
  	envObj.CL[i] = (envObj.CL_agar - envObj.CL_water)*(envObj.MEDIUM) + envObj.CL_water;
  	envObj.CN[i] = (envObj.CN_agar - envObj.CN_water)*(envObj.MEDIUM) + envObj.CN_water;
  }
		

  // Place envObj.Objects in UserEnviron (if using them)  
  if(envObj.N_objects > 0){
	if(LAYOUT == 0){
  		//Square post array ala Park et al.
		//Adjust these parameters to modify configuration
  		float post_radius = 0.1e-3;
  		float post_spacing = 0.4e-3;
  		ofstream object_outfile;
  		object_outfile.open("Objects.csv");
  		for(int i = 0; i < envObj.root_N; ++i){
			for(int j = 0; j < envObj.root_N; ++j){
				envObj.Objects[envObj.root_N*i + j][0] = post_spacing*j - 0.75*(envObj.root_N-1)*post_spacing;
				envObj.Objects[envObj.root_N*i + j][1] = post_spacing*i - 0.25*(envObj.root_N-1)*post_spacing;
				envObj.Objects[envObj.root_N*i + j][2] = post_radius;
				object_outfile << envObj.Objects[envObj.root_N*i + j][0] << ", " << envObj.Objects[envObj.root_N*i + j][1] << ", " << envObj.Objects[envObj.root_N*i + j][2] << "\n";
			}
  		}
  		object_outfile.close();	
	}
	else if(LAYOUT == 1){
		//Hexagonal post array ala Lockery et al.
		//Adjust these parameters to modify configuration
  		float post_radius = 0.1e-3;
  		float min_gap = 0.2e-3;
  		float horizontal_post_spacing = min_gap + 2*post_radius;
  		float vertical_post_spacing = sqrt(pow(horizontal_post_spacing,2) - pow(0.5*horizontal_post_spacing,2));
  		ofstream object_outfile;
  		object_outfile.open("Objects.csv");
  		for(int i = 0; i < envObj.root_N; ++i){
			for(int j = 0; j < envObj.root_N; ++j){
				envObj.Objects[envObj.root_N*i + j][0] = horizontal_post_spacing*j - 0.85*(envObj.root_N-1)*horizontal_post_spacing + 0.5*horizontal_post_spacing*(i%2 == 0);
				envObj.Objects[envObj.root_N*i + j][1] = vertical_post_spacing*i - 0.5*(envObj.root_N-1)*vertical_post_spacing;
				envObj.Objects[envObj.root_N*i + j][2] = post_radius;
				object_outfile << envObj.Objects[envObj.root_N*i + j][0] << ", " << envObj.Objects[envObj.root_N*i + j][1] << ", " << envObj.Objects[envObj.root_N*i + j][2] << "\n";
			}
  		}
  		object_outfile.close();
	}
	else if(LAYOUT == 2){
  		//Create random LAYOUT of envObj.Objects
		//Adjust these parameters to modify configuration
  		float X_lim = 1.5e-3;
  		float Y_lim = 1.5e-3;
  		float min_gap = 0.08e-3;
  		float R_min = 0.04e-3;
  		float R_max = 0.25e-3;
  		bool fits;
  		ofstream object_outfile;
  		object_outfile.open("Objects.csv");
  		for(int i = 0; i < envObj.N_objects; ++i){
			do{
				fits = true;
				envObj.Objects[i][0] = -1.5*X_lim + (2*X_lim)*(rand()/(RAND_MAX*1.0));
				envObj.Objects[i][1] = -Y_lim + (2*Y_lim)*(rand()/(RAND_MAX*1.0));
				envObj.Objects[i][2] = R_min + (R_max - R_min)*(rand()/(RAND_MAX*1.0));
	
				for (int j = 0; j < i; ++j){
					float dist = sqrt(pow((envObj.Objects[i][0] - envObj.Objects[j][0]),2) + pow((envObj.Objects[i][1] - envObj.Objects[j][1]),2));
					if(dist < (envObj.Objects[i][2] + envObj.Objects[j][2] + min_gap)){
						fits = false;
					}
				}
				for(int j = 0; j < NBAR; ++j){
					float dist = sqrt(pow((envObj.Objects[i][0] - yval[j*3]),2) + pow((envObj.Objects[i][1] - yval[j*3+1]),2));
					if(dist < (envObj.Objects[i][2] + 0.05e-3)){
						fits = false;
					}
				}
			}
			while(fits == false);
			object_outfile << envObj.Objects[i][0] << ", " << envObj.Objects[i][1] << ", " << envObj.Objects[i][2] << "\n";
  		}
  		object_outfile.close();
	} 
  	cerr << "\nObjects set up";
  }

  // Integration start time 
  t0 = RCONST(0.0);
  
  // Call IDACreate and IDAMalloc to initialize IDA memory (Copied from Sundials examples)
  mem = IDACreate();
  if(check_flag((void *)mem, (char *)"IDACreate", 0)) return(1);

  //Pass the user data (environmental variables and such)
  retval = IDASetUserData(mem, &envObj);
  if(check_flag(&retval, (char *)"IDASetUserData", 1)) return(1);

// IT LOOKS LIKE IDAMalloc IS AN OLD FUNCTION: SHOULD BE REPLACED WITH IDAInit

//  retval = IDAMalloc(mem, resrob, t0, yy, yp, IDA_SV, rtol, avtol);	
    retval = IDAInit(mem, resrob, t0, yy, yp);	
  if(check_flag(&retval, (char *)"IDAMalloc", 1)) return(1);

//Set the tolerances; note that this is a separate function in 2.6.0
    retval = IDASStolerances(mem, rtol, rtol);
  if(check_flag(&retval, (char *)"IDASStolerances", 1)) return(1);
  
  // Free avtol (Copied from Sundials examples) 
  N_VDestroy_Serial(avtol);

  // Call IDADense and set up the linear solver (Copied from Sundials examples) 
  retval = IDADense(mem, NEQ);
  if(check_flag(&retval, (char *)"IDADense", 1)) return(1);

  // Open output file  
  outfile.open(envObj.outfileName);
  outfile_meta.open(envObj.outfileNameMeta); 
  
  // Integrator inputs
  iout = 0; 
  tout = DELTAT;

  // Save current state
  SaveOutput(mem,0,yy,outfile,&envObj);


  // Loop through model for specified amount of time
  while(1) {
   
    	// Call residual function (Copied from Sundials examples)
	// to update physical model (Sundials takes multiple steps)
    	retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);    	
    
    	// Call stretch receptor update function
    	update_SR(tout, &envObj);

    	// Call neural model update function
    	update_neurons(tout, outfile_meta, &envObj);

    	//Call muscle model update function
    	update_muscles(tout, &envObj);

    	// Save current state
    	SaveOutput(mem,tret,yy,outfile,&envObj);

	// Check integration went ok (Copied from Sundials examples) 
    	if(check_flag(&retval, (char *)"IDASolve", 1)){fprintf(failFile,"Simulation %s failed.\n",pramsFileName); return(1);}
    
	// Prepare to go to next step
    	if (retval == IDA_SUCCESS) {
      		iout++;
      		tout += DELTAT;
    	}

	// End once enough simulation time has passed
    	if (tout > envObj.DURATION) break;
  }

  // (Copied from Sundials examples) 
  //PrintFinalStats(mem);

  // Free memory (Copied from Sundials examples) 
  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  
  // Close output file
  outfile.close();
  outfile_meta.close();
  fclose(failFile);

  return(0);
  
}

/*
 *--------------------------------------------------------------------
 * Model Functions
 *--------------------------------------------------------------------
 */
  // Neural circuit function
  void update_neurons(realtype timenow, ofstream& outfile_meta, UserEnviron* env){
	  //Updates the neurons after each deltaT of ODE solver
	  // Inputs:
	  // 	timenow = current time
	  // 	outfile_meta = filestream for exporting metadata; default is to export wave of suppression-related data
	  // 	UserEnviron = struct of environmental and communication variables

   	// Neural paramaters
   	const int N_units = N_UNITS;		// Number of neural units
   	int N_seg_per_unit = (int)(NSEG/N_units);

	//TIME DEPENDENT CURRENT
	float SRpos = (env->SRvelocity)*timenow;

	//Different currents for each segment, potentially
	float I_AVB[N_units], I_AVA[N_units];

	//Current from the "command neurons"
	// 	Note that 0.675 is the default in the original simulation
	for(int i = 0; i < N_units; ++i){
		I_AVB[i] = 0.675*(1.0-env->I_B_sign*doubleSigmoid(timenow,env->I_B_tStart,env->I_B_tEnd,env->I_B_slope));
		I_AVA[i] = 0.675-I_AVB[i];
	}

   	// Set up neuromuscular junctions for the A-class (backwards motion) neurons   
	// Breaks anterior-posterior symmetry
	   //Time dependence
   	float NMJ_weight_A[NSEG];
	for(int i = 0; i < NSEG; ++i){
		NMJ_weight_A[i] =  0.7*(1.0 - (NSEG - i - 1) * 0.6/NSEG) *(1.0+(env->NMJ_weight_A_factor)*doubleSigmoid(timenow,i*(env->NMJ_A_vel)+env->NMJ_A_tStart,i*(env->NMJ_A_vel)+env->NMJ_A_tEnd,env->NMJ_A_slope));
	}
   	NMJ_weight_A[NSEG-1] /= 1.5;				// Helps to prevent excessive bending of TAIL

	   //ORIGINAL, i.e. B-class (forward motion) neurons
	   //Time dependence
	float NMJ_weight_B[NSEG];
	for(int i = 0; i < NSEG; ++i){
		NMJ_weight_B[i] =  0.7*(1.0 - i * 0.6/NSEG) *(1.0+(env->NMJ_weight_B_factor)*doubleSigmoid(timenow,i*(env->NMJ_B_vel)+env->NMJ_B_tStart,i*(env->NMJ_B_vel)+env->NMJ_B_tEnd,env->NMJ_B_slope));
	}
   	NMJ_weight_B[0] /= 1.5;						// Helps to prevent excessive bending of head 

	// If this is the first time update_neurons is called, initialize with all neurons on one side ON
   	if(!(env->neurons_initialized)){
		for(int i = 0; i < N_units; ++i){
			env->State_A[i][0] = 1;
			env->State_A[i][1] = 0;
			env->State_B[i][0] = 1;
			env->State_B[i][1] = 0;
		}
		cout << "Neuron states initialized." << endl;
		env->neurons_initialized = true;

		outfile_meta << "Time " << "Stretch receptor weights" << endl;
   	}   

   	// Stretch receptor variables   	
   	float I_SR_D_A[N_units];
   	float I_SR_V_A[N_units];
   	float SR_weight_A[N_units];
	float I_SR_D_B[N_units];
   	float I_SR_V_B[N_units];
   	float SR_weight_B[N_units];

   	// SR_weight is a global weighting for each unit, used to compensate for 
	   //curvature gradient induced by the NMJ gradient above
	   //Time dependence
   	for(int i = 0; i < N_units; ++i){	
		SR_weight_A[i] = (0.65* (1.0-doubleSigmoid(SRpos,i*(env->SR_A_vel)+env->SR_A_tStart,i*(env->SR_A_vel)+env->SR_A_tEnd,env->SR_A_slope)) )*(0.4 + 0.08*(N_units-i-1))*(N_units/12.0)*(2.0/N_seg_per_unit);
		SR_weight_B[i] = (0.65* (1.0-doubleSigmoid(SRpos,i*(env->SR_B_vel)+env->SR_B_tStart,i*(env->SR_B_vel)+env->SR_B_tEnd,env->SR_B_slope)) )*(0.4 + 0.08*i)*(N_units/12.0)*(2.0/N_seg_per_unit);		
   	} 
  
  	// Reverse the direction of proprioception
	// Breaks anterior-posterior symmetry
   	// Add up stretch receptor contributions from all body segments in receptive field for each neural unit 
	   //NOTE: the mechanics of the stretch receptors don't differ for A- and B- class neurons
   	for(int i = env->N_SR; i < N_units; ++i){
		I_SR_D_A[i] = 0.0;
		I_SR_V_A[i] = 0.0;
		for(int j = 0; j < env->N_SR; ++j){
			I_SR_D_A[i] += env->I_SR[(i-j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*(env->I_SR)[(i-j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*(env->I_SR)[(i-j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i-j)*N_seg_per_unit + 3][0];	
			I_SR_V_A[i] += env->I_SR[(i-j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*(env->I_SR)[(i-j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*(env->I_SR)[(i-j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i-j)*N_seg_per_unit + 3][1];				
		}		
   	}

	for(int i = 0; i <= N_units - env->N_SR; ++i){
		I_SR_D_B[i] = 0.0;
		I_SR_V_B[i] = 0.0;
		for(int j = 0; j < env->N_SR; ++j){
			I_SR_D_B[i] += env->I_SR[(i+j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*(env->I_SR)[(i+j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*(env->I_SR)[(i+j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i+j)*N_seg_per_unit + 3][0];	
			I_SR_V_B[i] += env->I_SR[(i+j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*(env->I_SR)[(i+j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*(env->I_SR)[(i+j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i+j)*N_seg_per_unit + 3][1];				
		}		
   	}
	
	//  For units near the head, fewer segments contribute (because the body ends)
   	int tmp_N_SR = env->N_SR;
   	for(int i = env->N_SR-1; i >= 0; --i){			
		tmp_N_SR --;
		I_SR_D_A[i] = 0.0;
		I_SR_V_A[i] = 0.0;
		for(int j = 0; j < tmp_N_SR; ++j){
			I_SR_D_A[i] += env->I_SR[(i-j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*(env->I_SR)[(i-j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*(env->I_SR)[(i-j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i-j)*N_seg_per_unit + 3][0];	
			I_SR_V_A[i] += env->I_SR[(i-j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*(env->I_SR)[(i-j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*(env->I_SR)[(i-j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i-j)*N_seg_per_unit + 3][1];		}			
   	}	

   	// Compensate for the ANTERIOR segments with shorter processes   	
   	for(int i = env->N_SR-1; i >= 0; --i){ 
		I_SR_D_A[i] *= sqrt((env->N_SR/(i+1)));
		I_SR_V_A[i] *= sqrt((env->N_SR/(i+1)));	
   	}

	tmp_N_SR = env->N_SR;
   	for(int i = (N_units - env->N_SR + 1); i < N_units; ++i){			
		tmp_N_SR --;
		I_SR_D_B[i] = 0.0;
		I_SR_V_B[i] = 0.0;
		for(int j = 0; j < tmp_N_SR; ++j){
			I_SR_D_B[i] += env->I_SR[(i+j)*N_seg_per_unit][0] + (N_seg_per_unit >= 2)*(env->I_SR)[(i+j)*N_seg_per_unit + 1][0] + (N_seg_per_unit >= 3)*(env->I_SR)[(i+j)*N_seg_per_unit + 2][0]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i+j)*N_seg_per_unit + 3][0];	
			I_SR_V_B[i] += env->I_SR[(i+j)*N_seg_per_unit][1] + (N_seg_per_unit >= 2)*(env->I_SR)[(i+j)*N_seg_per_unit + 1][1] + (N_seg_per_unit >= 3)*(env->I_SR)[(i+j)*N_seg_per_unit + 2][1]  + (N_seg_per_unit >= 4)*(env->I_SR)[(i+j)*N_seg_per_unit + 3][1];
		}			
   	}	

   	// Compensate for the POSTERIOR segments with shorter processes   	
   	for(int i = (N_units - env->N_SR + 1); i < N_units; ++i){ 
		I_SR_D_B[i] *= sqrt(-(env->N_SR/(i-N_units)));
		I_SR_V_B[i] *= sqrt(-(env->N_SR/(i-N_units)));	
   	}

   	// Variables for total input current to each A- and B-class motorneuron
   	float I_D_A[N_units];
   	float I_V_A[N_units];
	float I_D_B[N_units];
   	float I_V_B[N_units];    

   	// Current bias to compensate for the fact that neural inhibition only goes one way
	//Time dependence
	float I_bias = (0.5-0.5*(env->I_bias_sign)*doubleSigmoid(timenow,env->I_bias_tStart,env->I_bias_tEnd,env->I_bias_slope));

	//Can prevent VD inhibition
	float tyramine = env->tyr_sign*doubleSigmoid(timenow,env->tyr_tStart,env->tyr_tEnd,env->tyr_slope); 

	/******** Combine AVB and AVA currents, stretch receptor current, neural inhibition and bias ********/

   	for(int i = 0; i < N_units; ++i){
		if(!(env->isInh_I_SR)){
			//NOTE: I_bias_D_inh=0 by default
			//Note that in the default program, the stretch receptor current can be inhibitory
			//	To compensate for this, there's a factor of two that would multiply the now strictly positive stretch receptor current
			//	Thus, the "ELSE" CONDITION HOLDS BY DEFAULT
			I_D_B[i] = (I_bias - env->State_B[i][1])*(env->I_bias_D_inh) + I_AVB[i] + SR_weight_B[i]*(env->SR_weight_factor_B)*fmax(I_SR_D_B[i],0.0);
			I_V_B[i] = (I_bias - env->State_B[i][0])*(1.0-tyramine) + I_AVB[i] + SR_weight_B[i]*(env->SR_weight_factor_B)*fmax(I_SR_V_B[i],0.0);
		
			I_D_A[i] = (I_bias - env->State_A[i][1])*(env->I_bias_D_inh) + I_AVA[i] + SR_weight_A[i]*(env->SR_weight_factor_A)*fmax(I_SR_D_A[i],0.0);
			I_V_A[i] = (I_bias - env->State_A[i][0])*(1.0-tyramine) + I_AVA[i] + SR_weight_A[i]*(env->SR_weight_factor_A)*fmax(I_SR_V_A[i],0.0);
		}
		else{
			I_D_B[i] = (I_bias - env->State_B[i][1])*(env->I_bias_D_inh) + I_AVB[i] + SR_weight_B[i]*I_SR_D_B[i];
			I_V_B[i] = (I_bias - env->State_B[i][0])*(1.0-tyramine) + I_AVB[i] + SR_weight_B[i]*I_SR_V_B[i];
			

			I_D_A[i] = (I_bias - env->State_A[i][1])*(env->I_bias_D_inh) + I_AVA[i] + SR_weight_A[i]*I_SR_D_A[i];
			I_V_A[i] = (I_bias - env->State_A[i][0])*(1.0-tyramine) + I_AVA[i] + SR_weight_A[i]*I_SR_V_A[i];
		}
   	}

	if(env->Gap_coupling>0.0){
		// Add gap junction currents if they are being used (typically Gap_coupling = 0)
		I_D_B[0] += (env->State_B[1][0] - env->State_B[0][0])*env->Gap_coupling;
		I_V_B[0] += (env->State_B[1][1] - env->State_B[0][1])*env->Gap_coupling;
		I_D_A[0] += (env->State_A[1][0] - env->State_A[0][0])*env->Gap_coupling;
		I_V_A[0] += (env->State_A[1][1] - env->State_A[0][1])*env->Gap_coupling;

		for(int i = 1; i < N_units-1; ++i){
			I_D_B[i] += ((env->State_B[i+1][0] - env->State_B[i][0]) + (env->State_B[i-1][0] - env->State_B[i][0]))*env->Gap_coupling;
			I_V_B[i] += ((env->State_B[i+1][1] - env->State_B[i][1]) + (env->State_B[i-1][1] - env->State_B[i][1]))*env->Gap_coupling;
			I_D_A[i] += ((env->State_A[i+1][0] - env->State_A[i][0]) + (env->State_A[i-1][0] - env->State_A[i][0]))*env->Gap_coupling;
			I_V_A[i] += ((env->State_A[i+1][1] - env->State_A[i][1]) + (env->State_A[i-1][1] - env->State_A[i][1]))*env->Gap_coupling;
		}

		I_D_B[N_units-1] += (env->State_B[N_units-2][0] - env->State_B[N_units-1][0])*env->Gap_coupling;
		I_V_B[N_units-1] += (env->State_B[N_units-2][1] - env->State_B[N_units-1][1])*env->Gap_coupling;
		I_D_A[N_units-1] += (env->State_A[N_units-2][0] - env->State_A[N_units-1][0])*env->Gap_coupling;
		I_V_A[N_units-1] += (env->State_A[N_units-2][1] - env->State_A[N_units-1][1])*env->Gap_coupling;
	}
   
	int id = omp_get_thread_num();

   	// Update state for each bistable B-class neuron
   	for(int i = 0; i < N_units; ++i){
   		if(I_D_A[i] > (0.5 + env->Hyst/2.0 - env->Hyst*env->State_A[i][0])){
			env->State_A[i][0] = 1;}
  		else{
			env->State_A[i][0] = 0;}

   		if(I_V_A[i] > (0.5 + env->Hyst/2.0 - env->Hyst*env->State_A[i][1])){
			env->State_A[i][1] = 1;}
   		else{
			env->State_A[i][1] = 0;}

		if(I_D_B[i] > (0.5 + env->Hyst/2.0 - env->Hyst*env->State_B[i][0])){
			env->State_B[i][0] = 1;}
  		else{
			env->State_B[i][0] = 0;}

   		if(I_V_B[i] > (0.5 + env->Hyst/2.0 - env->Hyst*env->State_B[i][1])){
			env->State_B[i][1] = 1;}
   		else{
			env->State_B[i][1] = 0;}
   	}

   	// Compute effective input to each muscle including A- and B-class excitation and contralateral env->D-class inhibition
	   //For now, use a weighting to show the percentage of each circuit  
	float forward_circuit;
	float backward_circuit;
	int thisN;

   	for(int i = 0; i < NSEG; ++i){
		thisN = (int)(i*N_units/NSEG); //Translate between neurons (12) and muscles (48)

		if(I_AVB[thisN]+I_AVA[thisN] > 0.1){ //Make sure it's not a 0 overall current!
			forward_circuit = I_AVB[thisN]/abs(I_AVB[thisN]+I_AVA[thisN]);
			backward_circuit = I_AVA[thisN]/abs(I_AVB[thisN]+I_AVA[thisN]);
	   }
	   else{
		   forward_circuit = 0;
		   backward_circuit = 0;
	   }

		NMJ_weight_A[i] *= backward_circuit;
		NMJ_weight_B[i] *= forward_circuit;

		// Add in a tyramine factor that might turn off VD (i.e. allow VB to activate for a long time)
		env->V_neuron[i][1] = NMJ_weight_A[i]*env->State_A[thisN][1] + NMJ_weight_B[i]*env->State_B[thisN][1] - tyramine*(NMJ_weight_A[i]*env->State_A[thisN][0] + NMJ_weight_B[i]*env->State_B[thisN][0]);

		env->V_neuron[i][0] = NMJ_weight_A[i]*env->State_A[thisN][0] + NMJ_weight_B[i]*env->State_B[thisN][0] - (NMJ_weight_A[i]*env->State_A[thisN][1] + NMJ_weight_B[i]*env->State_B[thisN][1]);
   	}


	//Save the SR_weight_B values for testing

	float FPS = 25;	// Frame rate at which output should be logged
  	// Save data to file
  	if(env->tCount_meta == 1.0/(DELTAT*FPS)){
  
  		outfile_meta << timenow;	
  		for(int i = 0; i < N_units; ++i){
			outfile_meta << ", " << SR_weight_B[i];
  		}
  		outfile_meta << "\n";
		env->tCount_meta = 0;
  	}
  	++env->tCount_meta;

  }

  // Update the stretch receptors (for each segment). These
  // are weighted and combined as input to the neural units in function "update_neurons"
  void update_SR(realtype timenow, UserEnviron* env){
   	for(int i = 0; i < NSEG; ++i){	
		// Bilinear SR function on one side to compensate for asymmetry and help worm go straight
		env->I_SR[i][0] = env->SR_shape_compensation[i]*((env->L_SR[i][0] - env->L0_P[i])/env->L0_P[i]*((env->L_SR[i][0] > env->L_seg) ? 0.8:1.2));
		env->I_SR[i][1] = env->SR_shape_compensation[i]*((env->L_SR[i][1] - env->L0_P[i])/env->L0_P[i]);
   	}
  }

  // Update the simple muscle "model" (electronic)
  void update_muscles(realtype timenow, UserEnviron* env){
  	//Muscle transfer function is just a simple LPF
  	for(int i = 0; i < NSEG; ++i){
		for(int j = 0; j < 2; ++j){
			realtype dV = (env->V_neuron[i][j] - env->V_muscle[i][j])/env->T_muscle;
			env->V_muscle[i][j] += dV*DELTAT;		
		}
  	}
  }

  // System residual function which implements physical model (Based on Sundials examples) 
  int resrob(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *rdata)
  {
  	// Import data from vectors
  	realtype *yval, *ypval, *rval;
  	yval = NV_DATA_S(yy); 
  	ypval = NV_DATA_S(yp); 
  	rval = NV_DATA_S(rr);

	//Declare user data struct
	UserEnviron* env = (UserEnviron*) rdata;

  	//Declare variables
  	realtype CoM[NBAR][3];
  	realtype V_CoM[NBAR][3];
  	realtype term[NBAR][2][2];  		// Nseg, env->D/v, x/y
  	realtype V_term[NBAR][2][2];
  	realtype dy,dx,dVy,dVx,F_even,F_odd;
  	realtype F_term[NBAR][2][2];
  	realtype F_term_rotated[NBAR][2][2];
  	realtype V_CoM_rotated[NBAR][3];

  	realtype L[NSEG][2];				// Nseg, env->D/v
  	realtype Dir[NSEG][2][2];			// Nseg, env->D/v, x/y
  	realtype S[NSEG][2];
  	realtype L_D[NSEG][2];			// Nseg, \,/   <- these are the angles of the diagonals
  	realtype Dir_D[NSEG][2][2];  			// Nseg, \,/ , x/y
  	realtype S_D[NSEG][2];

  	realtype L0_AE, T, F_AE, F_PE, F_PD;
  	realtype F_H[NSEG][2];  
  	realtype F_D[NSEG][2];

  	realtype L_EXT;
  	realtype F_EXT[2];

  	realtype F_object[NBAR][2][2];
	// Initialize all object forces to zero incase env->Objects are not being used
	for(int i = 0; i < NBAR; ++i){
		for(int j = 0; j < 2; ++j){
			for(int k = 0; k < 2; ++k){
				F_object[i][j][k] = 0;
			}
		}
	}
  	
  	for(int i = 0; i < NBAR; ++i){
		// Extract CoM of each solid rod from vectors
		int three_i = i*3;	
		CoM[i][0] = yval[three_i];
		CoM[i][1] = yval[three_i + 1];
		CoM[i][2] = yval[three_i + 2];

		// Calculate positions of env->D/V points based on CoM, angle and radius
		dx = env->R[i]*cos(CoM[i][2]);
		dy = env->R[i]*sin(CoM[i][2]);

		// Neighbor positions?
		term[i][0][0] = CoM[i][0] + dx;
		term[i][0][1] = CoM[i][1] + dy;
		term[i][1][0] = CoM[i][0] - dx;
		term[i][1][1] = CoM[i][1] - dy;

		// Extract CoM velocities of each solid rod from vectors
		V_CoM[i][0] = ypval[three_i];
		V_CoM[i][1] = ypval[three_i + 1];
		V_CoM[i][2] = ypval[three_i + 2];

		// Calculate velocity of env->D/V points based on CoM velocity, rate of rotation and radius
		realtype V_arm = env->R[i]*V_CoM[i][2];		
		dVx = V_arm*cos(CoM[i][2] + HALFPI);  		
		dVy = V_arm*sin(CoM[i][2] + HALFPI);  		

		V_term[i][0][0] = V_CoM[i][0] + dVx;
		V_term[i][0][1] = V_CoM[i][1] + dVy;
		V_term[i][1][0] = V_CoM[i][0] - dVx;
		V_term[i][1][1] = V_CoM[i][1] - dVy;		
  	}	

  
  	// Get Horizontal/Diagonal element lengths and lengthening/shortening velocities
  	for(int i = 0; i < NSEG; ++i){

		// Strange format for efficiency
		int iplus1 = i+1;
	
		//Calculate the difference between CoM points, using the neighborhood points calculated above ('term')
		Dir[i][0][0] = (term[iplus1][0][0] - term[i][0][0]);//Consecutive CoM points +dx
		Dir[i][0][1] = (term[iplus1][0][1] - term[i][0][1]);//"                    " +dy
		//Normalize the distances
		L[i][0] = sqrt(pow(Dir[i][0][0],2.0) + pow(Dir[i][0][1],2.0));
		Dir[i][0][0] /= L[i][0];
		Dir[i][0][1] /= L[i][0];
		//Difference is consecutive velocity neighborhood points multiplied by difference in consecutive CoM neighborhood points
		S[i][0] = (V_term[iplus1][0][0] - V_term[i][0][0])*Dir[i][0][0] + (V_term[iplus1][0][1] - V_term[i][0][1])*Dir[i][0][1];

		//Same as above but for -dx and -dy
		Dir[i][1][0] =  (term[iplus1][1][0] - term[i][1][0]);
		Dir[i][1][1] =  (term[iplus1][1][1] - term[i][1][1]);
		L[i][1] = sqrt(pow(Dir[i][1][0],2.0) + pow(Dir[i][1][1],2.0));
		Dir[i][1][0] /= L[i][1];
		Dir[i][1][1] /= L[i][1];
		S[i][1] = (V_term[iplus1][1][0] - V_term[i][1][0])*Dir[i][1][0] + (V_term[iplus1][1][1] - V_term[i][1][1])*Dir[i][1][1];

		//Same as two above blocks, but for the diagonal portion
		Dir_D[i][0][0] =  (term[iplus1][1][0] - term[i][0][0]);
		Dir_D[i][0][1] =  (term[iplus1][1][1] - term[i][0][1]);
		L_D[i][0] = sqrt(pow(Dir_D[i][0][0],2.0) + pow(Dir_D[i][0][1],2.0));
		Dir_D[i][0][0] /= L_D[i][0];
		Dir_D[i][0][1] /= L_D[i][0];
		S_D[i][0] = (V_term[iplus1][1][0] - V_term[i][0][0])*Dir_D[i][0][0] + (V_term[iplus1][1][1] - V_term[i][0][1])*Dir_D[i][0][1];

		Dir_D[i][1][0] =  (term[iplus1][0][0] - term[i][1][0]);
		Dir_D[i][1][1] =  (term[iplus1][0][1] - term[i][1][1]);
		L_D[i][1] = sqrt(pow(Dir_D[i][1][0],2.0) + pow(Dir_D[i][1][1],2.0));
		Dir_D[i][1][0] /= L_D[i][1];
		Dir_D[i][1][1] /= L_D[i][1];
		S_D[i][1] = (V_term[iplus1][0][0] - V_term[i][1][0])*Dir_D[i][1][0] + (V_term[iplus1][0][1] - V_term[i][1][1])*Dir_D[i][1][1];  

  		// Calculate force contributions on each env->D/V point

  		//Dorsal forces due to horizontal elements
		//If the muscle is fully activated, we are at the minimum length 
  		L0_AE = env->L0_P[i] - fmax(env->V_muscle[i][0],0)*(env->L0_P_minus_L_min[i]);  	
  		//Muscle activation spring term
		F_AE = env->k_AE*fmax(env->V_muscle[i][0],0)*(L0_AE - L[i][0]);
		//Passive spring terms: Include a ()^4 term if the length is greater than the rest length
  		F_PE = env->k_PE*((env->L0_P[i] - L[i][0]) + ((L[i][0]-env->L0_P[i]) > RCONST(0.0))*pow(RCONST(2.0)*(L[i][0]-env->L0_P[i]),4)); //NOTE: in the paper, this 2.0 is outside the ^4 term
  		//Damping term
		F_PD = (env->D_PE + fmax(env->V_muscle[i][0],0)*(env->D_AE))*S[i][0];

  		F_H[i][0] = F_PE + F_AE - F_PD;
	  
  		//Ventral forces due to horizontal elements
  		L0_AE = env->L0_P[i] - fmax(env->V_muscle[i][1],0)*(env->L0_P_minus_L_min[i]);  	
  		
  		F_AE = env->k_AE*fmax(env->V_muscle[i][1],0)*(L0_AE - L[i][1]);
  		F_PE = env->k_PE*((env->L0_P[i] - L[i][1]) + ((L[i][1]-env->L0_P[i]) > RCONST(0.0))*pow(RCONST(2.0)*(L[i][1]-env->L0_P[i]),4));
  		F_PD = (env->D_PE + fmax(env->V_muscle[i][1],0)*(env->D_AE))*S[i][1];

  		F_H[i][1] = F_PE + F_AE - F_PD;
	
  		//Diagonal forces due to diagonal elements
  		F_D[i][0] = (env->L0_D[i] - L_D[i][0])*(env->k_DE) - env->D_DE*S_D[i][0];
  		F_D[i][1] = (env->L0_D[i] - L_D[i][1])*(env->k_DE) - env->D_DE*S_D[i][1];
  	}

  	// If using env->Objects, check for object collisions and calculate associated forces
	if(env->N_objects > 0){
  		realtype P_x,P_y,Distance,magF,D_scale,magF_P1,magF_P2;
  		env->ContactForce = 0;
  		for(int i = 0; i < NBAR; ++i){
			for(int j = 0; j < 2; ++j){
				// First ensure they contain zeros
				F_object[i][j][0] = 0;
				F_object[i][j][1] = 0;
				P_x = term[i][j][0];
				P_y = term[i][j][1];

				// Now check proximity to each object
				for(int k = 0; k < env->N_objects; ++k){
					if((P_x<(env->Objects[k][0]+env->Objects[k][2]))&&(P_x>(env->Objects[k][0]-env->Objects[k][2]))&&(P_y<(env->Objects[k][1]+env->Objects[k][2]))&&(P_y>(env->Objects[k][1]-env->Objects[k][2]))){

						//This means the point is within the bounding box of the object, so now we must compute the force (if any)
						dx = P_x - env->Objects[k][0];
						dy = P_y - env->Objects[k][1];
						Distance = sqrt(pow(dx,2) + pow(dy,2));
						D_scale = 0.01*(env->Objects)[k][2];

						if(Distance < env->Objects[k][2]){
							magF = env->k_Object*(env->Objects[k][2] - Distance) + D_scale*(env->k_Object)*(pow((env->Objects[k][2] - Distance)/D_scale,2));
							F_object[i][j][0] += (dx/Distance)*magF;
							F_object[i][j][1] += (dy/Distance)*magF;
							env->ContactForce += magF;
						}
					}			
				}
			}
  		}
	}

  	// Add up force contributions for each env->D/V point
  	F_term[0][0][0] = -F_H[0][0]*Dir[0][0][0] - F_D[0][0]*Dir_D[0][0][0] + F_object[0][0][0];
  	F_term[0][0][1] = -F_H[0][0]*Dir[0][0][1] - F_D[0][0]*Dir_D[0][0][1] + F_object[0][0][1];

  	F_term[0][1][0] = -F_H[0][1]*Dir[0][1][0] - F_D[0][1]*Dir_D[0][1][0] + F_object[0][1][0];
  	F_term[0][1][1] = -F_H[0][1]*Dir[0][1][1] - F_D[0][1]*Dir_D[0][1][1] + F_object[0][1][1];

  	for(int i = 1; i < NSEG; ++i){
		int i_minus_1 = i-1;

		F_term[i][0][0] = F_H[i_minus_1][0]*Dir[i_minus_1][0][0] - F_H[i][0]*Dir[i][0][0] + F_D[i_minus_1][1]*Dir_D[i_minus_1][1][0] - F_D[i][0]*Dir_D[i][0][0] + F_object[i][0][0];
		F_term[i][0][1] = F_H[i_minus_1][0]*Dir[i_minus_1][0][1] - F_H[i][0]*Dir[i][0][1] + F_D[i_minus_1][1]*Dir_D[i_minus_1][1][1] - F_D[i][0]*Dir_D[i][0][1] + F_object[i][0][1];

		F_term[i][1][0] = F_H[i_minus_1][1]*Dir[i_minus_1][1][0] - F_H[i][1]*Dir[i][1][0] + F_D[i_minus_1][0]*Dir_D[i_minus_1][0][0] - F_D[i][1]*Dir_D[i][1][0] + F_object[i][1][0];
		F_term[i][1][1] = F_H[i_minus_1][1]*Dir[i_minus_1][1][1] - F_H[i][1]*Dir[i][1][1] + F_D[i_minus_1][0]*Dir_D[i_minus_1][0][1] - F_D[i][1]*Dir_D[i][1][1] + F_object[i][1][1];
  	}

  	F_term[NSEG][0][0] = F_H[NSEG_MINUS_1][0]*Dir[NSEG_MINUS_1][0][0] + F_D[NSEG_MINUS_1][1]*Dir_D[NSEG_MINUS_1][1][0] + F_object[NSEG][0][0];
  	F_term[NSEG][0][1] = F_H[NSEG_MINUS_1][0]*Dir[NSEG_MINUS_1][0][1] + F_D[NSEG_MINUS_1][1]*Dir_D[NSEG_MINUS_1][1][1] + F_object[NSEG][0][1];

  	F_term[NSEG][1][0] = F_H[NSEG_MINUS_1][1]*Dir[NSEG_MINUS_1][1][0] + F_D[NSEG_MINUS_1][0]*Dir_D[NSEG_MINUS_1][0][0] + F_object[NSEG][1][0];
  	F_term[NSEG][1][1] = F_H[NSEG_MINUS_1][1]*Dir[NSEG_MINUS_1][1][1] + F_D[NSEG_MINUS_1][0]*Dir_D[NSEG_MINUS_1][0][1] + F_object[NSEG][1][1];
  
  	// Convert net forces on env->D/V points to force and torque	acting on rod CoM
  	for(int i = 0; i < NBAR; ++i){
		realtype cos_thi = cos(CoM[i][2]);
		realtype sin_thi = sin(CoM[i][2]);
		for(int j = 0; j < 2; ++j){
			//Note: all *_rotated terms are in radial coordinates	
			F_term_rotated[i][j][0] = F_term[i][j][0]*cos_thi + F_term[i][j][1]*sin_thi;	// This is Fperp
			F_term_rotated[i][j][1] = F_term[i][j][0]*sin_thi - F_term[i][j][1]*cos_thi;    // THis is Fparallel
		}

		V_CoM_rotated[i][0] = (F_term_rotated[i][0][0] + F_term_rotated[i][1][0])/env->CN[i]; //Perp forces sum simply

		F_even = (F_term_rotated[i][0][1] + F_term_rotated[i][1][1]);	//Took out the /2 ...//Linear motion
		F_odd = (F_term_rotated[i][1][1] - F_term_rotated[i][0][1])/RCONST(2.0);	//TORQUE

		V_CoM_rotated[i][1] = (F_even)/env->CL[i];				//Allowing me to take out *2
		V_CoM[i][2] = (F_odd/env->CL[i])/(M_PI*2.0*(env->R)[i]);

		V_CoM[i][0] = V_CoM_rotated[i][0]*cos_thi + V_CoM_rotated[i][1]*sin_thi;
		V_CoM[i][1] = V_CoM_rotated[i][0]*sin_thi - V_CoM_rotated[i][1]*cos_thi;
	
		int three_i = i*3;

		rval[three_i] = V_CoM[i][0] - ypval[three_i];
		rval[three_i+1] = V_CoM[i][1] - ypval[three_i+1];
		rval[three_i+2] = V_CoM[i][2] - ypval[three_i+2];
  	}

  	// Store old lengths for Stretch Receptors 
  	for(int i = 0; i < NSEG; ++i){
		env->L_SR[i][0] = L[i][0];
		env->L_SR[i][1] = L[i][1];
  	}

  	return(0);
  }


  float doubleSigmoid(float timenow, float tstart, float tend, float slope)
  {
	  //This function impliments a float sigmoidal function (with a tanh approximation), 
	  //shifting from 'OFF' at tstart to 'ON' to back 'OFF' at tend

	if(slope<0.001){
		return 0.0;
	}
	else{
		float tstartSig = tstart + 2.0/slope; //We want the CENTER of the first sigmoid to be after tstart
		float tendSig   = tend - 2.0/slope; //Similarly with the second sigmoid
		
		//Sum of an inverted and a regular sigmoid (approximated by tanh)
		return( 0.5*( tanh( slope*(timenow-tstartSig )) ) - 0.5*tanh( slope*(timenow-tendSig) ) );
	}
  }



/*
 *--------------------------------------------------------------------
 * Function for writing prams files
 *--------------------------------------------------------------------
 */

int write_prams_file(int n, ...)
{
	//pramsFile.open(pramsFileName.c_str());

	va_list userPrams;
	va_start(userPrams,n);

	char* argName, *argValue;
	argName = va_arg(userPrams, char*);
	argValue = va_arg(userPrams, char*);

	if(strcmp("FILENAME",argName) == 0){
		//NOTE: doing this the C way because the va_arg function wants me to
		char* extStr = (char *)".csv";
		char* preStr = (char *)"prams_";
		char* pramsFullFileName = (char *) malloc(1 + strlen(preStr) + strlen(argValue)+ strlen(extStr) );
		strcpy(pramsFullFileName, preStr);
		strcat(pramsFullFileName, argValue);
		strcat(pramsFullFileName, extStr);

		FILE *pramsFile;

		if (ifstream(pramsFullFileName))
		{
     		cout << "Parameter file already exists; using old file." << endl;
    	 	return 2;
		}
		else{
			pramsFile = fopen(pramsFullFileName, "w");
			cout << "Writing file: " << pramsFullFileName << endl;
		}

		fprintf(pramsFile,"%s ",argName);
		fprintf(pramsFile,"%s",argValue);
		fprintf(pramsFile,"\r\n");

		for(int j=0;j<(n/2-1);++j)
		{
			argName = va_arg(userPrams, char*);
			argValue = va_arg(userPrams, char*);
			
			fprintf(pramsFile,"%s ",argName);
			fprintf(pramsFile,"%s",argValue);
			fprintf(pramsFile,"\r\n");
		}

		va_end(userPrams);
		fclose(pramsFile);
	}
	else{
		cout << "Error: the second argument must be 'FILENAME' and the third must be the file's name.\n Not:" << argName << ' ' << argValue << endl;
		return(1);
	}

	return(0);

}



/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */
  // Save output to csv file for visualization
  static void SaveOutput(void *mem, realtype t, N_Vector y, ofstream& outfile, UserEnviron* env){
  
  	float FPS = 25;	// Frame rate at which output should be logged
  	realtype *yval;
  	int retval, kused;
  	long int nst;
  	realtype hused;

  	yval  = NV_DATA_S(y);
  
  	// Save data to file
  	if(env->tCount == 1.0/(DELTAT*FPS)){
  		outfile << t;	
  		for(int i = 0; i < NBAR; ++i){
			outfile << ", " << yval[i*3] << ", " << yval[i*3+1] << ", " << yval[i*3+2];
  		}
  		outfile << "\n";
		env->tCount = 0;
  	}

  	++env->tCount;
  	return;  
  }

  double randn(double mu, double sigma) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double polar, rsquared, var1, var2;
	
	//	If no deviate has been stored, the polar Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose pairs of uniformly distributed deviates, discarding those 
		//	that don't fall within the unit circle
		do {
			var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);
		
		//	calculate polar tranformation for each deviate
		polar=sqrt(-2.0*log(rsquared)/rsquared);
		
		//	store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=true;
		
		//	return second deviate
		return var2*polar*sigma + mu;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
  }



  // Print final integrator statistics (Copied from Sundials examples)
  static void PrintFinalStats(void *mem){
  	int retval;
  	long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  	retval = IDAGetNumSteps(mem, &nst);
  	check_flag(&retval, (char *)"IDAGetNumSteps", 1);
  	retval = IDAGetNumResEvals(mem, &nre);
  	check_flag(&retval, (char *)"IDAGetNumResEvals", 1);
  	//retval = IDADenseGetNumJacEvals(mem, &nje);
  	//check_flag(&retval, "IDADenseGetNumJacEvals", 1);
  	retval = IDAGetNumNonlinSolvIters(mem, &nni);
  	check_flag(&retval, (char *)"IDAGetNumNonlinSolvIters", 1);
  	retval = IDAGetNumErrTestFails(mem, &netf);
  	check_flag(&retval, (char *)"IDAGetNumErrTestFails", 1);
  	retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  	check_flag(&retval, (char *)"IDAGetNumNonlinSolvConvFails", 1);
  	//retval = IDADenseGetNumResEvals(mem, &nreLS);
  	//check_flag(&retval, "IDADenseGetNumResEvals", 1);
  	retval = IDAGetNumGEvals(mem, &nge);
  	check_flag(&retval, (char *)"IDAGetNumGEvals", 1);

  	printf("\nFinal Run Statistics: \n\n");
  	printf("Number of steps                    = %ld\n", nst);
  	printf("Number of residual evaluations     = %ld\n", nre+nreLS);
  	printf("Number of Jacobian evaluations     = %ld\n", nje);
  	printf("Number of nonlinear iterations     = %ld\n", nni);
  	printf("Number of error test failures      = %ld\n", netf);
  	printf("Number of nonlinear conv. failures = %ld\n", ncfn);
  	printf("Number of root fn. evaluations     = %ld\n", nge);
  }

/*
 * Check function return value... (Copied from Sundials examples)
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, 
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, 
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", 
              funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    fprintf(stderr, 
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n", 
            funcname);
    return(1);
  }

  return(0);
}
