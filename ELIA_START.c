#include <windows.h>
#include <stdio.h>
#include <tchar.h>
#include <conio.h>
#include <string.h>
#include <math.h>
#include<io.h>

// adress of executed program
#define FILEPATH "d:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\ELIMA_programming\\ELiMA_git\\ELiMA\\ELiMA\\ELiMA\\x64\\Debug\\ELiMA.exe"
//#define FILEPATH "D:\\01_Khelemelia\\02_IAP\\07_Programming\\01_MPI_C\\01_MPI_Projects\\HelloWorld\\x64\\Debug\\HelloWorld.exe"
										
// procedure dataFileWrite
// income:
	// file_name -- name of file
	// parameterCount -- count of parameters
	// parametersName -- names of parameters
	// parameterValue -- values of parameters
	// msg -- pointer of message
// outcome:
	// msg -- message of result of execution
// produce:
	// names of parameters are written to file 
	// values of parameters are written to file 
// return 1  - the work have been interrupted 
// return 0  - the work finished normally 
int dataFileWrite(char *fileName, int parameterCount, char** parametersName, double *parameterValue, char *msg[]);

// procedure dataFileWrite
// income:
	// unionVector -- array of values union vector
	// vectorValue -- absolute value of vector
	// vector -- poiner of vector
	// msg -- pointer of message
// outcome:
	// vvector -- value of vector components  
	// msg -- message of result of execution
// produce:
	// value of vector components has been written to array "vector"
// return 1  - the work have been interrupted 
// return 0  - the work finished normally 
int vectorDefine(double unionVector[3], double vectorValue, double *vector, char *msg[]);

// procedure vectorDirection
// income:
	// angle between vector direction and x-axis
	// unionVector -- pointer of array of values union vector
	// msg -- pointer of message
// outcome:
	// unionVector -- value of union vector components  
	// msg -- message of result of execution
// produce:
	// transform value jf angle to union vector components
// return 1  - the work have been interrupted 
// return 0  - the work finished normally 

int vectorDirection(const int angle, double *unionVector, char *msg[]);

int main( int argc, TCHAR *argv[] )
{
	int i;
	double velocityValue = 1.2;
	double magnFieldValue = 0.;
	const int angleVelocity = 0;	// set(0, 30, 45, 60, 90)
	const int angleMagnField = 0; // set(0, 30, 45, 60, 90)
	double tempPerp = 0.001;
	double tempParall = 0.001;
	char *msg;						// pattern for message of result of execution
	int countParmtrs = 10;			// count of task parameters
	double *arrayParmtrValue;		// array of task parameters values
	char **arrayParmtrName;			// array of task parameters names
	double unionVelVector[3];		// union velocity vector component
	double vectorVelIon[3];			// array for velocity vector components
	double unionMagnVector[3];		// union magnetic vector component
	double vectorMagnStrength[3];	// array for velocity vector components

	/*char *adrDateFile = "D:\\01_Khelemelia\\02_IAP\\07_Programming\\01_MPI_C\\01_MPI_Projects\\HelloWorld\\data.dat";*/
	char* adrDateFile = "d:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\ELIMA_programming\\ELIMA_git\\ELiMA\\ELiMA\\data.dat";

										//file to record parameters  xx
	char *parmtrsExecProg = " -n 4 ";	// parameters of exec program;
	char* adrExecProg = "mpiexec";		// execute program
	char commandLine[1024] = "";		// command line

	arrayParmtrValue = (double*) malloc(countParmtrs*sizeof(double));
	arrayParmtrName = (char**) malloc(countParmtrs*6*sizeof(char)); 
									// less than 6 chars for abbreviation of the parameter 
	
	// assign parameter names
	arrayParmtrName[0] = "Vx";
	arrayParmtrName[1] = "Vy";
	arrayParmtrName[2] = "Vz";
	arrayParmtrName[3] = "|V|";
	arrayParmtrName[4] = "Te_per";
	arrayParmtrName[5] = "Te_par";
	arrayParmtrName[6] = "Hx";
	arrayParmtrName[7] = "Hy";
	arrayParmtrName[8] = "Hz";
	arrayParmtrName[9] = "|H|";

	// Union velocity vector nV definition
	vectorDirection(angleVelocity, unionVelVector, &msg);

	// Union magnetic field strength vector nH definition
	vectorDirection(angleMagnField, unionMagnVector, &msg);

	// define components of velocity vector 
	if(vectorDefine(unionVelVector, velocityValue, vectorVelIon, &msg)) { // Error 
		fprintf(stdout,"%s", msg);
		_getch();
		return 1;
	}
	// define components of magnetic field strength vector 
	if (vectorDefine(unionMagnVector, magnFieldValue, vectorMagnStrength, &msg)) { // Error 
		fprintf(stdout, "%s", msg);
		_getch();
		return 1;
	}

	// pack parameters
	arrayParmtrValue[0] = vectorVelIon[0];		// Vx
	arrayParmtrValue[1] = vectorVelIon[1];		// Vy
	arrayParmtrValue[2] = vectorVelIon[2];		// Vz
	arrayParmtrValue[3] = velocityValue;		//|V|
	arrayParmtrValue[4] = tempPerp;				// Te_perp, eV
	arrayParmtrValue[5] = tempParall;			// Te_parallel, eV
	arrayParmtrValue[6] = vectorMagnStrength[0];// Hx
	arrayParmtrValue[7] = vectorMagnStrength[1];// Hy
	arrayParmtrValue[8] = vectorMagnStrength[2];// Hz
	arrayParmtrValue[9] = magnFieldValue;		//|H|
	
	//display parameters
	for(i = 0; i < countParmtrs; i++) {
		fprintf(stdout, "%s\t%e\n", arrayParmtrName[i], arrayParmtrValue[i]);
	}

	//write parameters to a file
	if (dataFileWrite(adrDateFile, countParmtrs, arrayParmtrName, arrayParmtrValue, &msg)) {
		fprintf(stdout,"%s", msg);
		_getch();
		return 1;
	}

	//create comandline
	strcat(commandLine, adrExecProg);
	strcat(commandLine, parmtrsExecProg);
	strcat(commandLine, FILEPATH);


	//fprintf(stdout, "%s\n", commandLine);
	
	//run executed program
	system(commandLine);
	
	//clear memory
	free(arrayParmtrName);
	free(arrayParmtrValue);
	
	// signal of the end
	fprintf(stdout, "\a");fflush(stdout);
	
	// wait a tap 
	_getch();
}

int dataFileWrite(char *fileName, int countParameter, char** parametersName, double *parameterValue, char *msg[]) {
	FILE *dataFile;
	int i;

	if(fopen_s(&dataFile, fileName, "w+")){ // ERROR of creating of a file
		*msg = "dataFileWrite(Errorcode_0AC): Can't create datafile\n";
		return 1;
	}
	for (i = 0; i < countParameter; i++) {
		fprintf(dataFile, "%s\t", parametersName[i]);
	}
	fprintf(dataFile, "\n");
	for (i = 0; i < countParameter; i++) {
		fprintf(dataFile, "%.6e\t", parameterValue[i]); 
	}											
	fprintf(dataFile, "\n"); 
	fclose(dataFile);
	*msg = "dataFileWrite: program finished OK\n";
	return 0;
}

int vectorDefine(double unionVector[3], double vectorValue, double *vVector, char *msg[]) {
	int i;
	double precision = 0.000000001;
	double valueUnionV = 0;
	for (i = 0; i < 3; i++) {
		valueUnionV += unionVector[i]*unionVector[i];
	}
	if ((fabs(sqrt(valueUnionV) - 1) > precision)) {
		*msg = "vectorDefine (Errorcode_0BA): Uncorrect value of parameters. Union velocity vector\n";
		return 1;
	}
	if (vectorValue < 0) {
		*msg = "vectorDefine (Errorcode_0BB): Uncorrect value of parameters. Value of velocity vector\n";
		return 1;
	}
	for (i = 0; i < 3; i++) {
		vVector[i] = unionVector[i]*vectorValue;
	}
	*msg = "vectorDefine: program finished ok";
	return 0;
}

int vectorDirection(const int angle, double *unionVector, char *msg[]) {
	if (angle < 0 || angle > 180) {
		*msg = "vectorDirection(Errorcode_0CA): Uncorrect value of the angle parameter.\n";
		return 1;
	}
	switch (angle) {
		case 0: // parallel
			unionVector[0] = 0.;			//nHx
			unionVector[1] = 0.;			//nHy
			unionVector[2] = 1.;			//nHz
			break;
		case 30:// angle = 30 degrees
			unionVector[0] = 0.;			//nHx
			unionVector[1] = sqrt(3.) / 2.; //nHy
			unionVector[2] = sqrt(1.) / 2.; //nHz
			break;
		case 45:// angle = 45 degrees
			unionVector[0] = 0.;			//nHx
			unionVector[1] = sqrt(2.) / 2.; //nHy
			unionVector[2] = sqrt(2.) / 2.; //nHz
			break;
		case 60:// angle = 30 degrees
			unionVector[0] = 0.;			//nHx
			unionVector[1] = sqrt(1.) / 2.; //nHy
			unionVector[2] = sqrt(3.) / 2.; //nHz
			break;
		case 90:// angle = 30 degrees
			unionVector[0] = 0.;			//nHx
			unionVector[1] = 1;				//nHy
			unionVector[2] = 0;				//nHz
			break;
		default:
			*msg = "vectorDirection(Errorcode_0CA): Value of the angle parameter is not present in set.\n";
			return 1;
	}
	*msg = "vectorDirection: program finished ok";
	return 0;
}