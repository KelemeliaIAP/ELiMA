//////////////////////////////////////////////////////////////////////
// Calculate integral of function ELiMA.
// ELiMA is contained in 
// 1) function.c: double delta (TePerp) need to change parameters to the HESR experiment
///////////////////////////////////
/// Integrand is taken from:
///  Akhiezer. Sov.Phys.JETP		//
/// //////////////////////////////////////////////////////////////////////////////////////////////
/// Equation is three dimantional integral
///////////////////////////////////////////////////////////////////
/// athor: O.V.Khelemelia
/// project start date: 17.03.2016
///////////////////////////////////////////////////////////////////
/// You need to do next steps to include the "mpi.h" libruary:
/// 1) project->properties->VC++ Directories: in Include Derictories add adress to mpich2/include
/// 2) project->properties->VC++ Directories: in Library Derictories add adress to mpich2/lib
/// 3) project->properties->linker->input: in Additional Dependencies add mpi.h
/// 4) sometimes 
////	project->properties->lincer->general->Incremental(NO)

///////////////////////////////////////////////////////////////////
// 1) timing is missed
// 2) parseFile function in readNLine.c work uncorrect: msg veriable cannot be assigned

#include "MPIIntegration.h"
#include "Boundaries.h"
#include "Function.h"
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <stdlib.h> //malloc
#include <time.h>
#include <string.h> // atof
 

# ifndef PRECISION
#define PRECISION 0.0001
#endif

// numbers of task parameters
#define COUNTPARMTRS 10	// 3 - the union 3-dim velocity vector
			// 1 - an absolutly value of the velocity
			// 2 - the 2-dim electron Maxwell temperature (parallel and perpendicular component)
			// 3 - the union 3-dim magnetic field vector
			// 1 - an absolutly value of the magnetic field
// Program name
#define PROG_NAME "ELiMA"
// Program version
#define PROG_VERSION ".v.7.01"
// format of file for calculated data
#define FILE_FORMAT ".dat"
// Names od income parameters
// Tranverse temperature
#define TEMPERATURAPer "Te+"
// Longitudinal temperature
#define TEMPERATURAPar "Te||"
// Magnetic field strength
#define STRENGTHMagnetic "H"
// The angle of velocity of ions according to the direction of the electron beam 
#define ANGLEVelocity "aV"
// The angle of magnetic field according to the direction of the electron beam 
#define ANGLEMagnetic "aH"

// Output file is written here:
#define OUT_FILEPATH "D:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\MPI_Results\\"
#define IN_FILEPATH "D:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\ELIMA_programming\\ELiMA_git\\ELiMA\\ELiMA\\data.dat"
#define PARAMETERS_COUNT 10
#define PARAMETERS_NAME_LENGTH 6

int velocityDefine(double unionvector[3], double vectorvalue, double *vvector, char *msg[]);

//!!!add count of external parameters in data.dat as income argue
int main (int argc, char **argv) {
	
	/* LOGation of performing */
	// log file declaration
	FILE *LOGfile = NULL;
	// log message declaration
	char** msg;				// a message informs how the program works
	// locate memory on log message
	msg = (char**)malloc(100 * sizeof(char));

	/*Timing of performing*/
	struct tm * tmptime;

	/*MPI rutine initialization*/
	int MPIErrorCode = 1;

	// MPI initialization
	MPI_Init(&argc, &argv); 

	// test MPI initialization
	int flagMPI_Init;
	MPI_Initialized(&flagMPI_Init);
	
	if (!flagMPI_Init) {
		fprintf(stderr, "CUSTOM ERROR: MPI haven't been initializated. \n"); fflush(stdout);
		fprintf(stderr, "Failed in file %s at line # %d\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}

	// Create a new custom MPI communicator (unic for a each new running of procedure)
	MPI_Comm  comm;
	MPI_Group groupWorld;
	char* commName = "customComm";

	MPIErrorCode = MPI_Comm_group(MPI_COMM_WORLD, &groupWorld);
	MPIErrorCode = MPI_Comm_create(MPI_COMM_WORLD, groupWorld, &comm);
	MPIErrorCode = MPI_Comm_set_name(comm, commName);

	// Assigning an id (rank) for coprocessors
	int processorRank, processorSetSize;
	//define a number of work processors
	MPI_Comm_size(comm, &processorSetSize);
	//define indexes of work processors
	MPI_Comm_rank(comm, &processorRank); // index number of processor

	// Program is valid when more than one processors
	if (processorRank == 0) {
		if(processorSetSize <= 1) {
			MPI_Finalize();
			fprintf(stderr, "CUSTOM ERROR: Number of processes have to be more than 1. Calculation have been performed for multiprocessors\n"); fflush(stdout);
			fprintf(stderr, "Failed in file %s at line # %d\n", __FILE__, __LINE__);
			exit(EXIT_FAILURE);
		}
	}

	////////////////////////////
	// Главная часть программы//
	////////////////////////////

	//Greeting
	if (processorRank == 0) {
		fprintf(stdout, "Welcome to our program!\n");
		fprintf(stdout, "Calculation is performed on %d processors within %s .\n", processorSetSize, commName);
	}

	
//	int series = 1;			// index of program's start

	// open and read file
		// create pointers
		// char** headers
		// double **parameters
		//locate memory on tmp message

	char** setParmtrsName;		//pointer to an array with parameter names
	//locate memory on parameter names array
	setParmtrsName = (char**)malloc(PARAMETERS_COUNT * PARAMETERS_NAME_LENGTH * sizeof(char));

	double* setParmtrsValue;	//pointer to an array with parameter values
	//locate memory on parameter values array
	setParmtrsValue = (double*)malloc(PARAMETERS_COUNT * sizeof(double));

	if (processorRank == 0) {
		FILE* InFile = NULL;
		if (fopen_s(&InFile, IN_FILEPATH, "r")) { // true is error of opening
			fprintf(stdout, "File is not open");
			// don't forget to get free all pointers
			MPI_Finalize();
			exit(EXIT_FAILURE);
		}
		parseFile(&InFile, PARAMETERS_COUNT, &setParmtrsName, &setParmtrsValue, &msg);
		for (int i = 0; i < PARAMETERS_COUNT; i++) {
			fprintf(stdout, "%s\t", setParmtrsName[i]);
		}
		fprintf(stdout, " \n");
		for (int i = 0; i < PARAMETERS_COUNT; i++) {
			fprintf(stdout, "%.2e", setParmtrsValue[i]);
		}
		fprintf(stdout, " \n");
	}


	////create new MPI type to send and receive parametrs
	//MPI_Datatype sndrcvparmtrs;					//name of new MPI data type
	//MPI_Datatype oldtypes[1 + 1];	//array of old parameters will be included in the new MPI data type
	//MPI_Aint offsets[1+1];		//array of blocks in the new MPI data type
	//int blklens[1 + 1];				//array of element counts in each block
	//
	////old data types
	//oldtypes[0] = MPI_INT;									
	//oldtypes[1] = MPI_DOUBLE;

	////number of elements in each block					
	//blklens[0] = 1;														
	//blklens[1] = PARAMETERS_COUNT;
	//
	////offset of each block calculated in bytes					
	//offsets[0] = 0;														
	//MPI_Type_size(MPI_INT, &offsets[1]);						

	////creating new MPI structure
	//MPI_Type_create_struct(1 + 1, blklens, offsets, oldtypes, &sndrcvparmtrs);		

	//// commiting new MPI data type
	//MPI_Type_commit(&sndrcvparmtrs);										


	// run integration procedure

	// summ subresults from coprocessors
		// send all subresults to root coprocessor
		// summ results by 

	// delete pointers
		// *headers, *parameters
	// 


	//if (processorRank == 0) {
	//	FILE* inputFile;
	//}
	
	//{
	//	unionvectorV[0] = atof(argv[2]); //nVx
	//	unionvectorV[1] = atof(argv[3]); //nVy
	//	unionvectorV[2] = atof(argv[4]); //nVz

	//	if(velocityDefine(unionvectorV, atof(argv[5]), vectorV, msg)) { // error
	//		fprintf(stdout,"%s", msg); 
	//		return 1;
	//	}

	//	set_parmtrs[0] = vectorV[0];	// Vx
	//	set_parmtrs[1] = vectorV[1];	// Vy
	//	set_parmtrs[2] = vectorV[2];	// Vz
	//	set_parmtrs[3] = atof(argv[5]);	// |V|
	//	set_parmtrs[4] = atof(argv[6]);	// Te_per
	//	set_parmtrs[5] = atof(argv[7]);	// Te_par

	//// !!!!не забути підключити ці параметри до розрахунків у наступних файлах!!!!
	//	set_parmtrs[6] = delta(set_parmtrs[4]);	// Delta
	//	set_parmtrs[7] = tau_perpendicular();	
	//	set_parmtrs[8] = tau_parallel(atof(argv[7]), atof(argv[6])); //(Te_par, Te_per)
	//
	//}


	//////////////////////////////////////////////////////////////////
	// Call intagration function
	//////////////////////////////////////////////////////////////////
	//fprintf(stdout, "Wellcome to integation\n");
	double result = 0.;		// result of calculation

		//PreIntegration(rank, size, series, set_parmtrs, &result);
		//////////////////////////////////////////////////////////////////
//
//
//		/// Открытие файла для вывода результатов вычислений
//		if (processorRank == 0) {
//			double temperature_coeficient; 
//			// создание имени файла, который будет 
//			// содержать результаты вычислений
//			char *file_name;			// массив под имя файла
//			int length_name;			// длинна имени файла
//			// Result of calculation are written to the file
//			FILE *OUTResult;		// Выводится результат вычислений
//
//			// V_e/V_o = sqrt(T_e * C^2/(m_eC^2)) = sqrt(T_e) * 41.9382450;
//			// [V_e/V_o]/[T_e] = cm/eV/s
//			temperature_coeficient = 41.94;
//			// Длинна имени файла
////			length_name = sizeof(PROG_NAME)+sizeof(PROG_VERSION)+sizeof(FILE_FORMAT)
////							+sizeof(TEMPERATURAPer)+sizeof(TEMPERATURAPar)+sizeof(double)+sizeof(FILEPATH);
//			length_name = sizeof(PROG_NAME)+sizeof(PROG_VERSION)+sizeof(FILE_FORMAT)
//							+sizeof(TEMPERATURAPer)+sizeof(TEMPERATURAPar)+sizeof(double);
//			// выделение памяти под массив, содержащий имя файла
//			file_name = (char*)malloc(length_name);
//			// образование имени файла
//			
//			fprintf(file_name, "%s%s%s%g%s%g%s", PROG_NAME, PROG_VERSION, 
//							 TEMPERATURAPer, (double)set_parmtrs[4], TEMPERATURAPar, (double)set_parmtrs[5], FILE_FORMAT);
//			
//			fopen_s(&OUTResult, file_name, "a+");
//			if (!OUTResult) {
//				int MPIerrorCode = 1; // File was not opened
//				fprintf(stdout, "ERROR: Can't open OUTResult-file\n");fflush(stdout);
//				MPI_Abort(MPI_COMM_WORLD, MPIerrorCode);
//			}
//			fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
//			fprintf(LOGfile, "ResultFile was opened successfully\n");fflush(stdout);
//
////			// Vi/V0	\t	Te \t	Delta \t	Tau \t 
////			// result \t precision \t program \t 
////			// bouX_ \t bouX^ \t bouK_ \t bouK^\n", 
//			fprintf(OUTResult, "%g\t%g\t%g\t%g\t %g\t%g\t%g\t %g\t%g\t%g\t%s%s\t  %g\t%g\t%g\t%g\t%g\t%g\n", 
//						set_parmtrs[0], set_parmtrs[1], set_parmtrs[2], set_parmtrs[3],
//						set_parmtrs[3]*temperature_coeficient*sqrt(set_parmtrs[4]), set_parmtrs[4], set_parmtrs[5], 
//						result, result/(temperature_coeficient*sqrt(set_parmtrs[4])), PRECISION, PROG_NAME, PROG_VERSION, 
//						LOWX, UPX, LOWY, UPY, LOWZ, UPZ);fflush(OUTResult);
//								
//			fclose(OUTResult);
//			if(!OUTResult) {
//				fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
//				fprintf(LOGfile, "File was not closed\n");		
//				fprintf(stdout, "File was not closed\n");		
//			} else {
//				fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
//				fprintf(LOGfile, "ResultFile was closed successfully\n");fflush(stdout);
//			}
//
//		}
//	/// Закрытие файлов для вывода (log-file, result-file)
//	if(rank == 0) {
//		if(LOGfile) {
//			time_t tt = time(NULL);
//			struct tm ptm;
//			char buf[26];
//			asctime_s(buf, sizeof buf, localtime_s(&tt, &buf));
//			printf("Current local time and date: %s", buf);
//			
//			fprintf(LOGfile, "%s", asctime_s(buf, sizeof buf, localtime_s(&tt, &buf)));
//			fprintf(LOGfile, "END\n");		
//			fclose(LOGfile);
//		}
//	}
//
//	fprintf(stdout, "proc %d finalize\n", rank);fflush(stdout);

	free(msg);

	//MPI-work are stoped
	//noone mpi_procedure cannot be initialized after
	MPI_Finalize();

	return 0;
}

int velocityDefine(double unionvector[3], double vectorvalue, double *vvector, char *msg[]) {
	int i;
	double precision = 0.000000001;
	// Àáñîëþòíîå çíà÷åíèå åäèíè÷íîãî âåêòîðà = 1
	double valueunionV = 0;
	for (i = 0; i < 3; i++) {
		valueunionV += unionvector[i]*unionvector[i];
	}
	if ((fabs(valueunionV - 1) > precision)) {
		*msg = "velocityDefine (Errorcode_0BA): Uncorrect value of parameters. Union velocity vector";
		return 1;
	}

	// Àáñîëþòíîå çíà÷åíèå âåêòîðà ñêîðîñòè áîëüøå 0
	if (vectorvalue <= 0) {
		*msg = "velocityDefine (Errorcode_0BB): Uncorrect value of parameters. Value of velocity vector";
		return 1;
	}

	// çàäàíèå çíà÷åíèé êîìïîíåíò ñêîðîñòè
	for (i = 0; i < 3; i++) {
		vvector[i] = unionvector[i]*vectorvalue;
	}
	*msg = "velocitydefine (errorcode_0BC): program finished ok";
	return 0;
}

