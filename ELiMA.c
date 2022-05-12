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
#include <gsl/gsl_sf_bessel.h> 

# ifndef PRECISION
#define PRECISION 0.0001
#endif

// numbers of task parameters
#define COUNTPARMTRS 11	// 3 - the union 3-dim velocity vector
			// 1 - an absolutly value of the velocity
			// 2 - the 2-dim electron Maxwell temperature (parallel and perpendicular component)
			// 1 - electron plasma frequency magnitude
			// 3 - the union 3-dim magnetic field vector
			// 1 - an absolutly value of the magnetic field
// Program name
#define PROG_NAME "ELiMA"
// Program version
#define PROG_VERSION ".v.git"
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
#define PARAMETERS_COUNT 11
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

	///*Timing of performing*/
	//struct tm * tmptime;

	/*MPI rutine initialization*/
	int MPIErrorCode = 1;

	/*MPI initialization*/
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
		fflush(stdout);
	}

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

	// Read parameters from file
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
			printf("par[%d]=%f\n", i, setParmtrsValue[i]);
		}
	}

	//share parameters value with coprocessors
	if (processorRank == 0)
	{
		for (int i = 1; i < processorSetSize; i++) {
			MPI_Send(setParmtrsValue, PARAMETERS_COUNT, MPI_DOUBLE, i, 123, comm);
		}
	}

	if (processorRank != 0)
	{
		MPI_Status status;
		MPI_Recv(setParmtrsValue, PARAMETERS_COUNT, MPI_DOUBLE, 0, 123, comm, &status);
		//for (int i = 0; i < PARAMETERS_COUNT; i++) {
		//	printf("par[%d]=%f\n", i, setParmtrsValue[i]);
		//}

	}

	// run integration procedure

	// summ subresults from coprocessors
		// send all subresults to root coprocessor
		// summ results by 

	// delete pointers
		// *headers, *parameters
	// 


	//////////////////////////////////////////////////////////////////
	// Call intagration function
	//////////////////////////////////////////////////////////////////
	double result = 0.;		// result of calculation

		PreIntegration(processorRank, processorSetSize, 1, setParmtrsValue, &result);
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
	//free(setParmtrsName);
	free(setParmtrsValue);
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

