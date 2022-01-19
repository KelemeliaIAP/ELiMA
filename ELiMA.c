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
#define FILEPATH "D:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\MPI_Results\\"

int velocityDefine(double unionvector[3], double vectorvalue, double *vvector, char *msg[]);

int main (int argc, char **argv) {
	// log-file
	FILE *LOGfile = NULL;

//	time_t timeW;
	struct tm * tmptime;
	double result = 0.;	// result of calculation
	int series = 1;		// index of program's start

			
	int MPIErrorCode = 1;
	int rank, size;
	int flagMPI_Init;

	double *setparmtrs;		// ìàññèâ ïîä ïàðàìåòðû çàäà÷è
	double unionvectorV[3];	// åäèíè÷íûé âåêòîð-ñêîðîñòè
	double vectorV[3];		// âåêòîð-ñêîðîñòè
	char ** msg;		// a message informs how the program works
	//int i;
	MPI_Comm MPI_COMM_External;
	double *set_parmtrs; //an array for the incoming parameters
	//to define a size of the array for the incoming parameters
	set_parmtrs = (double*)malloc(argc*sizeof(double));	
	msg =  (char**) malloc(100*sizeof(char));
	
	MPI_Init(&argc, &argv); // MPI initialization
	MPI_Initialized (&flagMPI_Init); // test MPI initialization
	
	//definition of number of work processors
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	//definition of indexes of work processors
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // index number of processor

	// Create new MPI communicator
	MPI_Comm_dup(MPI_COMM_WORLD, &MPI_COMM_External);
	//if (rank == 0) {
	//	//Define start time of program
	//	//Data is written to log.file in format
	//	//Day-Month-Year Hour(00-23):Minute(00-59):Second(00-61)
	//	time_t tt = time(NULL);
	//	struct tm ptm;
	//	char buf[26];
	//	asctime_s(buf, sizeof(buf), localtime_s(&tt, &buf));
	//	//printf("Current local time and date: %s", buf);
	////	for Windows
	//	if ( (fopen_s(&LOGfile ,"C:\\Program Files (x86)\\MPICH2\\bin\\log.dat", "w")) ) { // TRUE is equal "File was not opened"
	//		int MPIerrorCode = 1; // File was not opened
	//		fprintf(stdout, "Can't open log-file\n");fflush(stdout);
	//		MPI_Abort(MPI_COMM_WORLD, MPIerrorCode);
	//	}

	//	fprintf(LOGfile, "START\n");
	//	fprintf(LOGfile, "%s", asctime_s(buf, sizeof buf, localtime_s(&tt, &buf)));
	//	fflush(LOGfile);
	//}
	

	// Program is valid when more than one processors
	if (rank == 0) {
		if(size <= 1) {
			fprintf(stdout, "Number of processes have to be more than 1\n");fflush(stdout);
			fprintf(LOGfile, "Program was stoped. Number of processes have to be more than 1\n");fflush(LOGfile);
			MPI_Abort(MPI_COMM_WORLD, MPIErrorCode);
		}
	}

	////////////////////////////
	// Главная часть программы//
	////////////////////////////
	{
		unionvectorV[0] = atof(argv[2]); //nVx
		unionvectorV[1] = atof(argv[3]); //nVy
		unionvectorV[2] = atof(argv[4]); //nVz

		if(velocityDefine(unionvectorV, atof(argv[5]), vectorV, msg)) { // error
			fprintf(stdout,"%s", msg); 
			return 1;
		}

		set_parmtrs[0] = vectorV[0];	// Vx
		set_parmtrs[1] = vectorV[1];	// Vy
		set_parmtrs[2] = vectorV[2];	// Vz
		set_parmtrs[3] = atof(argv[5]);	// |V|
		set_parmtrs[4] = atof(argv[6]);	// Te_per
		set_parmtrs[5] = atof(argv[7]);	// Te_par

	// !!!!не забути підключити ці параметри до розрахунків у наступних файлах!!!!
		set_parmtrs[6] = delta(set_parmtrs[4]);	// Delta
		set_parmtrs[7] = tau_perpendicular();	
		set_parmtrs[8] = tau_parallel(atof(argv[7]), atof(argv[6])); //(Te_par, Te_per)
	
	}

		//////////////////////////////////////////////////////////////////
		// Call intagration function
		//////////////////////////////////////////////////////////////////
		PreIntegration(rank, size, series, set_parmtrs, &result);
		//////////////////////////////////////////////////////////////////


		/// Открытие файла для вывода результатов вычислений
		if (rank == 0) {
			double temperature_coeficient; 
			// создание имени файла, который будет 
			// содержать результаты вычислений
			char *file_name;			// массив под имя файла
			int length_name;			// длинна имени файла
			// Result of calculation are written to the file
			FILE *OUTResult;		// Выводится результат вычислений

			// V_e/V_o = sqrt(T_e * C^2/(m_eC^2)) = sqrt(T_e) * 41.9382450;
			// [V_e/V_o]/[T_e] = cm/eV/s
			temperature_coeficient = 41.94;
			// Длинна имени файла
//			length_name = sizeof(PROG_NAME)+sizeof(PROG_VERSION)+sizeof(FILE_FORMAT)
//							+sizeof(TEMPERATURAPer)+sizeof(TEMPERATURAPar)+sizeof(double)+sizeof(FILEPATH);
			length_name = sizeof(PROG_NAME)+sizeof(PROG_VERSION)+sizeof(FILE_FORMAT)
							+sizeof(TEMPERATURAPer)+sizeof(TEMPERATURAPar)+sizeof(double);
			// выделение памяти под массив, содержащий имя файла
			file_name = (char*)malloc(length_name);
			// образование имени файла
			
			fprintf(file_name, "%s%s%s%g%s%g%s", PROG_NAME, PROG_VERSION, 
							 TEMPERATURAPer, (double)set_parmtrs[4], TEMPERATURAPar, (double)set_parmtrs[5], FILE_FORMAT);
			
			fopen_s(&OUTResult, file_name, "a+");
			if (!OUTResult) {
				int MPIerrorCode = 1; // File was not opened
				fprintf(stdout, "ERROR: Can't open OUTResult-file\n");fflush(stdout);
				MPI_Abort(MPI_COMM_WORLD, MPIerrorCode);
			}
			fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
			fprintf(LOGfile, "ResultFile was opened successfully\n");fflush(stdout);

//			// Vi/V0	\t	Te \t	Delta \t	Tau \t 
//			// result \t precision \t program \t 
//			// bouX_ \t bouX^ \t bouK_ \t bouK^\n", 
			fprintf(OUTResult, "%g\t%g\t%g\t%g\t %g\t%g\t%g\t %g\t%g\t%g\t%s%s\t  %g\t%g\t%g\t%g\t%g\t%g\n", 
						set_parmtrs[0], set_parmtrs[1], set_parmtrs[2], set_parmtrs[3],
						set_parmtrs[3]*temperature_coeficient*sqrt(set_parmtrs[4]), set_parmtrs[4], set_parmtrs[5], 
						result, result/(temperature_coeficient*sqrt(set_parmtrs[4])), PRECISION, PROG_NAME, PROG_VERSION, 
						LOWX, UPX, LOWY, UPY, LOWZ, UPZ);fflush(OUTResult);
								
			fclose(OUTResult);
			if(!OUTResult) {
				fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
				fprintf(LOGfile, "File was not closed\n");		
				fprintf(stdout, "File was not closed\n");		
			} else {
				fprintf(LOGfile, "\\---------------------------------\\ \n");fflush(stdout);
				fprintf(LOGfile, "ResultFile was closed successfully\n");fflush(stdout);
			}

		}
	/// Закрытие файлов для вывода (log-file, result-file)
	if(rank == 0) {
		if(LOGfile) {
			time_t tt = time(NULL);
			struct tm ptm;
			char buf[26];
			asctime_s(buf, sizeof buf, localtime_s(&tt, &buf));
			printf("Current local time and date: %s", buf);
			
			fprintf(LOGfile, "%s", asctime_s(buf, sizeof buf, localtime_s(&tt, &buf)));
			fprintf(LOGfile, "END\n");		
			fclose(LOGfile);
		}
	}

	fprintf(stdout, "proc %d finalize\n", rank);fflush(stdout);
	//MPI-work are stoped//
	//noone mpi_procedure cannot be initialized after//
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

