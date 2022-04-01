// Интервалы интегрирования разбиваются на части,	//
// соответственно к количеству рабочих процессов	//
// После вызывает интегрирующий алгоритм,			//
// передавая ему как параметр масив интервалов

/* define boundaries*/
#include "Boundaries.h"  
/* Double Exponenta*/
/* L.Ye NUMERICAL QUADRATURE: THEORY AND COMPUTATION*/
#include "DoubleExponenta.h" 
/*Function of integration*/
#include "Function.h"

#include "InputData.h"
#include "integration_procedure.h"
#include <stdio.h>
#include <stdlib.h> //malloc
//#include <conio.h>
#include <string.h>	//memcpy
#include <mpi.h>
#include <math.h>
#include <time.h> //

void PreIntegration (int rank, int size, int series, double *set_parmtrs, double *result) {
	int MPIErrorCode=1;

	double subSumm=0.0, totalSumm = 0.0;
	double intervalX, intervalY, intervalZ;
	double statisticTime = 0.0, totalstatisticTime = 0.0;
	
	double bounds[6];		

	MyStr *list;			
	MyStr *sendlist;	

	int *step;
	int tempstep = 0;
	int totalstep = 0;
	int i, j, k;
	double lowX, lowY, lowZ;
	int count = 0;
	double startwtime = 0.0;
	double endwtime = 0.0;

	MPI_Status status;
	MPI_Request reqSumm[1] = {MPI_REQUEST_NULL};

	////////////////////////////
	//..Ñîçäàåì íîâûé òèï MPI...
	// Íàçâàíèå òèïà
	MPI_Datatype strtype;
	// Äàííûå ñ èäåíòèôèêàöèîííûì íîìåðîì id
	MPI_Datatype IdData;
	// ñîäåðæàíèå	
	MPI_Datatype type[6] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
	MPI_Datatype IdType[7] = {MPI_UNSIGNED_LONG, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

	// êîëè÷åñòâî ýëåìåíòîâ
	int blocklen[6] = {1, 1, 1, 1, 1, 1};
	int IdBlocklen[7] = {1, 1, 1, 1, 1, 1, 1};
	// ðàññòîÿíèå ìåæäó ýëåìåíòàìè (â ïàìÿòè)
	MPI_Aint disp[6] = {0, 0, 0, 0, 0, 0};
	MPI_Aint IdDisp[7] = {0, 0, 0, 0, 0, 0, 0};
	// êîíñòðóèðóåì òèï
	MPI_Type_create_struct(6, blocklen, disp, type, &strtype);
	MPI_Type_create_struct(7, IdBlocklen, IdDisp, IdType, &IdData);
	// ðååñòðàöèÿ òèïà
	MPI_Type_commit(&strtype);
	MPI_Type_commit(&IdData);
	
	// çàñåêàåì âðåìÿ ðàáîòû ïðîãðàììû
	startwtime = MPI_Wtime();	

	list = (MyStr *)malloc(size*size*size*sizeof(MyStr));//ñïèñîê äëÿ èíòåðâàëîâ èíòåãðèðîâàíèÿ
	sendlist = (MyStr *)malloc(size*size*sizeof(MyStr)); // ñïèñîê äëÿ îòñûëêè çàäàíèé (èíòåðâàëîâ èíòåãðèðîâàíèÿ)
	
	step = (int *)malloc(size*sizeof(int));			// ìàññèâ ñ÷åò÷èêîâ
														//, óêàçûâàþùèõ íà êîëè÷åñòâî èòåðàöèé, 
														// ïðîâåäåííûõ êàæäûì ïðîöåññîì
	for (i = 0; i <size; i++) {
		step[i] = 0;
	}
	// in ELI_La parameter[0-2] is a Direction of unit velocity vector 
	//parameter[0] =  sqrt(1.0 - wz*wz);
	//parameter[1] =  wz;

	// in ELI_La parameter[3-4] is a Temperature of electron gas, Dimentionless
	//parameter[2] =  TAU/3.;
	//parameter[3] =  TAU/3.;

	count = 0;
	if (rank == 0) {
		// ïðåîáðàçîâàíèå êîîðäèíàò ïî àëãîðèòìó äâîéíîé ýêñïîíåíòû
		list[count].lowX = - DEboundary();   
		list[count].upX = -list[count].lowX; 
		list[count].lowY = list[count].lowX;
		list[count].upY = -list[count].lowX;
		list[count].lowZ = list[count].lowX;
		list[count].upZ = -list[count].lowX;

		//fprintf(stdout, "%g\t%g\t%g\t%g\n", list[0].upZ, list[0].lowZ, list[0].upY, list[0].lowY);fflush(stdout);
		
		lowX = list[count].lowX ;
		lowY = list[count].lowY ;
		lowZ = list[count].lowZ ;

		// Ðàçáèåíèå èíòåðâàëà èíòåãðèðîâàíèÿ íà áîëåå ìåëêèå ÷àñòè
		// øàãè ðàçáèåíèÿ
		intervalX = fabs(list[count].upX - list[count].lowX) / size ;
		intervalY = fabs(list[count].upY - list[count].lowY) / size ;
		intervalZ = fabs(list[count].upZ - list[count].lowZ) / size ;
		for (j = 0; j < size; j++) {
			for (k = 0; k < size; k++) {
				for (i = 0; i < size; i++) {
					list[count].lowX = lowX + intervalX*k;
					list[count].upX = lowX + (intervalX)*(k+1);
					list[count].lowY = lowY + intervalY*j;
					list[count].upY = lowY + (intervalY)*(j+1);
					list[count].lowZ = lowZ + intervalZ*i;
					list[count].upZ = lowZ + (intervalZ)*(i+1);
					count++;
				}
			}
		}
		count = count/size; // êîëè÷åñòî ýëåìåíòîâ äëÿ îòïðàâêè
		// #WARNING
		// Ðàññûëêó ìîæíî çàìåíèòü ïðîñòî âûáîðî÷íûì öûêëîì,
		// êîãäà èç öûêëà êàæäûé ïðîöåññ âûáèðàåò ýëëåìåíò 
		// ÷åðåç îïðåäåëåííûé øàã
		for (i = 1; i < size; i++) {
			// Îòñûëàåì êàæäîìó ïðîöåññó çàäàíèÿ:
			// îòñûëêà êîëè÷åñòâà çàäàíèé
			MPI_Send(&count, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			// îòñûëêà çàäàíèé count øò.
			memcpy(sendlist, &list[i*(count)], count*sizeof(MyStr));
			MPI_Send(sendlist, 6*count, strtype, i, 0, MPI_COMM_WORLD);
		}
	} else {
		// Ïîëó÷åíèå çàäàíèé äëÿ ðàáîòû:
		// ïîëó÷åíèå êîëè÷åñòâà ïðèñûëàåìûõ çàäàíèé
		MPI_Recv(&count, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		// ïîëó÷åíèå çàäàíèé count øò.
		//MPI_Recv(list, 4*count, strtype, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(list, 6*count, strtype, 0, 0, MPI_COMM_WORLD, &status);
		// #WARNING 
		// Ïîëó÷þùèé ñïèñîê ìîæåò áûòü è ìåíüøèì, 
		// ðàçìåðà sendlist
	}
	bounds[0] = LOWX;
	bounds[1] = UPX;
	bounds[2] = LOWY;
	bounds[3] = UPY;
	bounds[4] = LOWZ;
	bounds[5] = UPZ;

	//if (rank == 0) {
	//	for(i = 0; i < 6; i++) {
	//	fprintf(stdout, "bound %d = %g\n", i, bounds[i]); fflush(stdout);
	//	}
	//}
	fprintf(stdout, "Proc %d(%d) start!\n", rank, size); fflush(stdout);
	///%%%%%%%%%%%%%%%%%%%%
	fprintf(stdout, "tau_per_norm = %f (%f eV), tau_par_norm = %f (%f eV) \n", 
		tau_perp_norm(set_parmtrs[4], set_parmtrs[5]), set_parmtrs[4],
		tau_parallel_norm(set_parmtrs[4], set_parmtrs[5]), set_parmtrs[5]); fflush(stdout);
	/////////////////////////////////////////////////////////////
	//Integration(ELiMA, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);
	//fprintf(stdout, "Proc %d(%d) start!\n", rank, size);
//	if (rank == 0) {
		//fprintf(stdout, "proc %d(%d): series %d, count %d, \n par[2] = %f, bounds[0] = %f\n", rank, size, series, count, set_parmtrs[2], bounds[0]);
//	}
	//Integration(ELIA, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);
	///%%%%%%%%%%%%%%%%%%%%	
	//Integration(ELI, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);
	//Integration(function3D, series, set_parmtrs, rank, size, count, list, &subSumm, &tempstep, bounds);

	///////////////////////////////////////////////////////////
	//Ïðîñóììèðóåì ïðîìåæóòî÷íûå çíà÷åíèÿ ñóììû//
	// Ñáîð ðåçóëüòàòîâ îò êàæäîãî ïðîöåññà
	if (rank != 0) {
		MPI_Send(&tempstep, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	} else {
		totalstep += tempstep;
		for(i = 1; i < size; i++) {	
			MPI_Recv(&step[i], 1, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
			totalstep +=step[i];
		}
	}
	if (rank != 0) {
		MPI_Send(&subSumm, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
	} else {
		totalSumm += subSumm;
		for(i = 1; i < size; i++) {	
			MPI_Recv(&subSumm, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
			totalSumm += subSumm;
		}
	}
	if (rank == 0) {
		step[0] = tempstep;
		for(i = 0; i < size; i++) {
			fprintf(stdout, "proc %d do %6d (%5.2f %%) from %d iterations\n", 
					i, step[i], (double)step[i]/(double)totalstep*100., totalstep);
		}
		fprintf(stdout, "summ = %e\n", totalSumm);
		fflush(stdout);
		*result = totalSumm;
		endwtime = MPI_Wtime();

		printf("wall clock time = %f\n", endwtime-startwtime);fflush(stdout);
	}


	// îñâîáîæäåíèå ïàìÿòè
	free(step);
	free(list);
	free(sendlist);
	//printf("proc %d****************\n", rank);fflush(stdout);
}
