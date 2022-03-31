#include "Boundaries.h"
#include "DoubleExponenta.h"
#include "Function.h"
#include "InputData.h"
//#include "utilprog.h"
#include <stdio.h>
#include <stdlib.h> //malloc
//#include <conio.h>
#include <string.h>	//memcpy
#include <mpi.h>
#include <math.h>
#include <time.h> //
#include "integration_procedure.h"
// точность вычислений


void Integration(double (*function3D)(double, double, double, double*), 
						int series, double *parameter, int myid, int size, int count, 
						MyStr *listINcom, double *subSumm, int *step, double *bounds) { 
	int KEY_exit = 1;	// key of exit from integration loop
	int incount= 0;		//count of income tasks 

	////////////////////////////////////////////////////////////////////////
	/// Tags-key of messages between coprocessors 
	// tag 10 - for count of send(receive) additional tasks
	// tag 11 - for send(receive) additional tasks
	// tag 12 - for send(receive) request on additional tasks
	// tag 13 - for send(receive) announce about no task to perfome
	
	///////////////////////////////////////////////////////////////
	// MPI statuses of received messages
	MPI_Status status10s;
	MPI_Status status10r;
	MPI_Status status11r;
	MPI_Status status12r;
	MPI_Status* status13r; // 13r - no tasks, received (array from all co-processors)

	/////////////////////////////////////////////////////////////
	// MPI requests 
	MPI_Request req10r = MPI_REQUEST_NULL;
	MPI_Request req12s = MPI_REQUEST_NULL;
	MPI_Request req12r;
	MPI_Request *req13r;
	
	///////////////////////////////////////////////////////////
	// flags
	int flag10s;
	int flag10r = 0;
	int flag11r = 1;
	int flag12s = 1;
	int flag12r = 1;
	int *flag13r;

	int k, i, j;			// iterators

	int countsend;		//count of sended messages 
	int *flagexitIn;	//coprocessor leave(no leave) procedure
	int *IncreaseProcId;//list of coprocessors for each processor
	//int *LeaveProcId;	
	//int LeaveN = 0;
	//int LeaveJ = 0;
	//int bye = 1;
	int InN;
	int InNTemp;
	int Inj;
	int Inprev;
	int flagExitProc = 0;
	int signexitIncom;
	int *signallIncom;
	int *signallAsk;

	MyStr *list;
	int *INcountcompare;	// preliminary information on the available count of tasks 
	
							////////////////////////////////////////
	// Simpson integration variables
	double subBubble;
	double subCells;
	double stepIterationX, stepIterationY, stepIterationZ; // a dinamic step of integration 
	double func;	// a calculated function value
	double lowX;	// a low boundary of dimension X
	double lowY;	// a low boundary of dimension Y
	double lowZ;	// a low boundary of dimension Z
	double X, Y, Z;	// variables of integration
	double t, u;	// Simpson`s coefficient for 1D
	double r, s;	// rize range of Simpson`s coefficient to 2D
	double o, p;	// rize range of Simpson`s coefficient to 3D
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//    MPI_Datatype sndrcvdata, sndrcvcount, sndrcvmsg;
    MPI_Datatype sndrcvcount, sndrcvmsg;
	MPI_Datatype Strtype, arhtypes, oldtypes[2];
    MPI_Aint offsets[2], arhsets;
    int blklens[2], arhblklens;

	/* sender */
    //////////////////////////////////////////////////////////////////////
	/*//       Create a strust MyStr struct type */						//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	arhtypes = MPI_DOUBLE;												//
		/*Количество элементов в каждом блоке*/							//
	arhblklens = 6;														//
		/*Смещение каждого блока, измеряемые в байтах*/					//
	arhsets = 0;														//
		/*Собственно создание новой структуры*/							//
    MPI_Type_create_struct( 1, &arhblklens, &arhsets, &arhtypes, &Strtype);	//
		/*Реестрация новой структуры*/									//
	MPI_Type_commit( &Strtype );										//
	//////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	/*//	  Create a send-recv struct type like Int-Char*/				//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	oldtypes[0] = MPI_INT;									//
    oldtypes[1] = MPI_CHAR;												//
        /*Количество элементов в каждом блоке*/							//
	blklens[0] = 1;														//
	blklens[1] = 1;														//
		/*Смещение каждого блока, измеряемые в байтах*/					//
	offsets[0] = 0;														//
    MPI_Type_size(MPI_INT, &offsets[1]);						//
		/*Собственно создание новой структуры*/							//
    MPI_Type_create_struct( 2, blklens, offsets, oldtypes, &sndrcvmsg );		//
		/*Реестрация новой структуры*/									//
	MPI_Type_commit( &sndrcvmsg );									//
	//////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////
	/*//	  Create a send-recv struct type like Int-Int*/				//
	//////////////////////////////////////////////////////////////////////
       /*Data types*/													//
	oldtypes[0] = MPI_INT;									//
    oldtypes[1] = MPI_INT;												//
        /*Количество элементов в каждом блоке*/							//
	blklens[0] = 1;														//
	blklens[1] = 1;														//
		/*Смещение каждого блока, измеряемые в байтах*/					//
	offsets[0] = 0;														//
    MPI_Type_size(MPI_INT, &offsets[1]);						//
		/*Собственно создание новой структуры*/							//
    MPI_Type_create_struct( 2, blklens, offsets, oldtypes, &sndrcvcount );		//
		/*Реестрация новой структуры*/									//
	MPI_Type_commit( &sndrcvcount );									//
	//////////////////////////////////////////////////////////////////////
	
	///сделать ограничения на переполнение(е. 10^5 -> 10^6);
	//////////////////////////////////////
	// ..Создание динамических массивов...
	// Выделение памяти на:
	// Стек заданий
	list = (MyStr *)malloc(100000*sizeof(MyStr));
	// очереди на сообщения oб отсутствии дополнительных заданий;;
	req13r = (MPI_Request*)malloc(size*sizeof(MPI_Request));
	// статусы полученых сообщений oб отсутствии дополнительных заданий;
	status13r = (MPI_Status*)malloc(size*sizeof(MPI_Status));
	// Список флагов о прикращении работы соседних процессов
	flagexitIn = (int*)malloc(size*sizeof(int));
	// Список флагов о входящих сигналах oб отсутствии дополнительных заданий;
	flag13r = (int*)malloc(size*sizeof(int));
	// Список доступных в прошлом количеств задач у соседних процессов;
	INcountcompare = (int*)malloc(size*sizeof(int));
	// Список соседних процессов для каждого процесса;
	IncreaseProcId = (int*)malloc(size*sizeof(int));
	// Список рабочих (или завершивших работу) процессов;
	//LeaveProcId = (int*)malloc(size*sizeof(int));
	// Список флагов на вхоящие задания;
	signallIncom = (int *)malloc(size*sizeof(int));
	// Список флагов запросов на дополните задания;
	signallAsk = (int *)malloc(size*sizeof(int));

	// копируем входящий список заданий
	memcpy(&list[1], listINcom, count*sizeof(MyStr));

	// начальная подсумма каждого процесса
	*subSumm = 0;

	// Изначально все соседние процессы равны
	// Предполагается, что у них нет ни одного задания
	for (k = 0; k < size; k++) {
		INcountcompare[k] = 0;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	/// example: processors 4, size of processors set is 6
	///			5(1), 0(2), 1(3), 2(4), 3(5),
	///         start indexing of co-processors from 1, 
	for(i = myid; i < size; i++) {
		IncreaseProcId[i-myid] = i;
	}
	for(i = 0; i < myid; i++) {
		IncreaseProcId[size-myid+i] = i;
	}
	InN = size-1;	// count of co-processors 
	InNTemp = InN;	// 
	Inj = 1;		// serial number of a source co-processor ,
	Inprev = 0;		// serial number of a previous source co-processor

	// assign start value 
	for (k = 1; k < size; k++) {																//
		signallIncom[k] = 1;																		//
	}

	// assign start value :
	// 1) flag-key: have co-processors left procedure?
	// 2) flag-key: no additional tasks
	// 3) request on co-processors procedure leaving 
	for (k = 1; k < size; k++) {																//
		flagexitIn[k] = 0;
		flag13r[k] = 0;
		req13r[k] = MPI_REQUEST_NULL;
	}
	
	//
	for (k = 1; k < size; k++) {																//
		signallAsk[k] = 1;																		//
	}																							//

	KEY_exit = 1;			// Key of exit

	// interception of artifacts
	for(i=1; i <= InN; i++) {
		MPI_Iprobe(IncreaseProcId[i], 10, MPI_COMM_WORLD, &flag10s, &status10s);
		if(flag10s) {
			int intbuf[2];
			MPI_Recv(intbuf, 1, sndrcvcount, IncreaseProcId[i], 10, MPI_COMM_WORLD, &status10s);			
		}
		MPI_Iprobe(IncreaseProcId[i], 10, MPI_COMM_WORLD, &flag10s, &status10s);
		if(flag10s) {
			int intbuf[2];
			MPI_Recv(intbuf, 1, sndrcvcount, IncreaseProcId[i], 10, MPI_COMM_WORLD, &status10s);			
		}
		MPI_Iprobe(IncreaseProcId[i], 13, MPI_COMM_WORLD, &flag13r[i], &status10s);
		if (flag13r[i]) {
			int byebuf = 0;
			MPI_Recv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
			flag13r[i] = 0;
		}
		if (flag13r[i]) {
			int byebuf = 0;
			MPI_Recv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
			flag13r[i] = 0;
		}
		if (flag13r[i]) {
			int byebuf = 0;
			MPI_Recv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
			flag13r[i] = 0;
		}
		if (flag13r[i]) {
			int byebuf = 0;
			MPI_Recv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
			flag13r[i] = 0;
		}
	}

	//	fprintf(stdout, " myid %d list[0].lowX = %f, list[0].upX = %f \n", myid, listINcom[0].lowX, listINcom[0].upX); fflush(stdout);
	//	fprintf(stdout, " myid %d list[0].upX = %f, list[0].upX = %f \n", myid, UPX, UPY); fflush(stdout);

	while (KEY_exit > 0) {
		int seriestmp12[1];
		while(count) {
			// count of iterations
			*step+=1;

			// elemental current interval of integration
			stepIterationX = fabs(list[count].upX - list[count].lowX) / 4.0;
			stepIterationY = fabs(list[count].upY - list[count].lowY) / 4.0;
			stepIterationZ = fabs(list[count].upZ - list[count].lowZ) / 4.0;

			//intermediate summs
			subBubble = 0.0; subCells = 0.0;

			//double o, p;	// Simpson`s coefficient for 1D
			for (i=0; i<=4; i++) {
				Z = list[count].lowZ + i*stepIterationZ; //Z
				if ((i == 0)||(i == 4)) {
					o=1;
					p=1;
				}
				if ((i == 2)) {
					o=2;
					p=4;
				}
				if ((i == 1)||(i == 3)) {
					o=4;
					p=0;
				}

				//	double r, s;	// rize range of Simpson`s coefficient to 2D
				for (j=0; j<=4; j++) {
					Y = list[count].lowY + j*stepIterationY ; //Y
					if ((j == 0)||(j==4)) {
						r=1;
						s=1;
					}
					if ((j == 2)) {
						r=2;
						s=4;
					}
					if ((j == 1)||(j==3)) {
						r=4;
						s=0;
					}


					//	double t, u;	// rize range of Simpson`s coefficient to 3D
					for (k=0; k<=4; k++) {
						X = list[count].lowX + k*stepIterationX ;	//X
						if ((k == 0)||(k==4)) {
							t=1;
							u=1;
						}
						if ((k == 2)) {
							t=2;
							u=4;
						}
						if ((k == 1)||(k==3)) {
							t=4;
							u=0;
						}
//						func = function(X, Y, Z, parameter);						
						func = DEtransformFunction3D(function3D, X, Y, Z, parameter, bounds);						
						subBubble += func*u*s*p;
						subCells += func*t*r*o;
					}
				}
			}
			subBubble = subBubble*stepIterationX*stepIterationY*stepIterationZ*8./27.;
			subCells = subCells*stepIterationX*stepIterationY*stepIterationZ/27.;//			subCells = subCells*stepIterationX/3;
		// УСЛОВИЕ НЕОБХОДИМОЙ ТОЧНОСТИ ВЫЧИСЛЕНИЙ
			//if (fabs(subCells-subBubble) > PRECISION*fabs(subCells)) { 	
			// subCells --- это значение функции после	DE-преобразования
			// Поскольку, DE-преобразование сильно подавляет значение функции,
			// то получаем, что на некоторых участках, далеких от вершини купола
			// DE-преобразования, значение subCells -> 0
			// Как результат, получаемая точность превышает нужную
			// Вывод: условие относительной ошибки сильно замедляет работу проограммы

			if (fabs(subCells-subBubble) > PRECISION) { // don't use '>='
				lowX = list[count].lowX;
				lowY = list[count].lowY;
				lowZ = list[count].lowZ;
				for (i = 0; i < 4; i++) {
					for (j = 0; j < 4; j++) {				
						for (k = 0; k < 4; k++) {
							list[count].lowX = lowX + (stepIterationX)*(k);
							list[count].upX = lowX + (stepIterationX)*(k+1);
							list[count].lowY = lowY + stepIterationY*(j);
							list[count].upY = lowY + (stepIterationY)*(j+1);
							list[count].lowZ = lowZ + stepIterationZ*(i);
							list[count].upZ = lowZ + (stepIterationZ)*(i+1);
							count++;
						}
					}
				}
				count--;
			} else 
			{
				*subSumm += subCells;
				count -= 1;
			}
			//if (0) {
			if (count >= 50) {
				if (flag12r) {
					MPI_Irecv(seriestmp12, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &req12r);
				}
				MPI_Test(&req12r, &flag12r, &status12r);
				if (flag12r) {
					//fprintf(stdout, "proc %d recv from proc %d ask for tasks with series %d\n", myid, status12r.MPI_SOURCE, series);
					//fflush(stdout);
					if (seriestmp12[0] == series) {
						MPI_Datatype sndrcvdata;
						double *bufdatasend;
						int position = 0;
						int intbuf[2];
						int mpiintsize, strtypesize;
						int seriesBuf[1];
						MyStr *sendlist;

						countsend = 0;

						if (count%2) { // непарное количество заданий
							countsend = (count-1)/2;
							count = countsend+1;
							sendlist = (MyStr*)malloc(countsend*sizeof(MyStr));
							memcpy(sendlist, &list[count+1], countsend*sizeof(MyStr));
						}
						else { // парное количество заданий
							countsend = (count)/2;																	////
							count = countsend;																		////
							sendlist = (MyStr*)malloc(countsend*sizeof(MyStr));
							memcpy(sendlist, &list[count+1], countsend*sizeof(MyStr));								////
						}
						if (countsend == 0) {
							int err = 001;
							fprintf(stdout, "ERROR_%d\n", err);fflush(stdout);
							MPI_Abort(MPI_COMM_WORLD, err);
						}

							//отсылка количества передаваемых данных с серийником
						intbuf[0] = series;			// серийник выполняемого общего задания
						intbuf[1] = countsend;		// количество частных заданий
						MPI_Send(intbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);

						//////////////////////////////////////////////////////////////////////
						/*//	  Create a send-recv struct type like MyStr*/				//
						//////////////////////////////////////////////////////////////////////
						   /*Data types*/													//
						oldtypes[0] = MPI_INT;									//
						oldtypes[1] = Strtype;												//
							/*Количество элементов в каждом блоке*/							//
						blklens[0] = 1;								//
						blklens[1] = countsend;												//
							/*Смещение каждого блока, измеряемые в байтах*/					//
						offsets[0] = 0;												//
						MPI_Type_size(MPI_INT, &offsets[1]);						//
							/*Собственно создание новой структуры*/							//
						MPI_Type_create_struct( 2, blklens, offsets, oldtypes, &sndrcvdata );		//
							/*Реестрация новой структуры*/									//
						MPI_Type_commit( &sndrcvdata );										//
						//////////////////////////////////////////////////////////////////////
						MPI_Type_size(MPI_INT, &mpiintsize);
						MPI_Type_size(Strtype, &strtypesize);
						bufdatasend = (double*)malloc(countsend*strtypesize+mpiintsize);
						//fprintf(stdout, "buffer size is %d - bufdatasend size is %d\n", countsend* strtypesize + mpiintsize,sizeof(bufdatasend)); fflush(stdout);
						position = 0;

						//fprintf(stdout, "11S>proc %d START pack COUNT tasks to proc %d %d(%d) tasks\n", myid, status12r.MPI_SOURCE, countsend, intbuf[1]);
						//fflush(stdout);


						//MPI_Pack(&intbuf[0], 1, MPI_INT, bufdatasend, countsend*strtypesize+mpiintsize, &position, MPI_COMM_WORLD);

						//fprintf(stdout, "10S>proc %d pack START TASKS to proc %d %d(%d) tasks\n", myid, status12r.MPI_SOURCE, countsend, intbuf[1]);
						//fflush(stdout);
						seriesBuf[0] = series;
						MPI_Pack(seriesBuf, 1, MPI_INT, bufdatasend, countsend * strtypesize + mpiintsize, &position, MPI_COMM_WORLD);
						MPI_Pack(sendlist, countsend, Strtype, bufdatasend, countsend*strtypesize+mpiintsize, &position, MPI_COMM_WORLD);

						//fprintf(stdout, "10S>proc %d pack FINISH TASKS to proc %d %d(%d) tasks. Size of data is %d\n", 
						//	myid, status12r.MPI_SOURCE, countsend, intbuf[1], sizeof(countsend*strtypesize+mpiintsize));
						//fflush(stdout);

						MPI_Send(bufdatasend, position, MPI_PACKED, status12r.MPI_SOURCE, 11, MPI_COMM_WORLD);


						//MPI_Send(bufdatasend, 1, sndrcvdata, status12r.MPI_SOURCE, 11, MPI_COMM_WORLD);
						//MPI_Send(sendlist, 6*countsend, Strtype, status12r.MPI_SOURCE, 11, MPI_COMM_WORLD);
						//fprintf(stdout, "proc %d HAVE already SEND to proc %d %d(%d) tasks\n", myid, status12r.MPI_SOURCE, countsend, intbuf[1]);
						//fflush(stdout);
						free(bufdatasend);
						free(sendlist);
						MPI_Type_free(&sndrcvdata);
					}
				}
			}
		}
		///////////////////////////////////////////////////////////////////
		//	if(count == 0) 	// е. у процесса закончились задания	///
		///////////////////////////////////////////////////////////////////
	//break;
		//fprintf(stdout, "proc %d: no tasks. count of coprocessors is %d\n", myid, InNTemp); fflush(stdout);
	while(!count) {
		//KEY_exit = 0;
		//break;
		int inbuf[2];
		int outbuf[2]; 
		//int byebuf=1;
		int byebuf[1];
		int seriestmp10;
///
///получить сообщения о прекращении работы других процессов
///
		for(i = 1; i <= InNTemp; i++) {
			if(!flagexitIn[i]){ //flag on request: is coprocessors leave procedure?
				// nonblock request on receiving leaving message
				//MPI_Irecv(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
				MPI_Irecv(byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD, &req13r[i]);
				// no requests more
				flagexitIn[i] = 1;
				//fprintf(stdout, "13>proc %d check leaved co-processors %d \n", myid, IncreaseProcId[i]); fflush(stdout);

			}	
			if(!flag13r[i]) { //TEST: Have leaving message already been received?
				MPI_Test(&req13r[i], &flag13r[i], &status13r[i]);
			}

			if (flag13r[i]) { // Leaving message have already been received
				//test is series correct (Is this message no artifact?);
				int seriestmp13;
				seriestmp13 = byebuf[0];
				if (seriestmp13 == series) {
					//fprintf(stdout, "proc %d know proc %d has left program\n", myid, IncreaseProcId[i]);
					//fflush(stdout);
					//decrease count of co-processors
					InNTemp --;
					//fprintf(stdout, "13R>proc %d knows that proc %d (%d) LEAVE procedure\n", myid, IncreaseProcId[i], InNTemp); fflush(stdout);
					// cancel request on count of additional tasks
					if (Inj == i) {
						flag11r = 1;
						if (!flag10r) {
							MPI_Cancel(&req10r);
							flag10r = 1;
						}
					}
					// remove co-processors from list 
					if(Inj > i) {
						Inj --;
					}
					for(j = i; j <= InNTemp; j++){
						IncreaseProcId[j] = IncreaseProcId[j+1];
						signallIncom[j] = signallIncom[j+1]; 
					}
					for (j = i; j<=InNTemp; j++) {
						req13r[j] = req13r[j+1];
						flag13r[j] = flag13r[j+1] ;
					}
				} else { //series is  uncorrect: This message is artifact.;
					// reassign flag value (ready for monitoring on leaving message)
					//fprintf(stdout, "13RW>proc %d recieve from proc %d msg with WRONG SERIES %d (%d) \n", myid, IncreaseProcId[i], seriestmp13,series); fflush(stdout);
					flagexitIn[i] = 0;
					flag13r[i]=0;
				}
			}
        }

		InN = InNTemp; // reassign count of working co-processors;

		//It's no working co-processors more
		if (InN == 0) {
			if (!flag12r) {
				//cancell request on additional tasks
				MPI_Cancel(&req12r);
			}
			if (!flag10r) {
				//cancell request on count of  additional tasks
		        MPI_Cancel(&req10r);
			}
			//fprintf(stdout, "proc %d leaving caused by no more processes\n", myid);	fflush(stdout);
			KEY_exit  = 0;
			break;
		}
		
		// ASK additional tasks
		if (flag11r) {
			int tmp[1];
			tmp[0] = series;
			MPI_Send(tmp, 1, MPI_INT, IncreaseProcId[Inj], 12, MPI_COMM_WORLD);
		}

		// TEST: Do we have an ask on additional rasks
        if (flag12r) {
            MPI_Irecv(seriestmp12, 1, MPI_INT, MPI_ANY_SOURCE, 12, MPI_COMM_WORLD, &req12r);										////
            flag12r = 0;																									////
        }
        MPI_Test(&req12r, &flag12r, &status12r);
 
		if (flag12r) { //ask have been received
			//fprintf(stdout, "proc %d has %d series and recv %d tmpser\n", myid, series, series);
			//fflush(stdout);
			if( seriestmp12[0] == series) {
				outbuf[0] = series;
				outbuf[1] = 0;
				//send msg: 'no task'
				MPI_Send(outbuf, 2, MPI_INT, status12r.MPI_SOURCE, 10, MPI_COMM_WORLD);
				//fprintf(stdout, "10S>proc %d send I'M EMPY to proc %d\n", myid, status12r.MPI_SOURCE);
				//fflush(stdout);
				countsend = 0;
				flag12r = 1;
			} 
        }

		//TEST: Have additional messages been send?
		if (flag11r) {
            MPI_Irecv(inbuf, 2, MPI_INT, IncreaseProcId[Inj], 10, MPI_COMM_WORLD, &req10r);
            flag11r = 0;
            flag10r = 0;
        }
        MPI_Test(&req10r, &flag10r, &status10r);

		//additional tasks is send
		if (flag10r) {
			flag11r = 1;    
			seriestmp10 = inbuf[0];
			if (seriestmp10 == series) 
			{
				incount = inbuf[1];
				//fprintf(stdout, "proc %d recv from proc %d %d tasks\n", myid, status10r.MPI_SOURCE, incount);
				//fflush(stdout);
				//test is set of task empty?
			   if (incount == 0) { 
				   //fprintf(stdout, "10R>proc %d recv from proc %d msg: NO TASK \n", myid, status10r.MPI_SOURCE, incount);
				   //fflush(stdout);
					i = 0;
					signallIncom[Inj] = 0; // sender does not have tasks
					//ready to ask next co-processor
					if ((Inj <= (InN)) && (Inj > 0)) {
						if (incount >= INcountcompare[Inprev] ) {
							Inprev = Inj;
							if (Inj == (InN)) {
								Inj = 1;
							} else {
								Inj++;
							}
						} else {
							Inprev = Inj;
							if (Inj == 1) {
								Inj = (InN);
							} else {
								Inj--;
							}
						}
					}
					//additional task is presented
			   } else if (incount > 0) {
					MPI_Datatype sndrcvdata;
					int seriestmp11[1];
					int position = 0;
					double *bufdatarecv;
					int bufferSize;
					//////////////////////////////////////////////////////////////////////
					/*//	  Create a send-recv struct type like MyStr*/				//
					//////////////////////////////////////////////////////////////////////
					   /*Data types*/													//
					oldtypes[0] = MPI_INT;									//
					oldtypes[1] = Strtype;												//
						/*Количество элементов в каждом блоке*/							//
					blklens[0] = 1;														//
					blklens[1] = incount;												//
						/*Смещение каждого блока, измеряемые в байтах*/					//
					offsets[0] = 0;														//
					MPI_Type_size(MPI_INT, &offsets[1]);						//
						/*Собственно создание новой структуры*/							//
					MPI_Type_create_struct( 2, blklens, offsets, oldtypes, &sndrcvdata );		//
						/*Реестрация новой структуры*/									//
					MPI_Type_commit( &sndrcvdata );										//
					MPI_Type_size(sndrcvdata, &bufferSize);
					//////////////////////////////////////////////////////////////////////
					bufdatarecv = (double*)malloc(incount*sizeof(MyStr)+sizeof(int));
					
					//bufdatarecv = (double*)malloc(incount * sizeof(MyStr) + sizeof(double));
					//fprintf(stdout, "10R>proc %d recv msg from %d: I send to you set with %d(%d) tasks. Size of Message is %d\n", 
					//	myid, status10r.MPI_SOURCE, incount, inbuf[1], sizeof(bufdatarecv));
					//fflush(stdout);

					MPI_Recv(bufdatarecv, bufferSize, MPI_PACKED, IncreaseProcId[Inj], 11, MPI_COMM_WORLD, &status11r);
					//MPI_Recv(bufdatarecv, 1, sndrcvdata , IncreaseProcId[Inj], 11, MPI_COMM_WORLD, &status11r);

					//fprintf(stdout, "10R> proc %d START UNPACK message from %d\n", myid, status10r.MPI_SOURCE);
					//fflush(stdout);

					MPI_Unpack(bufdatarecv, incount*sizeof(MyStr)+ offsets[1], &position, seriestmp11, 1, MPI_INT, MPI_COMM_WORLD);
					if (seriestmp11[0] == series) {
						MyStr *recvlist;
						recvlist =(MyStr*)malloc(incount*sizeof(MyStr));
						MPI_Unpack(bufdatarecv, incount*sizeof(MyStr)+sizeof(int), &position, recvlist, incount, Strtype, MPI_COMM_WORLD);   

						memcpy(&list[1], recvlist, incount*sizeof(MyStr));										//
						signallIncom[Inj] = 1;
						count = incount;
						//change the order of co-processors list	
						if ((Inj <= (InN)) && (Inj > 0)) {
							if (incount >= INcountcompare[Inprev] ) {
								Inprev = Inj;
								if (Inj == (InN)) {
									Inj = 1;
								} else {
									Inj++;
								}
							} else {
								Inprev = Inj;
								if (Inj == 1) {
									Inj = (InN);
								} else {
									Inj--;
								}
							}
						}
						free(recvlist);
					}
					//fprintf(stdout, "10R> proc %d FINISH UNPACK message from %d\n", myid, status10r.MPI_SOURCE);
					//fflush(stdout);

					MPI_Type_free(&sndrcvdata);
					free(bufdatarecv);
					break;
				} else if (incount < 0) {
					//fprintf(stdout, "\nproc %d !!!!!SystemError!!!!\n", myid);
					//fflush(stdout);
					MPI_Finalize();
				}
				INcountcompare[Inprev] = incount;	
			} 
        }

		// If all co-processors don't have tasks
		// It's a enough condition to leave procedure

		// send message: 'I leave procedure'

		if (!flagExitProc) {
            signexitIncom = 0;
			for (k = 1; k <= InN; k++) {		//
                signexitIncom += signallIncom[k];			//
            }
            if (!signexitIncom) {
                for (i = 1; i <= InNTemp; i++) {
					//byebuf = series;
					byebuf[0] = series;
					//MPI_Send(&byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD);
					MPI_Send(byebuf, 1, MPI_INT, IncreaseProcId[i], 13, MPI_COMM_WORLD);
					//fprintf(stdout, "Proc %d send to proc %d msg with series %d: \"I leave you\"\n", myid, IncreaseProcId[i], byebuf[0]);
					MPI_Cancel(&req13r[i]);
				}
				//fprintf(stdout, "proc %d leaving caused by no more tasks\n", myid);fflush(stdout);

				if (!flag12r) {
                    MPI_Cancel(&req12r);
				}
				if (!flag10r) {
                    MPI_Cancel(&req10r);
				}
                KEY_exit = 0; // сигнал выхода
                break;
            }
		}
	}


	}
	free(status13r);
	free(flagexitIn);
	free(req13r);
	free(flag13r);
	free(INcountcompare);
	free(IncreaseProcId);
	free(list);
	free(signallIncom);
	free(signallAsk);

	fprintf(stdout, "*\n");
	fflush(stdout);
}
