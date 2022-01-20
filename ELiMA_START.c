#include <windows.h>
#include <stdio.h>
#include <tchar.h>
#include <conio.h>
#include <string.h>
#include <math.h>
#include<io.h>
// Errorcode:
	// первая цифра --- глубина программ (0 - стартовая программа)
	// 1-й Символ --- порядковій номер подпрограммі
	// 2-й Символ --- порядковій номер исключения

// Пути к файлам
//#define FILEPATH "\"D:\\01_Khelemelia\\02_IAP\\07_Programming\\01_MPI_C\\01_MPI_Projects\\01_MPI_ELI\\ELiMA.v.5.04\\Debug\\ELiMA.v.5.04.exe\""

//#define FILEPATH "\"D:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\ELiMA_Programming\\ELiMA\\x64\\Debug\\ELiMA.exe\""
#define FILEPATH "\"D:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\ELiMA_Programming\\ELiMA\\x64\\Debug\\ELiMA.exe\""
//#define FILEPATH "\"D:\\01_Khelemelia\\02_IAP\\07_Programming\\01_MPI_C\\01_MPI_Projects\\01_MPI_ELI\\ELiMA.v.5.05\\Debug\\ELiMA.v.5.05.exe\""
//#define FILEPATH "\"D:\\01_Khelemelia\\02_IAP\\07_Programming\\01_MPI_C\\01_MPI_Projects\\01_MPI_ELI\\ELIA.v.1.04\\Debug\\ELIA.v.1.04.exe\""
#define DATAFILEPATH "D:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\MPI_Results\\data.dat"

// функция записи параметров в файл.
// Принимает:
	// file_name --- название файла 
	// up, low --- крайние возможные значения параметров
	// step --- интервал между значениями параметров
	// *parameter --- массив параметров
	// msg --- сообщение о состоянии выполненной работы
// Возвращает:
	// 0: корректное завершение работы
	// 1: некорректное завершение работы
	// msg: сообщение о состоянии выполненной работы
int dataFileWrite(char *file_name, double low, double up, double step, double *parameter, int countparameter, char *msg[]);

// Определение компонент вектора скорости
// Принимает: 
	// unionvector[3] --- 3 компоненті единичного вектора скорости, 
	//					которім задается направление движения частиці
	// vectorvalue --- Абсолютное значение вектора скорости
// Возвращает:
	// vvector[3] --- переопределенніе значения компонент скорости иона
	// булевское значение :
			// 1 - программа завершилась с проблемой
			// 0 - программа завершилась нормально
	// msg --- сообщенеи о ходе віполнения работі
int velocityDefine(double unionvector[3], double vectorvalue, double *vvector, char *msg[]);


int main( int argc, TCHAR *argv[] )
{
	int i;
	double a;
	double step = 1.;
	double low = 1.;
	double up = 1.1;
	char *msg;	// сообщение о ходе работі программі
	int countparmtrs = 6;	// количество параметров задачи
	double *setparmtrs;		// массив под параметрі задачи
	char **nameparmtrs;
	double unionvectorV[3];
	double vectorV[3];
	char executeProg[] = "mpiexec";
	char executePrmtrs[] = " -n 4 ";
	char progName[] = "D:\\01_Khelemelia\\02_IAP\\04_Projects\\ELiMA\\EliMA_WORK\\ELiMA_Programming\\ELiMA\\x64\\Debug\\ELiMA.exe";
	char cmdLine[BUFSIZ];
	// about TCHAR and LPTSTR read on:
		// http://habrahabr.ru/post/164193/
//	TCHAR *adrProgram = _T(" cpi.exe"); // an adress of triggered programm;
									
	TCHAR *adrProgram = _T(FILEPATH); 
	TCHAR *parameters = _T("\" -n 4 \"");	// parameters of exec program;
	TCHAR *adrExecProg = _T("\"C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe\"");
///	TCHAR* adrExecProg = _T("\"C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe\"");
										// an adress of the executable program;
	LPTSTR szCmdline;		// a full command name;

    STARTUPINFO si;
    PROCESS_INFORMATION pi;
	ZeroMemory( &si, sizeof(si) );
    si.cb = sizeof(si);
    ZeroMemory( &pi, sizeof(pi) );

	//////////////////////////////////////////////////////////////
	//// задание параметров задачи////////////////////////////////
	// Віделение памяти под параметрі
	setparmtrs = (double*) malloc(countparmtrs*sizeof(double));
	// Віделение памяти под название параметров
	nameparmtrs = (char**) malloc(countparmtrs*6*sizeof(char)); 
	// Names have to use to 6 characters
	nameparmtrs[0] = "Vx";
	nameparmtrs[1] = "Vy";
	nameparmtrs[2] = "Vz";
	nameparmtrs[3] = "|V|";
	nameparmtrs[4] = "Te_per";
	nameparmtrs[5] = "Te_par";
	//nameparmtrs[6] = "Hx";
	//nameparmtrs[7] = "Hy";
	//nameparmtrs[8] = "Hz";
	//nameparmtrs[9] = "|H|";

	// Union velocity vector nV
	//// perpendicular
	//unionvectorV[0] = 0.; //nVx
	//unionvectorV[1] = 1.; //nVy
	//unionvectorV[2] = 0.; //nVz
	// parallel
	unionvectorV[0] = 0.; //nVx
	unionvectorV[1] = 0.; //nVy
	unionvectorV[2] = 1.; //nVz
	//// angle = 45 degrees
	//unionvectorV[0] = 0.; //nVx
	//unionvectorV[1] = sqrt(2.)/2.; //nVy
	//unionvectorV[2] = sqrt(2.)/2.; //nVz
	//// angle = 30 degrees
	//unionvectorV[0] = 0.; //nVx
	//unionvectorV[1] = sqrt(3.)/2.; //nVy
	//unionvectorV[2] = sqrt(1.)/2.; //nVz

	// assign absolute value of velocity
	setparmtrs[3] =1;	// |V|

	// Define value of velocity vector component 
	if(velocityDefine(unionvectorV, setparmtrs[3], vectorV, &msg)) { // Error 
		fprintf(stdout,"%s", msg);
		_getch();
		return 1;
	}
	// copy value of velocity vector component to array of parameters
	setparmtrs[0] = vectorV[0]; // Vx
	setparmtrs[1] = vectorV[1];	// Vy
	setparmtrs[2] = vectorV[2];	// Vz
	// assing value of electron component temperature
	setparmtrs[4] = 0.001;	// Te_perp, eV
	setparmtrs[5] = 0.001;	// Te_parallel, eV
	
	////????add magnetic field components
	//setparmtrs[6] = 0; // Hx
	//setparmtrs[7] = 0;	// Hy
	//setparmtrs[8] = 0;	// Hz
	//setparmtrs[9] = 0; // |H|


	//print a name and a value of parameters
	for(i = 0; i < countparmtrs; i++) {
		fprintf(stdout, "%s\t%e\n", nameparmtrs[i], setparmtrs[i]);
	}

	if (dataFileWrite(DATAFILEPATH, low, up, step, setparmtrs, countparmtrs, &msg)) { 
		fprintf(stdout,"%s", msg);
		_getch();
		return 1;
	}

	strcpy_s(cmdLine,sizeof cmdLine, executeProg);
	strcat_s(cmdLine, sizeof cmdLine, executePrmtrs);
	strcat_s(cmdLine, sizeof cmdLine, progName);
	fprintf(stdout, "%s", cmdLine);
	system(cmdLine);
	szCmdline = (TCHAR *)malloc(400*sizeof(TCHAR));
	// Create a full command name consist of:
		// 1-st --- an adress of the executable program (adrExecProg)
		// 2-nd --- parameters of exec program (parameters)
		// 3-rd --- an adress of triggered programm (adrProgram)
	wsprintf(szCmdline, L"%s%s%s", adrExecProg, parameters, adrProgram);


	if (!CreateProcess( NULL,   // No module name (use command line)
		szCmdline ,      // Command line
		NULL,           // Process handle not inheritable
		NULL,           // Thread handle not inheritable
		FALSE,          // Set handle inheritance to FALSE
		0,              // No creation flags
		NULL,           // Use parent's environment block
		NULL,           // Use parent's starting directory 
		&si,            // Pointer to STARTUPINFO structure
		&pi )           // Pointer to PROCESS_INFORMATION structure
	) {
		fprintf(stdout, "%s", pi);
		Sleep(10);				// подождать
		//MessageBox(HWND_DESKTOP, "Unable to start program", "", MB_OK);
		TerminateProcess(pi.hProcess,NO_ERROR);	// убрать процесс
	}
	// Wait until child process exits.
	WaitForSingleObject( pi.hProcess, INFINITE );

	// Close process and thread handles. 
	CloseHandle( pi.hProcess );
	CloseHandle( pi.hThread );

	free(setparmtrs);
	free(szCmdline);
	fprintf(stdout, "\a"); fflush(stdout);

	_getch();
}

int dataFileWrite(char *file_name, double low, double up, double step, double *parameter, int countparameter, char *msg[]) {
	FILE *data_file;
	double a;
	if ( (up-low) <= 0 ) { // некоректные значения границ
		*msg = "dataFileWrite(Errorcode_0AA): Wrong data of parameter low and up\n";
		return 1;
	}

	if ( (step) <= 0 ) { // некоректное значение шага  
		*msg = "dataFileWrite(Errorcode_0AB): Uncorrect parameter of step\n";
		return 1;
	}

	if(fopen_s(&data_file, file_name, "w+")){ // не открывается файл для записи
		fprintf(stdout, "%s", file_name);
		*msg = "dataFileWrite(Errorcode_0AC): Can't create datafile\n";
		return 1;
	}
	// печать параметров в файл
	fprintf(data_file, "V0x\tV0y\t\tV0z\t\t|V|\tTe_Per\tTe_Par\n"); // определители параметров
	for (a = low; a <= up; a+=step){
		int i;
		for (i = 0; i < countparameter; i++) {
			fprintf(data_file, "%f\t", parameter[i]); // печать значений 
		}											// параметров в файл
		fprintf(data_file, "\n"); 
	}

	// закрытие файла
	fclose(data_file);
	*msg = "dataFileWrite(Errorcode_0AD): program has finished OK";
	return 0;
}

int velocityDefine(double unionvector[3], double vectorvalue, double *vvector, char *msg[]) {
	int i;
	double precision = 0.000000001;
	// Абсолютное значение единичного вектора = 1
	double valueunionV = 0;
	for (i = 0; i < 3; i++) {
		valueunionV += unionvector[i]*unionvector[i];
	}
	if ((fabs(valueunionV - 1) > precision)) {
		*msg = "velocityDefine (Errorcode_0BA): Uncorrect value of parameters. Union velocity vector";
		return 1;
	}

	// Абсолютное значение вектора скорости больше 0
	if (vectorvalue <= 0) {
		*msg = "velocityDefine (Errorcode_0BB): Uncorrect value of parameters. Value of velocity vector";
		return 1;
	}

	// задание значений компонент скорости
	for (i = 0; i < 3; i++) {
		vvector[i] = unionvector[i]*vectorvalue;
	}
	*msg = "velocitydefine (errorcode_0bc): program finished ok";
	return 0;
}
