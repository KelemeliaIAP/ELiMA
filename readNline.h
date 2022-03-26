#ifndef READNLINE_H
#define READNLINE_H

#include<conio.h>
#include<io.h>
#include <stdlib.h> // string
#include<stdio.h>
#include<string.h>


#define MSGLENGTH 500 
// функция readNline считывает 
// n-ую строку из входящего файла "in"
// возвращает массив типа double длинною count, 
// содержащийся в считаной строке
// е. в программе произошел сбой
//		возвращается значение 1
//		возращается запись ошибки msg
//	иначе:
//		возвращается значение 0
int readNline (FILE **in, unsigned int number, float *par, int count, char msg[]);

//parse file
// do:	- get headers from file
//		- get parameters from file
// in:	- file name
//		- pointer on a headers array
//		- pointer on a parameters array
//		- pointer on log message
// out: - array of headers (as pointer)
//		- array of parameters (as pointer)
//		- log msg (as pointer)
//		- return
//			<> 1 if procedure finishes as well as normal
//			<> 0 if procedure finishes with error
int parseFile(FILE** file, const int arraySize, char** headers[], double* parameters[], char *msg[]);
#endif