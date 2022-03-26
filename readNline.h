#ifndef READNLINE_H
#define READNLINE_H

#include<conio.h>
#include<io.h>
#include <stdlib.h> // string
#include<stdio.h>
#include<string.h>


#define MSGLENGTH 500 
// ������� readNline ��������� 
// n-�� ������ �� ��������� ����� "in"
// ���������� ������ ���� double ������� count, 
// ������������ � �������� ������
// �. � ��������� ��������� ����
//		������������ �������� 1
//		����������� ������ ������ msg
//	�����:
//		������������ �������� 0
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