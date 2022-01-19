  #include <unistd.h>
  #include <stdio.h>
  #include <stdlib.h>
  #include <string.h>

 int main(int argc, char **argv)
  {
	char *argu[5];
//	char *argu1 = "-np";
//	char *argu3 = "/home/xvdm/01_Khelemelia/01_IAP/ELiMA/06_Programming/ELiMA(Linux).v.3.01/elima";
	
	// an instruction 
	argu[0] = "-np";
	// a number of processors that need to produce calculation
	argu[1] = argv[1];
	// a program to calculation
	argu[2] = "/home/xvdm/01_Khelemelia/01_IAP/ELiMA/06_Programming/ELiMA(Linux).v.3.01/elima";
	//number of task parameters	
	argu[3] = "7"; 
	// VELOCITY
		// an union velocity vector 
	argu[4]	= "0."; //x - component
	argu[5] = "0.";	//y - component
	argu[6] = "1."; //z - component (parallel to MF component)
		// absolutly value of velocity (in units of Ve)
	argu[7] = "1.";
	// ELECTRON TEMPERATURE (in eV)
		// a parallel component
	argu[8] = "0.01";
		// a perpendicular component
	argu[9] = "0.01";
	// MAGNETIC FIELD
		// absolutly value of magnetic field (in units of plasma frequency)
	argu[10] = "100.";	
	
//	char *argu = " -np 2 /home/xvdm/01_Khelemelia/01_IAP/ELiMA/06_Programming/./elima";
	//fprintf(stdout, "%s\n", argv[1]);
        //execlp("/home/xvdm/01_Khelemelia/01_IAP/ELiMA/06_Programming/elima", "elima", NULL);
	execlp("/usr/bin/mpirun", "mpirun",argu[0], argu[1], argu[2],argu[3], argu[4],argu[5], argu[6], argu[7],argu[8],argu[9],argu[10], NULL);
        printf("Return not expected. Must be an execlp error.n");
	return 0;
  }
