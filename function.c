#include "Function.h"
#include "InputData.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h> //malloc
#include <gsl/gsl_sf_bessel.h>

// ÏÐÈÍÈÌÀÅÒ: äåéñòâèòåëüíûé àðãóìåíò, òèïà double
// îöåíèâàåò ïîëó÷åííûé àðãóìåíò
// ïîñëå ÷åãî îïðåäåëÿåò íóæíîå ïðèáëèæåíèå
// ÂÎÇÂÐÀÙÀÅÒ: ïðèáëèæåííîå çíà÷åíèå êîìïëåêñíîé ôóíêöèè îøèáîê
//				îò äåéñòâèòåëüíîãî àðãóìåíòà, óìíîæåíîé íà 
//				exp(-x^2)*erfi(x);
double DispertionFunctionApproximation (double x) {
	const double a = 1.5;
	const double b = 2.6;
	if (fabs(x) <= a) {
		const double a0 = 2./sqrt(PI);
		const double a1 = 1.;
		const double a2 = 0.33333333333333333333333333333333; //1/3
		const double a3 = 0.1; // 1/10
		const double a4 = 0.02380952380952380952380952380952; //1/42
		const double a5 = 0.00462962962962962962962962962963; //1/216
		const double a6 = 7.5757575757575757575757575757576e-4; //1/1320
		const double a7 = 1.0683760683760683760683760683761e-4; //1/9360
		const double a8 = 1.3227513227513227513227513227513e-5; //1/75600

		return exp(-x*x)*a0*x*(a1 + a2*pow(x,2) + a3*pow(x,4) + a4*pow(x,6) + 
			a5*pow(x,8) + a6*pow(x,10) + a7*pow(x,12) + a8*pow(x,14));
	} else if (fabs(x) >= b) {
		const double a0 = 2./sqrt(PI);
		const double a1 = 0.5; //1/2
		const double a2 = 0.25; //1/4
		const double a3 = 0.375; // 3/8
		const double a4 = 0.9375; //15./16
		const double a5 = 3.28125; //35*24/256 -> 840/256
		const double a6 = 14.765625; //63*120/512 -> 7560/512
		const double a7 = 81.2109375; //231*720/2048 -> 166320/2048
		const double a8 = 527.87109375; //429*5040/4096 -> 2162160/4096

		return a0*(a1 + a2/pow(x,2) + a3/pow(x,4) + a4/pow(x,6) + 
			a5/pow(x,8) + a6/pow(x,10) + a7/pow(x,12) + a8/pow(x,14))/x;
	} else if (x > 0){
		const double a0 = 2./sqrt(PI);
		const double a4 = 0.092;
		const double a1 = 3.05; //1/2
		const double a2 = 0.0121; //1/4
		const double a3 = 0.2253; // 3/8

		return a0*(a4*pow((x-a1),2) - a2*x + a3);
	} else if (x < 0) {
		const double a0 = 2./sqrt(PI);
		const double a4 = 0.092;
		const double a1 = 3.05; //1/2
		const double a2 = 0.0121; //1/4
		const double a3 = 0.2253; // 3/8

		return a0*(-a4*pow((x+a1),2) - a2*x - a3);
	} else {
		fprintf(stdout, "Error value of incoming parameter.");
		fprintf(stdout, "Prog: double DispertionFunctionApproximation (double x).");
		fprintf(stdout, "File: function.c\n");
		exit(1);
	}
}

// ÏÐÈÍÈÌÀÅÒ: 3 äåéñòâèòåëüíûõ àðãóìåíòà, òèïà double
	// ñôåðè÷åñêèå êîîðäèíàòû: phi, theta, rho
// ÂÎÇÂÐÀÙÀÅÒ: çíà÷åíèå ôóíêöèè, ïðåäñòàâëåííîé â ðàáîòå 
	// È.À. Ëàðêèí. Ïðîõîæäåíèå ÷àñòèö ÷åðåç ïëàçìó ÆÝÒÔ , 37(1), 264-272 (1959).
	// ôîðìóëà 27, îáåçðàçìåðåííàÿ ñ âîñòàíîâëåííûìè åäèíèöàìè e=hbar=m_e=1
double Larkin27(double x, double k, double* parameter){
	double delta = parameter[2];
	double tau = parameter[3]; 
	double Wk;
	double RePhi, ImPhi;
	double UpP, LowS, PreP, PreS; 
	double Wx1, Wx2;
	double A_Khel = 1.;			// Polarization coefficient by Khelemelia
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
	double k2;
	Wk = (x) ;//+ delta*k*mM; // defined as part of: \vec{k}*\vec{n_{1}} - DELTA*k^2*m_{e}/M_{1}
	PreP = A_Khel*sqrt(PI)/sqrt(2.*tau)/delta/k/2.;		// Palarization operator coefficient
	Wx1 = (Wk-delta*k)/sqrt(2.*tau);					// arguments of dispertion function
	Wx2 = (Wk+delta*k)/sqrt(2.*tau);
	UpP = exp(-Wx1*Wx1);
	ImPhi = PreP*UpP*(1.-exp(-2.*delta*k*Wk/tau));		// Imaginary Part of Polarization Operator
	RePhi = PreP*(DispertionFunctionApproximation(Wx2) 
		- DispertionFunctionApproximation(Wx1));		// Real Part of Polarization Operator
														
	PreS = A_Khel/sqrt(2*PI*tau)/delta;					// Energy Losses Integral Coefficient
	k2=k*k;
	LowS = 1/((k2+RePhi)*(k2+RePhi)+ImPhi*ImPhi);

	return PreS*Wk*k2*UpP*LowS;
}


// ÏÐÈÍÈÌÀÅÒ: 2 äåéñòâèòåëüíûõ àðãóìåíòà, òèïà double
	// ñôåðè÷åñêèå êîîðäèíàòû: cos(theta) = -1..1, rho = 0..infinity
// ÂÎÇÂÐÀÙÀÅÒ: çíà÷åíèå ôóíêöèè, ïðåäñòàâëåííîé â ðàáîòå 
	// È.À. Ëàðêèí. Ïðîõîæäåíèå ÷àñòèö ÷åðåç ïëàçìó ÆÝÒÔ , 37(1), 264-272 (1959).
	// ôîðìóëà 27, îáåçðàçìåðåííàÿ ñ âîñòàíîâëåííûìè åäèíèöàìè e=hbar=m_e=1
double ELI(double x, double k, double* parameter){
	double delta = parameter[2];
	double tau = parameter[3]; 
	double Wk;
	double RePhi, ImPhi;
	double UpP, LowS, PreP, PreS; 
	double Wx1, Wx2;
	double A_Khel = 1.;			// Polarization coefficient by Khelemelia
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
	double k2;
	Wk = (x) ;//+ delta*k*mM; // defined as part of: \vec{k}*\vec{n_{1}} - DELTA*k^2*m_{e}/M_{1}
//	Wk = (x) + delta*k*mM; // defined as part of: \vec{k}*\vec{n_{1}} - DELTA*k^2*m_{e}/M_{1}
	PreP = A_Khel*sqrt(PI)/sqrt(2.*tau)/delta/k/2.;		// Palarization operator coefficient
	Wx1 = (Wk-delta*k)/sqrt(2.*tau);					// arguments of dispertion function
	Wx2 = (Wk+delta*k)/sqrt(2.*tau);
	UpP = exp(-Wx1*Wx1);
	ImPhi = PreP*UpP*(1.-exp(-2.*delta*k*Wk/tau));		// Imaginary Part of Polarization Operator
	//ImPhi = PreP*UpP*(1.);		// Imaginary Part of Polarization Operator
	RePhi = PreP*(DispertionFunctionApproximation(Wx2) 
		- DispertionFunctionApproximation(Wx1));		// Real Part of Polarization Operator
														
	PreS = A_Khel/sqrt(2*PI*tau)/delta;					// Energy Losses Integral Coefficient
	k2=k*k;
	LowS = 1/((k2+RePhi)*(k2+RePhi)+ImPhi*ImPhi);

	//return 2*PreS*Wk*k2*UpP*LowS;
	return PreS*Wk*k2*ImPhi*LowS/PreP; // áåç òåìïåðàòóðíîãî ñëàãàåìîãî
}

double function2D (double x, double y, double *par) {
	//return 1.;
	return 605*y/((1+120*(1-y))*((1+120*(1-y))*(1+120*(1-y))+25*x*x*y*y));//(0..1, 0..1, 0..1) = 0.1047591113142868D  01
	//return 1/(x+y+z)/(x+y+z);//(0..1, 0..1, 0..1) = 0.8630462173553432D   00
	//return cos(x + y); // -0.4 (0..3PI)(0..3PI)
	//return exp(-(x*x+y*y+z*z))*(x+y+z+1)/((x+y+z+1)*(x+y+z+1)+1);
	//return x/pow((1-x*x), 0.5);
	//return pow(x, 0.5)*pow(z, 0.5)*pow(y, 0.5)*log(x*y*z); // -.5925925926 (0..1)
	//return pow(x, 0.5)*log(x); // -0.444444 (0..1)
	//return exp(x)*cos(x); //1.905238691 (0..PI/2)
}

double function3D (double x, double y, double z, double *par) {
	//return gsl_sf_bessel_In(1, 1.);
	//return 605*y/((1+120*(1-y))*((1+120*(1-y))*(1+120*(1-y))+25*x*x*y*y));//(0..1, 0..1, 0..1) = 0.1047591113142868D  01
	//return 1/(x+y+z)/(x+y+z);//(0..1, 0..1, 0..1) = 0.8630462173553432D   00
	//return cos(x + y); // -0.4 (0..3PI)(0..3PI)
	//return exp(-(x*x+y*y+z*z))*(x+y+z+1)/((x+y+z+1)*(x+y+z+1)+1);
	//return x/pow((1-x*x), 0.5);
	//return pow(x, 0.5)*pow(z, 0.5)*pow(y, 0.5)*log(x*y*z); // -.5925925926 (0..1)
	//return pow(x, 0.5)*log(x); // -0.444444 (0..1)
	//return exp(x)*cos(x); //1.905238691 (0..PI/2)
}
// Âîçðàùàåò çíà÷åíèå ïàðàìåòðà delta_0, ðàâíîãî
// delta_0 = hbar*omega_P/(2*T_e,perp)*tau_perp
// pre = hbar*omega_P/2
// Plank`s const, devided on 2*Pi
	// hbar = 6.582 E-16 eV*c
// Plasma frequency - for HESR 
	// omega_P = 2.9 E8 c^(-1)	
// const value - velocity of the light
	// c = 2.99792 E10 cm/c
// const value  - energy of electron
	// Ee = m_e*c^2 =  5.11 E5 eV
//// Ê ðàñ÷åòàì ïî Ïàðõîì÷óêó
//// ïëîòíîñòü n_e = 10^7 cm^(-3)
//double delta (double Te_per) {
//	double delta_0 = 5.89089*pow(10.,-8);	// Pre ïîäñ÷èòàí ïðåäâàðèòåëüíî
//	return delta_0/Te_per;
//}
// Ê ðàñ÷åòàì ïî HESR
double delta (double Te_per) {
	double delta_0 = 9.5439*pow(10.,-8);	// Pre ïîäñ÷èòàí ïðåäâàðèòåëüíî
	return delta_0/Te_per;
}

// Âîçðàùàåò çíà÷åíèå ïàðàìåòðà tau_i, ðàâíîãî
// tau_i = tau_0*(V_0/V_i)^2
// tau_0 = 1/(m_e*c^2)*c^2/V_0^2
// ñêîðîñòü ñâåòà â âàêóóìå
	// ñ = 2.99792 E10 cm/c
// ýíåðãèÿ ïîêîÿ ýëåêòðîíà
	// Ee = m_e*c^2 =  5.11 E5 eV
// íîðìèðîâî÷íàÿ ñêîðîñòü
	// V_0 = 10 E6 cm/c
double tau (double ViV0, double Te) {
	double tau_0 = 1.758811E3;	// ïîäñ÷èòàí ïðåäâàðèòåëüíî
	return tau_0*Te/ViV0/ViV0;
}

// Âîçðàùàåò çíà÷åíèå ïàðàìåòðà tau_per, ðàâíîãî
// tau_per = Te_per/me/Ve_per^2
// Ïîñêîëüêó ïðåäïîëàãàåòñÿ, ÷òî 
// Te_per = me*Ve_per^2,
// òîãäà tau_perpendicular = 1 !!âñåãäà!!
double tau_perpendicular () {
	return 1.;
}
// Âîçðàùàåò çíà÷åíèå ïàðàìåòðà tau_per, ðàâíîãî
// tau_par = Te_par/me/Ve_per^2
double tau_parallel (double Te_par, double Te_per) {
	return Te_par/Te_per;
}


// ÏÐÈÍÈÌÀÅÒ: 3 äåéñòâèòåëüíûõ àðãóìåíòà, òèïà double
// ÂÎÇÂÐÀÙÀÅÒ: çíà÷åíèå ôóíêöèè, ïðåäñòàâëåííîé â ðàáîòå 
	// Õåëåìåëè À.Â. Âëèÿíèå àíèçîòðîïèè íà ýíåðãåòè÷åñêèå ïîòåðè... 
double ELIA(double gx, double gy, double gz, double* parmtr) {
	// êîìïîíåíòû áåçðàçìåðíîãî èìïóëüñà íàëåòàþùåé ÷àñòèöû
	double sx = parmtr[0];
	double sy = parmtr[1];
	double sz = parmtr[2];

	// àáñîëþòíîå çíà÷åíèå âåêòîðà èìïóëüñà áåçðàçìåðíîãî
	double V = parmtr[3];

	// ïîïåðå÷íàÿ è ïðîäîëüíàÿ òåïåðàòóðû ýëåêòðîííîãî ãàçà
	double tau_per = tau_perpendicular();
	double tau_par = tau_parallel(parmtr[5], parmtr[4]);

	double Delta = delta(parmtr[4]);

	double W;
	double ge, g;
	double ReKappa, ImKappa;
	//	double UpP, LowS, PreP, PreS; 
	double Xi1, Xi2;
	double A_Khel = 1.;			// Polarization coefficient by Khelemelia
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
	double k2;
	double PreKappa;
	double PreInt, PostInt, Up;
	g = sqrt(gx * gx + gy * gy + gz * gz);
	ge = sqrt(gx * gx + gy * gy + gz * gz * tau_par / tau_per);

	W = (gx * sx + gy * sy + gz * sz) - Delta * mM * g * g;

	if ((gx == 0) && (gy == 0) && (gz == 0)) {
		return 0;
	}
	else {
		// àðãóìåòí³ äèñïåðñèîíí³õ ôóíêöèé
		Xi1 = (W + Delta * g * g) / sqrt(2 * tau_per) / ge; // Xi+;
		Xi2 = (W - Delta * g * g) / sqrt(2 * tau_per) / ge; // Xi+;

		Up = exp(-Xi2 * Xi2);
		PreKappa = A_Khel * sqrt(PI / 2 / tau_per) / Delta / (g * g) / ge / 2;
		ImKappa = PreKappa *
			Up * (1 - exp(-2 * Delta * W * g * g / ge / ge / tau_per));
		ReKappa = PreKappa *
			(DispertionFunctionApproximation(Xi1)
				- DispertionFunctionApproximation(Xi2));
		PreInt = 1 / Delta / sqrt(2 * PI * tau_per) / (2 * PI);
		PostInt = 1 / ((1 + ReKappa) * (1 + ReKappa) + ImKappa * ImKappa);
		fprintf(stdout, "num %d funct is %e\n", 1, PreInt);
		return PreInt * W * Up * PostInt / (g * g * g * g) / (ge);
	}
}

// ÏÐÈÍÈÌÀÅÒ: 3 äåéñòâèòåëüíûõ àðãóìåíòà, òèïà double
// ÂÎÇÂÐÀÙÀÅÒ: çíà÷åíèå ôóíêöèè, ïðåäñòàâëåííîé â ðàáîòå 
	// Õåëåìåëè À.Â. Ïîòåðè â çàìàãíè÷åííîì ýëåêòðîííîì ãàçå... 

	// h=Omega_H/Omega_P is dimless parameter of magnetic field
	// beta = 2*delta_0*h/tau_perp
	
	// äèýëåêòðè÷åñêàÿ âîñïðèèì÷èâîñòü ïðåäñòàâëåíà â âèäå ðÿäà
	// ïî ïàðàìåòðó a^2 = q_t^2 * Delta/h 

double ELiMA(double gx, double gy, double gz, double* parmtr){
		// it needs to obtaine this parameter for all new jobs 
	
	// a difference between Landau levels
	double s ; 
	int N = 4;

	// components of the union velocity vector V = V_i/V_0
	double Sx = parmtr[0]; 
	double Sy = parmtr[1];
	double Sz = parmtr[2];

	// value of the union velocity vector V = V_i/V_0
	double V = parmtr[3];

	// a perpendicular and a parallel components of the temperature of the electron gas
	// tau_per = tau_per / tau_per; = 1; 
	double tau_per = tau_perpendicular();
	double tau_par = tau_parallel(parmtr[5], parmtr[4]);

//	double h = 20000;
	double h = 100000;

	double tau;
	double W;		// dimmles frequency
	double beta;	// magnet argument of exponent
	double a2;		// parameter of Lambda-function, 
						// it is always used in form a^{2N}	
	double PreSusc;	// prefix before susseptibility
	double Delta = delta(parmtr[4]);		// parameter delta = hbar W_0 / 2m_eV_0^2 	
	double g;
	double ReKappa, ImKappa;
	double Xi1, Xi2;

//	double A_Khel = 1.;			// Polarization coefficient by Khelemelia (h=0)
//	double A_Lark = -1.;			// Polarization coefficient by Larkin
//	double A_Akhi = 2.;			// Polarization coefficient by Akhiezer
//	double AB = 4/PI;				// Polarization coefficient by Khelemelia (h != 0)

	double PreInt, PostInt;
	double  Up = 0, Up0;
	double Pi = 3.14159265;
	double TemperatureFactor;		// a temperature factor 1-exp(-2 w delta0/tau)
	double argBess, preBess;
	double ImMagnetPart = 0.;
	double ImMagnetPartTilde = 0.;
	double ReMagnetPart = 0.;
	double preScaledBess; 
	double exp_beta, exp_beta2;
	double ScaledBesselI_0, ScaledBesselI_1;
	double ScaledBesselInApprox;
	double ScaledBessel2 = 0, ScaledBessel1=0;
	double tmpRe0, tmpIm0, tmpImTilde0;

	//for (int i = 0; i < 10; i++) {
	//	fprintf(stdout, "par[%d] = %f\n", i, parmtr[i]); fflush(stdout);
	//}
	////////////////////////////////////////
	// General position
	g = sqrt(gx*gx + gy*gy + gz*gz);		// a dimless wave vector	
	W = (gx*Sx + gy*Sy + gz*Sz) - Delta*mM*g*g;		// a dimless frequency, determined as hw/hw_{0} = (E_{p} - E_{p-hk})/hw_{0}
	tau = tau_par/3. + 2.0*tau_per/3.0;		// temperature

	if (h == 0) { // procedure calculate non magnetic case (h = 0)
		return ELIA(gx, gy, gz, parmtr);
	} 

	if ((gz == 0) || ((gx == 0) && (gy == 0) && (gz == 0))) {  // to exclude numerical divergence
		return 0;
	} else {
		//PreSusc = 1/4.*sqrt(Pi)/sqrt(2*tau_par);  //a prefix of susceptibility. It's no variables of integration
		PreSusc = 1/4./Delta*sqrt(Pi)/sqrt(2*tau_par);  //a prefix of susceptibility. It's no variables of integration
		PreInt = PreSusc/Pi/Pi*10000;	

		Xi1 = (W +  Delta*gz*gz)/sqrt(2*tau_par)/gz; // variable 1 of exponenta;
		Xi2 = (W -  Delta*gz*gz)/sqrt(2*tau_par)/gz; // variable 2 of exponenta;

		a2 = (gx*gx + gy*gy)*Delta/h;	// an argument of LambdaFunction
											// the parameter is inverse to the parameter of the magnetic field h = W_H/W_0
		beta = 2.*Delta*h/tau_per;		// an argument of exp
											// the argument is proportinal to the parameter of the magnetic field h = W_H/W_0  
		exp_beta2 = exp(-beta/2.);
		exp_beta = exp(-beta);

		// for GS 
		//argBess = 2.*a2*(exp(-beta/2.)/(1.-exp(-beta)));
		argBess = 2.*a2*(exp_beta2/(1.-exp_beta));
		//preBess	= exp(-a2*(1.+exp(-beta))/(1.-exp(-beta)));	
		//preScaledBess = exp(-a2*(1.-exp(-beta/2.))/(1.+exp(-beta/2.)));	
		preScaledBess	= exp(-a2*(1.-exp_beta2)/(1.+ exp_beta2));	

		//Scaled modified bessel function exp(-|x|)*BesselI(n, x)
		//ScaledBesselI_0 = gsl_sf_bessel_I0_scaled(argBess);
		//ScaledBesselI_1 = gsl_sf_bessel_I1_scaled(argBess);
		
		for (s = 0; s<=10; s++) {
			double theta1Hp, theta1Hm, theta2Hp, theta2Hm;
			double XiSH;
			double RS;
//			double BesselInApprox;
			XiSH = s*h/sqrt(2.*tau_par)/gz;
			theta1Hp = Xi1 + XiSH;
			theta1Hm = Xi1 - XiSH;
			theta2Hp = Xi2 + XiSH;
			theta2Hm = Xi2 - XiSH;			

			RS = preScaledBess * gsl_sf_bessel_In_scaled(s, argBess);

			 {
				ImMagnetPart += - RS *
						(exp(-s*beta/2.)*(exp(-theta1Hm*theta1Hm)  - exp(-theta2Hp*theta2Hp)) +
						exp(s*beta/2.)*(exp(-theta1Hp*theta1Hp)  - exp(-theta2Hm*theta2Hm)));
				ReMagnetPart += RS *
					(exp(-s*beta/2)*(DispertionFunctionApproximation(theta1Hm) -
								DispertionFunctionApproximation(theta2Hp)) +
					exp(s*beta/2)*(DispertionFunctionApproximation(theta1Hp)-
								DispertionFunctionApproximation(theta2Hm)));
				ImMagnetPartTilde += RS * 
					(exp(s * beta / 2.)*exp(-theta2Hp * theta2Hp) +
						exp(-s * beta/2.) * exp(-theta2Hm * theta2Hm));

				//ImMagnetPart += -RS * exp(s*beta/2.)*
				//		(exp(-s*beta)*(exp(-theta1Hm*theta1Hm)  - exp(-theta2Hp*theta2Hp)) +
				//		(exp(-theta1Hp*theta1Hp)  - exp(-theta2Hm*theta2Hm)));

				//ImMagnetPartTilde +=  RS * exp(s*beta/2.)*
				//		(exp(-theta2Hp*theta2Hp) +
				//		 exp(-s*beta)*exp(-theta2Hm*theta2Hm));

				//ReMagnetPart += RS * exp(s*beta/2)*
				//	(exp(-s*beta)*(DispertionFunctionApproximation(theta1Hm) -
				//				DispertionFunctionApproximation(theta2Hp)) +
				//	(DispertionFunctionApproximation(theta1Hp)-
				//				DispertionFunctionApproximation(theta2Hm)));
				/*if (fabs(W) <= 1e-16) {
					ImMagnetPartTilde +=  (tau/tau_par)* RS *(
						exp(-s*beta/2.)*exp(-theta2Hm*theta2Hm)*(1+s*h/gz/gz/Delta) +
						exp(s*beta/2.)*exp(-theta2Hp*theta2Hp)*(1-s*h/gz/gz/Delta));
				} else*/ /*{
					ImMagnetPartTilde = ImMagnetPart;
				}*/
				if (s == 0) {
					tmpIm0=ImMagnetPart;
					tmpRe0=ReMagnetPart;
					tmpImTilde0=ImMagnetPartTilde;
				}
			}

		}
		/*if (fabs(W) <= 1e-16) {
			TemperatureFactor = 1.;
		} else*/{
			TemperatureFactor = 1/(1.-exp(-2.*Delta*W/tau));
		}
		//////////////////////////////////////////
		//to get the PostInt part
		// An imaginary part of susceptibility is
		ImMagnetPart=ImMagnetPart-tmpIm0/2.;
		ReMagnetPart=ReMagnetPart-tmpRe0/2.;
		ImMagnetPartTilde=ImMagnetPartTilde-tmpImTilde0/2.;

		ImKappa = PreSusc*ImMagnetPart/g/g/fabs(gz);

		// an real part of susceptibility is
		ReKappa = PreSusc*ReMagnetPart/g/g/fabs(gz);

		PostInt = 1/((1+ReKappa)*(1+ReKappa) + ImKappa*ImKappa);

		return PreInt*ImMagnetPartTilde*((gx*Sx + gy*Sy + gz*Sz) - Delta*mM*g*g)/(g*g*g*g*fabs(gz))*PostInt;
		//fprintf(stdout, "pre %e and up %e and post %e\n", PreInt, ImMagnetPart, PreInt* ((gx* Sx + gy * Sy + gz * Sz) - Delta * mM * g * g) / (g * g * g * g * gz) * PostInt);

		//return 1;
	}
}



double velocity_ion (double velocity_ion_norm) {
	return velocity_ion_norm * V_0();
}

double velocity_electron(double temperature) {
	// T_e = m_e v_e^2 -> v_e = sqrt(T_e / (m_e * c^2)) * c = sqrt(T_e) * c / sqrt(m_e c^2) 
	// [T_e] = eV 
// const value - velocity of the light
	// c = 2.99792458 E10 cm/s
// const value  - energy of electron
	// Ee = m_e*c^2 =  0.51099895 E6 eV
	// v_e = 4.19382881 e7 * sqrt(T_e) cm/s
	return 4.19382881e7 * sqrt(temperature);
}

double V_0() {
	//V_0 = 1e6 cm/s 
	// It responds electron temperature velocity = 1.33e-6 cm/s (T_e \sim 1e-3 eV - experimental temperature of electrons)
	return 1e6; 
}

double average_velocity_square(const double temperature_Eperp, const double temperature_Eparall){
	return pow(velocity_electron(temperature_Eperp), 2) + pow(velocity_electron(temperature_Eparall), 2);
}

double average_frequency (const double frequency_Eplasma, const double frequency_Ecyclotron){
	return sqrt(pow(frequency_Eplasma, 2) + pow(frequency_Ecyclotron, 2));
}

double frequency_norm(double frequency, const double frequency_Eplasma, const double frequency_Ecyclotron) {
	return frequency / average_frequency(frequency_Eplasma, frequency_Ecyclotron);
};

double average_tau(const double temperature_Eperp, const double temperature_Eparall) {
	return sqrt(temperature_Eperp * temperature_Eperp + temperature_Eparall * temperature_Eparall);
}

double tau_perp_norm(double tau_perp, double tau_parall) {
	return tau_perp / average_tau(tau_perp, tau_parall);
}

double tau_parallel_norm(double tau_perp, double tau_parall) {
	return tau_parall / average_tau(tau_perp, tau_parall);
}

/*recieve const variable :
 delta = hbar*omega_average/(2*m_e*c^2)*(c^2/V_average^2)
 reduced Plank`s const, devided on 2*Pi
	 hbar = 6.58211928 E-16 eV*s
 const value - velocity of the light
	 c = 2.99792458 E10 cm/s
 const value  - energy of electron
	 Ee = m_e*c^2 =  0.51099895 E6 eV
	  p = hbar/2/(m_e c^2)*c^2  = 5.78838155 E-1 cm^2/s  
	 delta = p * omega_average / V_average^2
*/
double delta_norm(const double temperature_Eperp, const double temperature_Eparall, const double frequency_Eplasma, const double frequency_Ecyclotron) {
	const double pre = 5.78838189e-5; // pre = hbar/2/(m_e) = 0.0578838189244059*10^(-3) J*s/kg

	return pre * average_frequency(frequency_Eplasma, frequency_Ecyclotron) /
		average_velocity_square(temperature_Eperp, temperature_Eparall);
}

double ELiMA_H(double gx, double gy, double gz, double* set_parmtr) {
	// gx, gy, gz - 3D coords (components of a dimentionless wave vector)  
	// set_parmtr[0] - Vx/Vo
	// set_parmtr[1] - Vy/Vo
	// set_parmtr[2] - Vz/Vo
	// set_parmtr[3] - Vi/Vo
	// set_parmtr[4] - Te,perp [eV]
	// set_parmtr[5] - Te, parall [eV]
	// set_parmtr[6] - Wp
	// set_parmtr[7] - WHx
	// set_parmtr[8] - WHy
	// set_parmtr[9] - WHz
	// set_parmtr[10] - WH

	double result = 0;
	// Velocity components
	double Vxn, Vyn, Vzn, Vn;
	Vxn = set_parmtr[0];
	Vyn = set_parmtr[1];
	Vzn = set_parmtr[2];
	Vn = set_parmtr[3];
	
	// temperatures of an electron gas 
	double tauPerp_eV, tauParall_eV;
		tauPerp_eV = set_parmtr[4];
	tauParall_eV = set_parmtr[5];
	//normalized
	double tauPerpN, tauParallN;
	tauPerpN = tau_perp_norm(tauPerp_eV, tauParall_eV);
	tauParallN = tau_parallel_norm(tauPerp_eV, tauParall_eV);

	//plasma frequency
	double Wp, Wh;  // plasma and cyclotron frequencies
	Wp = set_parmtr[6];
	Wh = set_parmtr[10];
	//normalized
	double Wpn = frequency_norm(Wp, Wp, Wh);
	double Whn = frequency_norm(Wh, Wp, Wh);

	// a dimentionless wave vector
	double g = sqrt(gx * gx + gy * gy + gz * gz);

	// quantum parameter of task
	double Delta = delta_norm(Vn, tauPerp_eV, tauParall_eV, Wp, Wh);

	// dimentionless frequency determined as hw/hw_{0} = (E_{p} - E_{p-hk})/hw_{0}
	double W = gx * Vxn + gy * Vyn + gz * Vzn - Delta * mM * g * g;	
	
	// computation infinitesimal value
	double mathPrecision = 1e-16;

	//partial cases:
	/*//1) anisotropic nonmagnetic 
	if (Whn <= mathPrecision) { // procedure calculate non magnetic case (h = 0)
		return ELIA(gx, gy, gz, set_parmtr);
	}*/
	/*// 2) isotropic nonmagnetic
	if (Whn <= mathPrecision && (TauPerpN - TauParallN) <= mathPrecision) {
		return ELI(g, set_parmtr); // don't performed
	}*/

	// see toNumCalc.pdf
	if ((gz <= mathPrecision) || (g<=mathPrecision)) {  // to exclude numerical divergence
		return 0;
	}
	
	//a prefix of susceptibility. It does not include variables of integration
	double PreSusc = Wpn * Wpn * sqrt(PI) / 4. / Delta / sqrt(2 * tauParallN); 
	
	//a prefix of dE/dt integral. There is no variable of integration 
	double PreInt = PreSusc / PI / PI;

	double Xi1, Xi2;
	Xi1 = (W + Delta * gz * gz) / sqrt(2 * tauParallN) / gz; // argument 1 of exponenta;
	Xi2 = (W - Delta * gz * gz) / sqrt(2 * tauParallN) / gz; // argument 2 of exponenta;

	// an argument of LambdaFunction
	// the parameter is inverse to the parameter of the magnetic field h = W_H/W_0
	double a2 = (gx * gx + gy * gy) * Delta / Whn;

	// an argument of exp
	// the argument is proportinal to the parameter of the magnetic field h = W_H/W_0  
	double beta = 2. * Delta * Whn / tauPerpN;		
	double exp_beta2 = exp(-beta / 2.);
	double exp_beta = exp(-beta);

	// magnetic field part

	// an argument of Bessel Function (GNU GSL)
	double argBess = 2. * a2 * (exp_beta2 / (1. - exp_beta));
	// prefix const of scaled modified Bessel function, In(z)*exp(-|z|);
	double preScaledBess = exp(-a2 * (1. - exp_beta2) / (1. + exp_beta2));

	
	double ImMagnetPart = 0;
	double ReMagnetPart = 0;

	int N = 0;
	//s is difference between Ladau levels (index of summation)
	for (int s = 0; s <= N; s++) {
		
		double XiSH;
		XiSH = s * Whn / sqrt(2. * tauParallN) / gz;

		double theta1Hp, theta1Hm, theta2Hp, theta2Hm;
		theta1Hp = Xi1 + XiSH;
		theta1Hm = Xi1 - XiSH;
		theta2Hp = Xi2 + XiSH;
		theta2Hm = Xi2 - XiSH;

		double RS = gsl_sf_bessel_In_scaled(s, argBess) * preScaledBess;

		ImMagnetPart += RS * 
			(exp(-s * beta/2.) * (exp(-theta1Hp * theta1Hp) - exp(-theta2Hm * theta2Hm)) +
				exp(s * beta / 2.) * (exp(-theta1Hm * theta1Hm) - exp(-theta2Hp * theta2Hp)));
		ReMagnetPart += RS * 
			(exp(-s * beta / 2.) * (DispertionFunctionApproximation(theta1Hp) -
				DispertionFunctionApproximation(theta2Hm)) +
				exp(s * beta / 2) * (DispertionFunctionApproximation(theta1Hm) -
					DispertionFunctionApproximation(theta2Hp)));
	}
	double TDenominator = 1 / (1 - exp(-2 * Delta * W / average_tau(tauPerpN, tauParallN)));
	double ImSucc = - PreSusc * 1 / g / g / gz * ImMagnetPart;
	double ReSucc = PreSusc * 1 / g / g / gz * ReMagnetPart;
	double postInt = 1 / ((1 - ReSucc) * (1 - ReSucc) + ImSucc * ImSucc);
	result = PreInt * W / g / g / g / g / gz * TDenominator * ImMagnetPart * postInt;
}
