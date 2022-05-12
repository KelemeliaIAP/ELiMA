// Function.h �������� ����������� ��������� ��������������� �������
#ifndef FUNCTION_H
#define FUNCTION_H

// Test 2-dimentional function, to verify calculation algorithm
// x,y - variables of integration
// par[] - set of integral external parameters
double function2D (double x, double y, double *par);

// Test 3-dimentional function, to verify calculation algorithm
// x,y,z - variables of integration
// par[] - set of integral external parameters
double function3D (double x, double y, double z, double *par);

// ���������: �������������� ��������, ���� double
	// ��������� ���������� ��������
	// ����� ���� ���������� ������ �����������
// ����������: ������������ �������� ����������� ������� ������
	// �� ��������������� ���������, ��������� �� 
	// exp(-x^2)*erfi(x);
double DispertionFunctionApproximation (double x);

// ���������: 3 �������������� ���������, ���� double
	// ����������� ����������: phi, theta, rho
// ����������: �������� �������, �������������� � ������ 
	// �.�. ������. ����������� ������ ����� ������ ���� , 37(1), 264-272 (1959).
	// ������� 27, ��������������� � ��������������� ��������� e=hbar=m_e=1
double Larkin27(double x, double k, double* parameter);

// ���������: 2 �������������� ���������, ���� double
	// ����������� ����������: cos(theta) = -1..1, rho = 0..infinity
// ����������: �������� �������, �������������� � ������ 
	// �.�. ������. ����������� ������ ����� ������ ���� , 37(1), 264-272 (1959).
	// ������� 27, ��������������� � ��������������� ��������� e=hbar=m_e=1
double ELI(double x, double k, double* parameter);

// ���������: �������� �������� ������������ �������� ����
// V_i/V_0
// ��������� �������� ��������� delta_0, �������
// delta_0 = hbar*omega_P/(2*m_e*c^2)*c^2/V_e^2
// pre = hbar*omega_P/(2*m_e*c^2)*c^2
// ���������� ������, ������������� �� 2*Pi
	// hbar = 6.582 E-16 eV*c
// ���������� �������, ����� �� ������������ HESR 
	// omega_P = 2.9 E8 c^(-1)	
// �������� ����� � �������
	// � = 2.99792 E10 cm/c
// ������� ����� ���������
	// Ee = m_e*c^2 =  5.11 E5 eV
double delta (double Te_per);

double delta_norm (const double velosity_ion, const double velosity_e,  const double frequency_Eplasma, const double frequency_Ecyclotron);

// ���������: �������� ������������ �������� ����, V_i/V_0
	//		� �������� ����������� ������������ ���, Te
// ���������: �������� ��������� tau_i, �������
// tau_i = tau_0*(V_0/V_i)^2
// tau_0 = 1/(m_e*c^2)*c^2/V_0^2
// �������� ����� � �������
	// � = 2.99792 E10 cm/c
// ������� ����� ���������
	// Ee = m_e*c^2 =  5.11 E5 eV
// ������������� ��������
	// V_0 = 10 E6 cm/c
double tau (double ViV0, double Te);

// ��������� �������� ��������� tau_per, �������
// tau_per = Te_per/me/Ve_per^2
double tau_perpendicular() ;

// ��������� �������� ��������� tau_per, �������
// tau_par = Te_par/me/Ve_per^2
double tau_parallel (double Te_par, double Te_per) ;

// ���������: 3 �������������� ���������, ���� double
// ����������: �������� �������, �������������� � ������ 
	// �������� �.�. ������� ����������� �� �������������� ������... 
double ELIA(double gx, double gy, double gz, double* parmtr);

// ���������: 3 �������������� ���������, ���� double
// ����������: �������� �������, �������������� � ������ 
	// �������� �.�. ������ � ������������� ����������� ����... 
double ELiMA(double gx, double gy, double gz, double* parmtr);

double V_0();

double velocity_ion(double velocity_ion_norm);

double velocity_electron(double temperature);

double average_velocity_square(const double temperature_Eperp, const double temperature_Eparall);

double average_frequency(const double frequency_Eplasma, const double frequency_Ecyclotron);

double frequency_norm(double frequency, const double frequency_Eplasma, const double frequency_Ecyclotron);

double average_tau(double tau_perp, double tau_parall);

double tau_perp_norm(double tau_perp, double tau_parall);

double tau_parallel_norm(double tau_perp, double tau_parall);

double delta_norm(const double temperature_Eperp, const double temperature_Eparall, const double frequency_Eplasma, const double frequency_Ecyclotron_norm);

double ELiMA_H(double gx, double gy, double gz, double* set_parmtr);
#endif
