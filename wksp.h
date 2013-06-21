#ifndef WKSP_H
#define WKSP_H

#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <omp.h>


class WKSP  // wksp = workspace
{
public:
	WKSP();  // initailize, memory allocation, define constant;
	~WKSP(); // memory free

	// setting
	void initial_read_setting(void);
	int N_radial; 
	double kc;
	double h_radial; // kc/N_radial;
	int N_theta;
	double h_theta; // 2pi/N_theta;
	double c_theta; // 1/N_theta;
	double kf;
	int numofthread;
	int N;  // num of layer
	int N2; // dim of system.
	// N2= 2*N for full band model, N2 = 2 for 2band model

	//constant
	void initial_define_constant(void);
	double indc_nlayer_arr[17]={1,2,4,5,6,8,10,15};
	double gamma0;
	double gamma1;
	double gamma2;
	double gamma3;
	double gamma4;
	double gamma5;
	double gamma6;	
	double epsilon; // dielectric constant
	double esq; //e^2

	double dielec;
	double e_mass;
	double imp;
	double c;
	double kB;

	double h;
	double hbar;
	double hev;
	double hevbar;

	double a; //lattice constant 2.46
	double d; // distance between layers 3.35 angstrom
	double dbar; // d/a (unitless)
	double esq_a;  //e^2/a
	double esq_ea; //e^2/a/epsilon
	double alpha;
	int tot_th;

	/// memory
	void initial_malloc(void);
	gsl_eigen_hermv_workspace** ws;
	gsl_vector** eval;

	gsl_matrix_complex*** eigen_state; //matrix* eigenstate[radial][theta]
	gsl_matrix_complex**** epho; //density_matrix //matrix* epho[band][radial][theta]
	double*** energy; // energy[band][radial][theta]
	void set_en(void); // calculate electron density
	double* en; // electron density[layer]

	gsl_matrix_complex*** H; //matrix* Hamailtonian[radial][theta]

	void diag_slfcssnt(int indc,double WFM,double nB,double nT,double Ef);
	void calc_opdc(double eta, double Ef);
	void exec_diag(int indc);

	void set_H(void);

	void set_H_A(void);
	//2band model
	void set_H_AB2(void);
	void set_H_ABC2(void);
	void set_H_ABCA2(void);

	// full band model
	void set_H_AB(void);
	
	void set_H_ABC(void);
	void set_H_ABCA(void);
	void set_H_ABCAB(void);

	void set_H_ABA(void);
	void set_H_ABAB(void);
	void set_H_ABABA(void);
	
	int calcul_pho(void); // wavefunction,eigensate -> density matrix
	void band_cal(void);

	void set_initial(char* log);
	void set_initial_n(char* log);

	double vg; // external field
	void print_epho(void);
	void print_epho44(void);
	void print_H(void);
	void print_en(void);
	void print_ken(void); //electron density ditribution[radial]
	void print_energy(void);

	bool hartree;
	bool fock;
};

#endif
