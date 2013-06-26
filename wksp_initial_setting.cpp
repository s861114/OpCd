#include "wksp.h"
#include "function.h"

WKSP::WKSP()
{

	initial_define_constant();
	initial_read_setting();
}

void WKSP::initial_define_constant(void)
{

	gamma0=3; // 3eV
	gamma1=0.3;//0.3eV
	esq = 14.39966058372; // eV*angstrom
	// e= 4.803205e-10 statC
	// (1statC)^2 = 1 erg cm
	// 1 erg = 6.241509e+11 eV
	a=2.46; // lattice constant 2.46 angstrom
	reala=2.46e-10;
	d=3.35; // distance between layers 3.35 angstrom
	hv_a = 0.5*sqrt(3)*gamma0; // eV
	esq_a=esq/a;  // eV
	alpha0=esq_a/hv_a; // e^2/(hbar v0) epsilon=1  alpha<=alpha0;
	dbar = d/a;
	epsilon=8.85418e-12;
	slf_const=1.602e-19*0.334e-9/(2*epsilon);
	set_alpha(1);

	// tot_th means num of logical CPU.
	//  = physical CPU # * 2 (hyperthreading)
	#pragma omp parallel
	{
		tot_th = omp_get_num_threads();
	}
}
void WKSP::set_alpha(double new_alpha)
{
	alpha = new_alpha;
	epsilon = alpha0/alpha;
	esq_ea = esq_a/epsilon;
}
void WKSP::set_epsilon(double new_epsilon)
{
	epsilon = new_epsilon;
	esq_ea = esq_a/epsilon;
	alpha = alpha0/epsilon;
}
void WKSP::initial_read_setting(void)
{
	FILE* fp = fopen("setup.set","r");
	char dump[100];
	char dump1[100];
	/*while(!feof(fp))
	{
		fgets(dump,100,fp);
		printf("%s",dump);
	}*/
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_radial);
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&kc);		
	fgets(dump,100,fp);	sscanf(dump,"%s  %d\n",dump1,&N_theta);	
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&kf);		
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&alpha);
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&Ef);
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&nB);
	fgets(dump,100,fp);	sscanf(dump,"%s  %lf\n",dump1,&nT);
	nB=nB*1e+16;
	nT=nT*1e+16;

	set_alpha(alpha);
	h_radial = kc / N_radial;
	h_theta = 2.0*M_PI/ N_theta;
	c_theta = 1.0/N_theta;	

}

void WKSP::initial_malloc(int Nnew, int N2new)
{
	printf("start initial_malloc()\n");
	umat<gsl_eigen_hermv_workspace*> mat_gsl_eigen_hermv_workspace_pointer;
	umat<gsl_matrix_complex*> mat_gsl_matrix_complex_pointer;
	umat<gsl_vector*> mat_gsl_vector_pointer;
	umat<double> mat_double;


	ws = mat_gsl_eigen_hermv_workspace_pointer.Tmatrix1(tot_th);
	eval = mat_gsl_vector_pointer.Tmatrix1(tot_th);
	for(int i=0; i<tot_th; i++)
	{
		ws[i] = gsl_eigen_hermv_alloc(N2new);
		eval[i] = gsl_vector_alloc(N2new);
	}

	eigen_state = mat_gsl_matrix_complex_pointer.Tmatrix2(N_radial,N_theta);
	H = mat_gsl_matrix_complex_pointer.Tmatrix2(N_radial,N_theta);
	epho = mat_gsl_matrix_complex_pointer.Tmatrix3(N2new,N_radial,N_theta);
	// N2 is # of bnads of system.
	energy = mat_double.Tmatrix3(N2new,N_radial,N_theta);
	en = mat_double.Tmatrix1(N2new);
	diag_term=mat_double.Tmatrix1(Nnew);

	for(int ii=0;ii<Nnew;ii++)
		diag_term[ii]=0;

	for(int i=0; i<N_radial; i++)
	{
		for(int j=0; j<N_theta; j++)
		{
			eigen_state[i][j]= gsl_matrix_complex_alloc(N2new,N2new);
			H[i][j]= gsl_matrix_complex_alloc(N2new,N2new);
			for(int e=0; e<N2new; e++)
			{
				epho[e][i][j] = gsl_matrix_complex_alloc(N2new,N2new);
			}
		}
	//	printf("malloc %d\n",i);
	}
	printf("finish initial_malloc()\n");
}


WKSP::~WKSP()
{
	printf("start memory free\n");
	umat<gsl_eigen_hermv_workspace*> mat_gsl_eigen_hermv_workspace_pointer;
	umat<gsl_matrix_complex*> mat_gsl_matrix_complex_pointer;
	umat<gsl_vector*> mat_gsl_vector_pointer;
	umat<double> mat_double;
	

	for(int i=0; i<tot_th; i++)
	{
		gsl_eigen_hermv_free(ws[i]);
		gsl_vector_free(eval[i]);
	}

	for(int i=0; i<N_radial; i++)
	{
		for(int j=0; j<N_theta; j++)
		{
			gsl_matrix_complex_free(eigen_state[i][j]);
			gsl_matrix_complex_free(H[i][j]);
			for(int e=0; e<N2; e++)
			{
				gsl_matrix_complex_free(epho[e][i][j]);
			}
		}
	}
	mat_gsl_matrix_complex_pointer.Tfree2(eigen_state);
	mat_gsl_matrix_complex_pointer.Tfree2(H);
	mat_gsl_matrix_complex_pointer.Tfree3(epho);
	mat_double.Tfree3(energy);
	mat_double.Tfree1(en);
	mat_double.Tfree1(diag_term);

	printf("finish memory free\n");
}


