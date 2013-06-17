#include "wksp.h"

void WKSP::set_H_A(void)
{
	#pragma omp parallel for
	for(int i=0; i<N_radial; i++)
	{
		double kk=(i+0.5)*h_radial;
		for(int j=0; j<N_theta; j++)
		{
			double theta = (j*0.5)*h_theta;
			double aa = hv_a;

			gsl_complex offdia01 = gsl_complex_polar(aa*kk,-theta);// = var1 exp(i var2)
			gsl_complex offdia10 = gsl_complex_conjugate(offdia01);

			gsl_matrix_complex_set_zero(H[i][j]);
			gsl_matrix_complex_set(H[i][j],0,0, gsl_complex_rect(vg,0) );
			gsl_matrix_complex_set(H[i][j],0,1,  offdia01);
			gsl_matrix_complex_set(H[i][j],1,0,  offdia10);
			gsl_matrix_complex_set(H[i][j],1,1, gsl_complex_rect(vg,0) );
		}
	}
}

void WKSP::set_H_AB(void)
{
	#pragma omp parallel for
	for(int i=0; i<N_radial; i++)
	{
		double kk=(i+0.5)*h_radial;
		for(int j=0; j<N_theta; j++)
		{
			double theta = (j*0.5)*h_theta;
			double aa = hv_a;

			gsl_complex offdia01 = gsl_complex_polar(aa*kk,-theta);// = var1 exp(i var2)
			gsl_complex offdia10 = gsl_complex_conjugate(offdia01);

			gsl_matrix_complex_set_zero(H[i][j]);
			gsl_matrix_complex_set(H[i][j],0,0, gsl_complex_rect(vg,0) );
			gsl_matrix_complex_set(H[i][j],0,1,  offdia01);
			gsl_matrix_complex_set(H[i][j],1,0,  offdia10);
			gsl_matrix_complex_set(H[i][j],1,1, gsl_complex_rect(vg,0) );
			gsl_matrix_complex_set(H[i][j],1,2, gsl_complex_rect(gamma1,0) );
			gsl_matrix_complex_set(H[i][j],2,1, gsl_complex_rect(gamma1,0) );
			gsl_matrix_complex_set(H[i][j],2,2, gsl_complex_rect(vg,0) );
			gsl_matrix_complex_set(H[i][j],2,3,  offdia01);
			gsl_matrix_complex_set(H[i][j],3,2,  offdia10);
			gsl_matrix_complex_set(H[i][j],3,3, gsl_complex_rect(vg,0) );
		}
	}
}
