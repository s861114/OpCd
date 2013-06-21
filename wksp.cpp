#include "wksp.h"

int WKSP::calcul_pho(void)
{
	#pragma omp parallel for
	for(int i=0; i<N_radial; i++)
	{
		int cur_th = omp_get_thread_num();
		// getting current thread num
		for(int j=0; j<N_theta; j++)
		{
			gsl_eigen_hermv(H[i][j],eval[cur_th],eigen_state[i][j],ws[cur_th]);			
			gsl_eigen_hermv_sort(eval[cur_th],eigen_state[i][j],GSL_EIGEN_SORT_VAL_ASC);
			// 0 is lowest energy band index
			for(int e=0; e<N2; e++)
				energy[e][i][j] = eval[cur_th]->data[e];

			for(int p=0; p<N2; p++)
			{
				for(int q=0; q<N2; q++)
				{
					for(int e=0; e<N2; e++)
					{
						gsl_complex newp = gsl_complex_conjugate(gsl_matrix_complex_get(eigen_state[i][j],p,e));
						gsl_complex newq = gsl_matrix_complex_get(eigen_state[i][j],q,e);

						gsl_complex newpq = gsl_complex_mul(newp,newq);

						gsl_matrix_complex_set(epho[e][i][j],p,q,newpq);
					}
				}
			}
		}
	}
	return 0;
}

void WKSP::diag_slfcssnt(int indc, double WFM, double nB, double nT, double Ef)
{
		exec_diag(indc);	

}



void WKSP::calc_opdc(double eta,double Ef)
{
}

void WKSP::exec_diag(int indc)
{
	switch (indc)
	{
		case 1:set_H_A();break;
		case 2:set_H_AB();break;
/*		case 4:set_H_ABC();break;
		case 5:set_H_ABA();break;
		case 6:set_H_ABCA();break;
		case 8:set_H_ABAB();break;
		case 10:set_H_ABCAB();break;
		case 15:set_H_ABABA();break;*/
	}
}
