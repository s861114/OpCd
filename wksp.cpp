#include "wksp.h"
#define PI 3.14159265358979 
#define realK(i) i/N_radial*kc/reala

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

void WKSP::band_cal(void)
{
	FILE *fp=fopen("output.txt","w");
	#pragma omp parallel for
	for(int i=0; i<N_radial; i++)
	{
		int cur_th = omp_get_thread_num();
		// getting current thread num
		for(int j=0; j<N_theta; j++)
		{
			gsl_eigen_hermv(H[i][j],eval[cur_th],eigen_state[i][j],ws[cur_th]);
			gsl_eigen_hermv_sort(eval[cur_th],eigen_state[i][j],GSL_EIGEN_SORT_VAL_ASC);
			
			for(int e=0; e<N2; e++)
				energy[e][i][j] = eval[cur_th]->data[e];			
		}
	}
				
	for(int i=0;i<N_radial;i++)
	{
		fprintf(fp,"%d ",i);
		for(int e=0;e<N2;e++)

		{	
			fprintf(fp,"%f ",energy[e][i][1]);
		}
		fprintf(fp,"\n");
	}


}

void WKSP::frmlvl_skr(void)
{
	double k1;
	double k2;
	int cntclc=-1;
	double drvtE;
	double hE;
	double lE;
	double delta=0;
	double density;
	int bandidx;
	int kcnt;
	int flag;
	for (double cyclcEf=Ef/2;cyclcEf<Ef;cyclcEf+=0.005)
	{
		cntclc++;
		density=0;

		//Finding the density corressponding to a given Ef	
		for (bandidx=2;bandidx==2;bandidx++)
		{
			k1=0;k2=0;delta=0;hE=0;lE=0;flag=0;


			for (kcnt=0;kcnt<N_radial-3;kcnt++)
			{
				hE=energy[bandidx][kcnt+1][0];
				lE=energy[bandidx][kcnt][0];
		//		printf("%d le=%f cyclcEf=%f\n",kcnt,hE,cyclcEf);
		//		getchar();
				
				if( (hE-cyclcEf)*(lE-cyclcEf)<0 )
				{
					drvtE=(hE-lE)/h_radial; //derivative of energy w.r.t k
					switch (flag)
					{
						case 0:
							k1=kcnt+(cyclcEf-hE)/drvtE;
							flag=1;
							break;
						case 1:
							k2=kcnt+(cyclcEf-hE)/drvtE;
							flag=2;
							goto EndOfLoop;
														
					}
				}
			
			}
EndOfLoop:
//			printf("%f\n,kcnt=%d, flag=%d\n",k1,kcnt,flag);
			switch (flag){
				case 0:delta=0;	break;	
				case 1:delta=(realK(k1)*realK(k1) )/PI;break;		
				case 2:delta=(realK(k2)*realK(k2)-realK(k1)*realK(k1) )/PI;	break;
			}

			if(bandidx>N-1)
				density+=delta;
			else
				density-=delta;				
		}
		//end of the part finding the density
		printf("%e\n",density);
		break;		
	}
		

}

void WKSP::slfcssnt(void)
{


}

void WKSP::opdc(void)
{
}
