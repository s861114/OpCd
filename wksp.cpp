#include "wksp.h"
#define PI 3.14159265358979 
#define realK(i) i/N_radial*kc/reala

int WKSP::calcul_pho(int N, int N2)
{
/*	#pragma omp parallel for
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
	return 0;*/
}

void WKSP::band_cal(int N,int N2)
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
				
/*	for(int i=0;i<N_radial;i++)
	{
		fprintf(fp,"%d ",i);
		for(int e=0;e<N2;e++)
		{	
			fprintf(fp,"%f ",energy[e][i][0]);
		}
		fprintf(fp,"\n");
	}
*/

}


//Searching for the fermi level with a charge carrier diesnty given
double WKSP::frmlvl_skr(int N, int N2,double density_tar)
{
	double k1;
	double k2;
	int cntclc=-1;
	double drvtE;
	double hE;
	double lE;
	double delta=0;
	double density;
	double density_old=density_tar;
	int bandidx;
	int kcnt;
	int flag;
	for (double cyclcEf=-Ef;cyclcEf<Ef;cyclcEf+=0.0001)
	{
		cntclc++;
		density=0;
		//Finding the density corressponding to a given Ef
		for (bandidx=0;bandidx<N2;bandidx++)
		{
			k1=0;k2=0;delta=0;hE=0;lE=0;flag=0;


			for (kcnt=0;kcnt<N_radial-3;kcnt++)
			{
				hE=energy[bandidx][kcnt+1][0];
				lE=energy[bandidx][kcnt][0];
				
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
//		printf("density %e | density_tar %e | density_old %e",density,density_tar,density_old);
//		getchar();
		if ((density-density_tar)*(density_old-density_tar)<0)
		{
			printf("real value : %f\n",cyclcEf);
			return cyclcEf;
		break;
		}
		density_old=density;
	}
	printf("out of range");
	getchar();
	
}

//Self-consistent Calculation
void WKSP::slfcssnt_Ef(int N, int N2)
{
	double nT,nB,sum1,sum2;
	double enOld[10]={0,0,0,0,0,0,0,0,0,0};
	int cnt=0;
	double alpha=0.01;
	double chksum=0;
	while(1)
	{
		cnt++;
		set_H_AB();
		band_cal(N,N2);
		find_n(N,N2);
		if (cnt>1){
				for(int i=0;i<N;i++)
					en[i]=en[i]*alpha+enOld[i]*(1-alpha);}

		nT=0;
		nB=-sum_en(N,N2);
		for(int i=0;i<N;i++)
			diag_term[i]=0;
		
				
		for(int i=0;i<N-1;i++)
		{
			sum1=0;
			sum2=0;
			for(int j=i+1;j<N;j++)
				sum1+=en[j];
			for(int j=0;j<=i;j++)
				sum2+=en[j];
		
			diag_term[i+1]=diag_term[i]+slf_const*((nT-nB)+sum1-sum2);
		}
		printf("\ndiag_term0 : %e, diag_term1 : %e\n",diag_term[0],diag_term[1]);
		
		chksum=0;
		for(int i=0;i<N;i++)
			chksum=chksum+en[i]-enOld[i];
//		if (chksum<1e+10)
//			break;



		for(int i=0;i<N;i++)
			enOld[i]=en[i];
		
		printf("\n%d",cnt);
		for(int i=0;i<N;i++)
		{
			printf("	%e",cnt,en[i]);
		}
	}	

}

double WKSP::sum_en(int N, int N2)
{
	double sum=0;
	for(int i=0;i<N;i++)
		sum=sum+en[i];
	return sum;
}

void WKSP::slfcssnt_density(int N, int N2)
{
}	

void WKSP::find_n(int N, int N2)
{
	double sum;
	double temp;
	double kRes=(kc/reala)/N_radial;
	for(int i=0;i<N;i++)
		en[i]=0;

	for (int kidx=0;kidx<N_radial;kidx++)
	{
		for (int i=0;i<N;i++) //here i refers to a layer
		{
			for(int j=0;j<N2;j++) //here j is a band index
			{
				if(energy[j][kidx][0]<=Ef)
				{
					sum=0;	
					for(int sumidx=2*i;sumidx<=2*i+1;sumidx++)
					{						
						temp=gsl_complex_abs2(gsl_matrix_complex_get(eigen_state[kidx][0],sumidx,j));
						sum=sum+2.0/PI*realK(kidx)*temp*kRes;
						//sum=sum+temp*4.0/(2*PI)/(2*PI)*h_radial;
					}
					en[i]=en[i]+sum;
				}	

			}
		}

		for (int j=0;j<N2/2;j++)
			en[j]=en[j]-2/PI*realK(kidx)*kRes;

	}
			


}
void WKSP::opdc(int N, int N2)
{
}
