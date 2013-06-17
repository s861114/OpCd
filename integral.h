#include <stdio.h>
#include <math.h>
#include <omp.h>
double gq( double (*fun)(double),double a, double b,int N);
double adaptive(double (*fun)(double), double a, double b);
double adaptive_heap(double (*fun)(double), double a, double b);
double adaptive_internal(double (*fun)(double),double a,double b
		,int step,double torr);
double gq2d( double (*fun)(double,double),double xa, double xb,double ya,double yb,int N);
double adaptive2d( double (*fun)(double,double),double xa, double xb,double ya,double yb);
double adaptive2d_internal( double (*fun)(double,double)
		,double xa, double xb,double ya,double yb
		,int step, double torr);

double gq2d_polar( double (*fun)(double,double,double)
		,double qq,double xa, double xb,double ya,double yb,int N);
double adaptive2d_polar( double (*fun)(double,double,double)
		,double qq,double xa, double xb,double ya,double yb);
double adaptive2d_polar_internal( double (*fun)(double,double,double)
		,double qq,double xa, double xb,double ya,double yb
		,int step, double torr);

double gq( double (*fun)(double),double a, double b,int N)
{
	double h = (b-a)/double(2*N);
	double total=0;
	for(int i=0; i<N; i++)
	{
		double x1 = a+(i*2+1+0.0000000000000000)*h;
		double x2 = a+(i*2+1-0.0586850543002595)*h;
		double x3 = a+(i*2+1+0.0586850543002595)*h;
		double x4 = a+(i*2+1-0.1171678090719551)*h;
		double x5 = a+(i*2+1+0.1171678090719551)*h;
		double x6 = a+(i*2+1-0.1752466621553257)*h;
		double x7 = a+(i*2+1+0.1752466621553257)*h;
		double x8 = a+(i*2+1-0.2327214037242726)*h;
		double x9 = a+(i*2+1+0.2327214037242726)*h;
		double x10= a+(i*2+1-0.2893939064516262)*h;
		double x11= a+(i*2+1+0.2893939064516262)*h;
		double x12= a+(i*2+1-0.3450688084957224)*h;
		double x13= a+(i*2+1+0.3450688084957224)*h;
		double x14= a+(i*2+1-0.3995541869539530)*h;
		double x15= a+(i*2+1+0.3995541869539530)*h;
		double x16= a+(i*2+1-0.4526622194618458)*h;
		double x17= a+(i*2+1+0.4526622194618458)*h;
		double x18= a+(i*2+1-0.5042098316571334)*h;
		double x19= a+(i*2+1+0.5042098316571334)*h;
		double x20= a+(i*2+1-0.5540193282770679)*h;
		double x21= a+(i*2+1+0.5540193282770679)*h;
		double x22= a+(i*2+1-0.6019190057137693)*h;
		double x23= a+(i*2+1+0.6019190057137693)*h;
		double x24= a+(i*2+1-0.6477437439165100)*h;
		double x25= a+(i*2+1+0.6477437439165100)*h;
		double x26= a+(i*2+1-0.6913355756013667)*h;
		double x27= a+(i*2+1+0.6913355756013667)*h;
		double x28= a+(i*2+1-0.7325442308075103)*h;
		double x29= a+(i*2+1+0.7325442308075103)*h;
		double x30= a+(i*2+1-0.7712276549255324)*h;
		double x31= a+(i*2+1+0.7712276549255324)*h;
		double x32= a+(i*2+1-0.8072524984168955)*h;
		double x33= a+(i*2+1+0.8072524984168955)*h;
		double x34= a+(i*2+1-0.8404945765458014)*h;
		double x35= a+(i*2+1+0.8404945765458014)*h;
		double x36= a+(i*2+1-0.8708392975582413)*h;
		double x37= a+(i*2+1+0.8708392975582413)*h;
		double x38= a+(i*2+1-0.8981820578754266)*h;
		double x39= a+(i*2+1+0.8981820578754266)*h;
		double x40= a+(i*2+1-0.9224286030428122)*h;
		double x41= a+(i*2+1+0.9224286030428122)*h;
		double x42= a+(i*2+1-0.9434953534644419)*h;
		double x43= a+(i*2+1+0.9434953534644419)*h;
		double x44= a+(i*2+1-0.9613096946231363)*h;
		double x45= a+(i*2+1+0.9613096946231363)*h;
		double x46= a+(i*2+1-0.9758102337149845)*h;
		double x47= a+(i*2+1+0.9758102337149845)*h;
		double x48= a+(i*2+1-0.9869470350233716)*h;
		double x49= a+(i*2+1+0.9869470350233716)*h;
		double x50= a+(i*2+1-0.9946819193080071)*h;
		double x51= a+(i*2+1+0.9946819193080071)*h;
		double x52= a+(i*2+1-0.9989899477763282)*h;
		double x53= a+(i*2+1+0.9989899477763282)*h;


		total+=0.0587187941511644*fun(x1);
		total+=0.0586175862327203*fun(x2);
		total+=0.0586175862327203*fun(x3);
		total+=0.0583143113622560*fun(x4);
		total+=0.0583143113622560*fun(x5);
		total+=0.0578100149917132*fun(x6);
		total+=0.0578100149917132*fun(x7);
		total+=0.0571064355362672*fun(x8);
		total+=0.0571064355362672*fun(x9);
		total+=0.0562059983817397*fun(x10);
		total+=0.0562059983817397*fun(x11);
		total+=0.0551118075239336*fun(x12);
		total+=0.0551118075239336*fun(x13);
		total+=0.0538276348687310*fun(x14);
		total+=0.0538276348687310*fun(x15);
		total+=0.0523579072298727*fun(x16);
		total+=0.0523579072298727*fun(x17);
		total+=0.0507076910692927*fun(x18);
		total+=0.0507076910692927*fun(x19);
		total+=0.0488826750326991*fun(x20);
		total+=0.0488826750326991*fun(x21);
		total+=0.0468891503407503*fun(x22);
		total+=0.0468891503407503*fun(x23);
		total+=0.0447339891036728*fun(x24);
		total+=0.0447339891036728*fun(x25);
		total+=0.0424246206345200*fun(x26);
		total+=0.0424246206345200*fun(x27);
		total+=0.0399690058435404*fun(x28);
		total+=0.0399690058435404*fun(x29);
		total+=0.0373756098034829*fun(x30);
		total+=0.0373756098034829*fun(x31);
		total+=0.0346533725835342*fun(x32);
		total+=0.0346533725835342*fun(x33);
		total+=0.0318116784590193*fun(x34);
		total+=0.0318116784590193*fun(x35);
		total+=0.0288603236178237*fun(x36);
		total+=0.0288603236178237*fun(x37);
		total+=0.0258094825107575*fun(x38);
		total+=0.0258094825107575*fun(x39);
		total+=0.0226696730570702*fun(x40);
		total+=0.0226696730570702*fun(x41);
		total+=0.0194517211076369*fun(x42);
		total+=0.0194517211076369*fun(x43);
		total+=0.0161667252566875*fun(x44);
		total+=0.0161667252566875*fun(x45);
		total+=0.0128260261442404*fun(x46);
		total+=0.0128260261442404*fun(x47);
		total+=0.0094412022849403*fun(x48);
		total+=0.0094412022849403*fun(x49);
		total+=0.0060242762269487*fun(x50);
		total+=0.0060242762269487*fun(x51);
		total+=0.0025916837205670*fun(x52);
		total+=0.0025916837205670*fun(x53);
	}
	total*=h;
	return total;
}

double gq2d( double (*fun)(double,double),double xa, double xb,double ya,double yb,int N)
{
	double hx = (xb-xa)/double(2*N);
	double hy = (yb-ya)/double(2*N);
	double total=0;
	double weight[54];
	weight[1 ]=0.0587187941511644;
	weight[2 ]=0.0586175862327203;
	weight[3 ]=0.0586175862327203;
	weight[4 ]=0.0583143113622560;
	weight[5 ]=0.0583143113622560;
	weight[6 ]=0.0578100149917132;
	weight[7 ]=0.0578100149917132;
	weight[8 ]=0.0571064355362672;
	weight[9 ]=0.0571064355362672;
	weight[10]=0.0562059983817397;
	weight[11]=0.0562059983817397;
	weight[12]=0.0551118075239336;
	weight[13]=0.0551118075239336;
	weight[14]=0.0538276348687310;
	weight[15]=0.0538276348687310;
	weight[16]=0.0523579072298727;
	weight[17]=0.0523579072298727;
	weight[18]=0.0507076910692927;
	weight[19]=0.0507076910692927;
	weight[20]=0.0488826750326991;
	weight[21]=0.0488826750326991;
	weight[22]=0.0468891503407503;
	weight[23]=0.0468891503407503;
	weight[24]=0.0447339891036728;
	weight[25]=0.0447339891036728;
	weight[26]=0.0424246206345200;
	weight[27]=0.0424246206345200;
	weight[28]=0.0399690058435404;
	weight[29]=0.0399690058435404;
	weight[30]=0.0373756098034829;
	weight[31]=0.0373756098034829;
	weight[32]=0.0346533725835342;
	weight[33]=0.0346533725835342;
	weight[34]=0.0318116784590193;
	weight[35]=0.0318116784590193;
	weight[36]=0.0288603236178237;
	weight[37]=0.0288603236178237;
	weight[38]=0.0258094825107575;
	weight[39]=0.0258094825107575;
	weight[40]=0.0226696730570702;
	weight[41]=0.0226696730570702;
	weight[42]=0.0194517211076369;
	weight[43]=0.0194517211076369;
	weight[44]=0.0161667252566875;
	weight[45]=0.0161667252566875;
	weight[46]=0.0128260261442404;
	weight[47]=0.0128260261442404;
	weight[48]=0.0094412022849403;
	weight[49]=0.0094412022849403;
	weight[50]=0.0060242762269487;
	weight[51]=0.0060242762269487;
	weight[52]=0.0025916837205670;
	weight[53]=0.0025916837205670;
	for(int i=0; i<N; i++)
	{
		double x[54];
		x[1 ]= xa+(i*2+1+0.0000000000000000)*hx;
		x[2 ]= xa+(i*2+1-0.0586850543002595)*hx;
		x[3 ]= xa+(i*2+1+0.0586850543002595)*hx;
		x[4 ]= xa+(i*2+1-0.1171678090719551)*hx;
		x[5 ]= xa+(i*2+1+0.1171678090719551)*hx;
		x[6 ]= xa+(i*2+1-0.1752466621553257)*hx;
		x[7 ]= xa+(i*2+1+0.1752466621553257)*hx;
		x[8 ]= xa+(i*2+1-0.2327214037242726)*hx;
		x[9 ]= xa+(i*2+1+0.2327214037242726)*hx;
		x[10]= xa+(i*2+1-0.2893939064516262)*hx;
		x[11]= xa+(i*2+1+0.2893939064516262)*hx;
		x[12]= xa+(i*2+1-0.3450688084957224)*hx;
		x[13]= xa+(i*2+1+0.3450688084957224)*hx;
		x[14]= xa+(i*2+1-0.3995541869539530)*hx;
		x[15]= xa+(i*2+1+0.3995541869539530)*hx;
		x[16]= xa+(i*2+1-0.4526622194618458)*hx;
		x[17]= xa+(i*2+1+0.4526622194618458)*hx;
		x[18]= xa+(i*2+1-0.5042098316571334)*hx;
		x[19]= xa+(i*2+1+0.5042098316571334)*hx;
		x[20]= xa+(i*2+1-0.5540193282770679)*hx;
		x[21]= xa+(i*2+1+0.5540193282770679)*hx;
		x[22]= xa+(i*2+1-0.6019190057137693)*hx;
		x[23]= xa+(i*2+1+0.6019190057137693)*hx;
		x[24]= xa+(i*2+1-0.6477437439165100)*hx;
		x[25]= xa+(i*2+1+0.6477437439165100)*hx;
		x[26]= xa+(i*2+1-0.6913355756013667)*hx;
		x[27]= xa+(i*2+1+0.6913355756013667)*hx;
		x[28]= xa+(i*2+1-0.7325442308075103)*hx;
		x[29]= xa+(i*2+1+0.7325442308075103)*hx;
		x[30]= xa+(i*2+1-0.7712276549255324)*hx;
		x[31]= xa+(i*2+1+0.7712276549255324)*hx;
		x[32]= xa+(i*2+1-0.8072524984168955)*hx;
		x[33]= xa+(i*2+1+0.8072524984168955)*hx;
		x[34]= xa+(i*2+1-0.8404945765458014)*hx;
		x[35]= xa+(i*2+1+0.8404945765458014)*hx;
		x[36]= xa+(i*2+1-0.8708392975582413)*hx;
		x[37]= xa+(i*2+1+0.8708392975582413)*hx;
		x[38]= xa+(i*2+1-0.8981820578754266)*hx;
		x[39]= xa+(i*2+1+0.8981820578754266)*hx;
		x[40]= xa+(i*2+1-0.9224286030428122)*hx;
		x[41]= xa+(i*2+1+0.9224286030428122)*hx;
		x[42]= xa+(i*2+1-0.9434953534644419)*hx;
		x[43]= xa+(i*2+1+0.9434953534644419)*hx;
		x[44]= xa+(i*2+1-0.9613096946231363)*hx;
		x[45]= xa+(i*2+1+0.9613096946231363)*hx;
		x[46]= xa+(i*2+1-0.9758102337149845)*hx;
		x[47]= xa+(i*2+1+0.9758102337149845)*hx;
		x[48]= xa+(i*2+1-0.9869470350233716)*hx;
		x[49]= xa+(i*2+1+0.9869470350233716)*hx;
		x[50]= xa+(i*2+1-0.9946819193080071)*hx;
		x[51]= xa+(i*2+1+0.9946819193080071)*hx;
		x[52]= xa+(i*2+1-0.9989899477763282)*hx;
		x[53]= xa+(i*2+1+0.9989899477763282)*hx;
		           
		for(int j=0; j<N; j++)
		{
			double y[54];
			y[1 ]= ya+(j*2+1+0.0000000000000000)*hy;
			y[2 ]= ya+(j*2+1-0.0586850543002595)*hy;
			y[3 ]= ya+(j*2+1+0.0586850543002595)*hy;
			y[4 ]= ya+(j*2+1-0.1171678090719551)*hy;
			y[5 ]= ya+(j*2+1+0.1171678090719551)*hy;
			y[6 ]= ya+(j*2+1-0.1752466621553257)*hy;
			y[7 ]= ya+(j*2+1+0.1752466621553257)*hy;
			y[8 ]= ya+(j*2+1-0.2327214037242726)*hy;
			y[9 ]= ya+(j*2+1+0.2327214037242726)*hy;
			y[10]= ya+(j*2+1-0.2893939064516262)*hy;
			y[11]= ya+(j*2+1+0.2893939064516262)*hy;
			y[12]= ya+(j*2+1-0.3450688084957224)*hy;
			y[13]= ya+(j*2+1+0.3450688084957224)*hy;
			y[14]= ya+(j*2+1-0.3995541869539530)*hy;
			y[15]= ya+(j*2+1+0.3995541869539530)*hy;
			y[16]= ya+(j*2+1-0.4526622194618458)*hy;
			y[17]= ya+(j*2+1+0.4526622194618458)*hy;
			y[18]= ya+(j*2+1-0.5042098316571334)*hy;
			y[19]= ya+(j*2+1+0.5042098316571334)*hy;
			y[20]= ya+(j*2+1-0.5540193282770679)*hy;
			y[21]= ya+(j*2+1+0.5540193282770679)*hy;
			y[22]= ya+(j*2+1-0.6019190057137693)*hy;
			y[23]= ya+(j*2+1+0.6019190057137693)*hy;
			y[24]= ya+(j*2+1-0.6477437439165100)*hy;
			y[25]= ya+(j*2+1+0.6477437439165100)*hy;
			y[26]= ya+(j*2+1-0.6913355756013667)*hy;
			y[27]= ya+(j*2+1+0.6913355756013667)*hy;
			y[28]= ya+(j*2+1-0.7325442308075103)*hy;
			y[29]= ya+(j*2+1+0.7325442308075103)*hy;
			y[30]= ya+(j*2+1-0.7712276549255324)*hy;
			y[31]= ya+(j*2+1+0.7712276549255324)*hy;
			y[32]= ya+(j*2+1-0.8072524984168955)*hy;
			y[33]= ya+(j*2+1+0.8072524984168955)*hy;
			y[34]= ya+(j*2+1-0.8404945765458014)*hy;
			y[35]= ya+(j*2+1+0.8404945765458014)*hy;
			y[36]= ya+(j*2+1-0.8708392975582413)*hy;
			y[37]= ya+(j*2+1+0.8708392975582413)*hy;
			y[38]= ya+(j*2+1-0.8981820578754266)*hy;
			y[39]= ya+(j*2+1+0.8981820578754266)*hy;
			y[40]= ya+(j*2+1-0.9224286030428122)*hy;
			y[41]= ya+(j*2+1+0.9224286030428122)*hy;
			y[42]= ya+(j*2+1-0.9434953534644419)*hy;
			y[43]= ya+(j*2+1+0.9434953534644419)*hy;
			y[44]= ya+(j*2+1-0.9613096946231363)*hy;
			y[45]= ya+(j*2+1+0.9613096946231363)*hy;
			y[46]= ya+(j*2+1-0.9758102337149845)*hy;
			y[47]= ya+(j*2+1+0.9758102337149845)*hy;
			y[48]= ya+(j*2+1-0.9869470350233716)*hy;
			y[49]= ya+(j*2+1+0.9869470350233716)*hy;
			y[50]= ya+(j*2+1-0.9946819193080071)*hy;
			y[51]= ya+(j*2+1+0.9946819193080071)*hy;
			y[52]= ya+(j*2+1-0.9989899477763282)*hy;
			y[53]= ya+(j*2+1+0.9989899477763282)*hy;
	
			for(int p=1; p<54; p++)
				for(int q=1; q<54; q++)
					total+=weight[p]*weight[q]*fun(x[p],y[q]);
		}
	}
	total*=hx*hy;
	return total;
}
double adaptive(double (*fun)(double), double a, double b)
{
	double torr = 1e-17;
	int th=32;
	double total=0;
	double h = (b-a)/th;
	for(int i=0; i<th; i++)
	{
		double init = a + h*i;
		double end = a+h*(i+1);
		total+=adaptive_internal(fun,init,end,0,torr);
	}
	return total;
}
double adaptive_internal(double (*fun)(double),double a,double b,int step,double torr)
{
	double res1 = gq(fun,a,b,4);
	double res2 = gq(fun,a,b,1);
	double mid = 0.5*(a+b);
//	printf("%e~\t%e\t\t%e\t%e\t%e\t%d\n",a,b,res1,res2,res1-res2,step);
	if( fabs(res1-res2)<torr || step>50 )
		return res1;
	else
		return adaptive_internal(fun,a,mid,step+1,torr)+adaptive_internal(fun,mid,b,step+1,torr);
}
double adaptive2d(double (*fun)(double,double),double xa,double xb,double ya,double yb)
{
	double torr = 1e-17;
	int th=32;
	double total=0;
	double hx = (xb-xa)/th;
	double hy = (yb-ya)/th;
	#pragma omp parallel for reduction (+:total)
	for(int index=0; index<th*th; index++)
	{
		int i = index/th;
		int j = index%th;
		double initx = xa + hx*i;
		double endx = xa+hx*(i+1);
		double inity = ya + hy*j;
		double endy = ya+hy*(j+1);
		total+=adaptive2d_internal(fun,initx,endx,inity,endy,0,torr);
	}
	return total;
}
double adaptive2d_internal(double (*fun)(double,double),double xa,double xb,double ya,double yb,int step,double torr)
{
	double res1 = gq2d(fun,xa,xb,ya,yb,4);
	double res2 = gq2d(fun,xa,xb,ya,yb,1);
	double midx = 0.5*(xa+xb);
	double midy = 0.5*(ya+yb);
//	printf("%f~\t%f\t\t%f\t%f\t%f\n",a,b,res1,res2,res1-res2);
	if( fabs(res1-res2)>torr  && step<50 )
		return adaptive2d_internal(fun,xa,midx,ya,midy,step+1,torr)
			+adaptive2d_internal(fun,midx,xb,ya,midy,step+1,torr)
			+adaptive2d_internal(fun,xa,midx,midy,yb,step+1,torr)
			+adaptive2d_internal(fun,midx,xb,midy,yb,step+1,torr);
	else
		return res1;
}
double gq2d_polar( double (*fun)(double,double,double)
		,double qq,double xa, double xb,double ya,double yb,int N)
{
	double hx = (xb-xa)/double(2*N);
	double hy = (yb-ya)/double(2*N);
	double total=0;
	double weight[54];
	weight[1 ]=0.0587187941511644;
	weight[2 ]=0.0586175862327203;
	weight[3 ]=0.0586175862327203;
	weight[4 ]=0.0583143113622560;
	weight[5 ]=0.0583143113622560;
	weight[6 ]=0.0578100149917132;
	weight[7 ]=0.0578100149917132;
	weight[8 ]=0.0571064355362672;
	weight[9 ]=0.0571064355362672;
	weight[10]=0.0562059983817397;
	weight[11]=0.0562059983817397;
	weight[12]=0.0551118075239336;
	weight[13]=0.0551118075239336;
	weight[14]=0.0538276348687310;
	weight[15]=0.0538276348687310;
	weight[16]=0.0523579072298727;
	weight[17]=0.0523579072298727;
	weight[18]=0.0507076910692927;
	weight[19]=0.0507076910692927;
	weight[20]=0.0488826750326991;
	weight[21]=0.0488826750326991;
	weight[22]=0.0468891503407503;
	weight[23]=0.0468891503407503;
	weight[24]=0.0447339891036728;
	weight[25]=0.0447339891036728;
	weight[26]=0.0424246206345200;
	weight[27]=0.0424246206345200;
	weight[28]=0.0399690058435404;
	weight[29]=0.0399690058435404;
	weight[30]=0.0373756098034829;
	weight[31]=0.0373756098034829;
	weight[32]=0.0346533725835342;
	weight[33]=0.0346533725835342;
	weight[34]=0.0318116784590193;
	weight[35]=0.0318116784590193;
	weight[36]=0.0288603236178237;
	weight[37]=0.0288603236178237;
	weight[38]=0.0258094825107575;
	weight[39]=0.0258094825107575;
	weight[40]=0.0226696730570702;
	weight[41]=0.0226696730570702;
	weight[42]=0.0194517211076369;
	weight[43]=0.0194517211076369;
	weight[44]=0.0161667252566875;
	weight[45]=0.0161667252566875;
	weight[46]=0.0128260261442404;
	weight[47]=0.0128260261442404;
	weight[48]=0.0094412022849403;
	weight[49]=0.0094412022849403;
	weight[50]=0.0060242762269487;
	weight[51]=0.0060242762269487;
	weight[52]=0.0025916837205670;
	weight[53]=0.0025916837205670;
	for(int i=0; i<N; i++)
	{
		double x[54];
		x[1 ]= xa+(i*2+1+0.0000000000000000)*hx;
		x[2 ]= xa+(i*2+1-0.0586850543002595)*hx;
		x[3 ]= xa+(i*2+1+0.0586850543002595)*hx;
		x[4 ]= xa+(i*2+1-0.1171678090719551)*hx;
		x[5 ]= xa+(i*2+1+0.1171678090719551)*hx;
		x[6 ]= xa+(i*2+1-0.1752466621553257)*hx;
		x[7 ]= xa+(i*2+1+0.1752466621553257)*hx;
		x[8 ]= xa+(i*2+1-0.2327214037242726)*hx;
		x[9 ]= xa+(i*2+1+0.2327214037242726)*hx;
		x[10]= xa+(i*2+1-0.2893939064516262)*hx;
		x[11]= xa+(i*2+1+0.2893939064516262)*hx;
		x[12]= xa+(i*2+1-0.3450688084957224)*hx;
		x[13]= xa+(i*2+1+0.3450688084957224)*hx;
		x[14]= xa+(i*2+1-0.3995541869539530)*hx;
		x[15]= xa+(i*2+1+0.3995541869539530)*hx;
		x[16]= xa+(i*2+1-0.4526622194618458)*hx;
		x[17]= xa+(i*2+1+0.4526622194618458)*hx;
		x[18]= xa+(i*2+1-0.5042098316571334)*hx;
		x[19]= xa+(i*2+1+0.5042098316571334)*hx;
		x[20]= xa+(i*2+1-0.5540193282770679)*hx;
		x[21]= xa+(i*2+1+0.5540193282770679)*hx;
		x[22]= xa+(i*2+1-0.6019190057137693)*hx;
		x[23]= xa+(i*2+1+0.6019190057137693)*hx;
		x[24]= xa+(i*2+1-0.6477437439165100)*hx;
		x[25]= xa+(i*2+1+0.6477437439165100)*hx;
		x[26]= xa+(i*2+1-0.6913355756013667)*hx;
		x[27]= xa+(i*2+1+0.6913355756013667)*hx;
		x[28]= xa+(i*2+1-0.7325442308075103)*hx;
		x[29]= xa+(i*2+1+0.7325442308075103)*hx;
		x[30]= xa+(i*2+1-0.7712276549255324)*hx;
		x[31]= xa+(i*2+1+0.7712276549255324)*hx;
		x[32]= xa+(i*2+1-0.8072524984168955)*hx;
		x[33]= xa+(i*2+1+0.8072524984168955)*hx;
		x[34]= xa+(i*2+1-0.8404945765458014)*hx;
		x[35]= xa+(i*2+1+0.8404945765458014)*hx;
		x[36]= xa+(i*2+1-0.8708392975582413)*hx;
		x[37]= xa+(i*2+1+0.8708392975582413)*hx;
		x[38]= xa+(i*2+1-0.8981820578754266)*hx;
		x[39]= xa+(i*2+1+0.8981820578754266)*hx;
		x[40]= xa+(i*2+1-0.9224286030428122)*hx;
		x[41]= xa+(i*2+1+0.9224286030428122)*hx;
		x[42]= xa+(i*2+1-0.9434953534644419)*hx;
		x[43]= xa+(i*2+1+0.9434953534644419)*hx;
		x[44]= xa+(i*2+1-0.9613096946231363)*hx;
		x[45]= xa+(i*2+1+0.9613096946231363)*hx;
		x[46]= xa+(i*2+1-0.9758102337149845)*hx;
		x[47]= xa+(i*2+1+0.9758102337149845)*hx;
		x[48]= xa+(i*2+1-0.9869470350233716)*hx;
		x[49]= xa+(i*2+1+0.9869470350233716)*hx;
		x[50]= xa+(i*2+1-0.9946819193080071)*hx;
		x[51]= xa+(i*2+1+0.9946819193080071)*hx;
		x[52]= xa+(i*2+1-0.9989899477763282)*hx;
		x[53]= xa+(i*2+1+0.9989899477763282)*hx;
		           
		for(int j=0; j<N; j++)
		{
			double y[54];
			y[1 ]= ya+(j*2+1+0.0000000000000000)*hy;
			y[2 ]= ya+(j*2+1-0.0586850543002595)*hy;
			y[3 ]= ya+(j*2+1+0.0586850543002595)*hy;
			y[4 ]= ya+(j*2+1-0.1171678090719551)*hy;
			y[5 ]= ya+(j*2+1+0.1171678090719551)*hy;
			y[6 ]= ya+(j*2+1-0.1752466621553257)*hy;
			y[7 ]= ya+(j*2+1+0.1752466621553257)*hy;
			y[8 ]= ya+(j*2+1-0.2327214037242726)*hy;
			y[9 ]= ya+(j*2+1+0.2327214037242726)*hy;
			y[10]= ya+(j*2+1-0.2893939064516262)*hy;
			y[11]= ya+(j*2+1+0.2893939064516262)*hy;
			y[12]= ya+(j*2+1-0.3450688084957224)*hy;
			y[13]= ya+(j*2+1+0.3450688084957224)*hy;
			y[14]= ya+(j*2+1-0.3995541869539530)*hy;
			y[15]= ya+(j*2+1+0.3995541869539530)*hy;
			y[16]= ya+(j*2+1-0.4526622194618458)*hy;
			y[17]= ya+(j*2+1+0.4526622194618458)*hy;
			y[18]= ya+(j*2+1-0.5042098316571334)*hy;
			y[19]= ya+(j*2+1+0.5042098316571334)*hy;
			y[20]= ya+(j*2+1-0.5540193282770679)*hy;
			y[21]= ya+(j*2+1+0.5540193282770679)*hy;
			y[22]= ya+(j*2+1-0.6019190057137693)*hy;
			y[23]= ya+(j*2+1+0.6019190057137693)*hy;
			y[24]= ya+(j*2+1-0.6477437439165100)*hy;
			y[25]= ya+(j*2+1+0.6477437439165100)*hy;
			y[26]= ya+(j*2+1-0.6913355756013667)*hy;
			y[27]= ya+(j*2+1+0.6913355756013667)*hy;
			y[28]= ya+(j*2+1-0.7325442308075103)*hy;
			y[29]= ya+(j*2+1+0.7325442308075103)*hy;
			y[30]= ya+(j*2+1-0.7712276549255324)*hy;
			y[31]= ya+(j*2+1+0.7712276549255324)*hy;
			y[32]= ya+(j*2+1-0.8072524984168955)*hy;
			y[33]= ya+(j*2+1+0.8072524984168955)*hy;
			y[34]= ya+(j*2+1-0.8404945765458014)*hy;
			y[35]= ya+(j*2+1+0.8404945765458014)*hy;
			y[36]= ya+(j*2+1-0.8708392975582413)*hy;
			y[37]= ya+(j*2+1+0.8708392975582413)*hy;
			y[38]= ya+(j*2+1-0.8981820578754266)*hy;
			y[39]= ya+(j*2+1+0.8981820578754266)*hy;
			y[40]= ya+(j*2+1-0.9224286030428122)*hy;
			y[41]= ya+(j*2+1+0.9224286030428122)*hy;
			y[42]= ya+(j*2+1-0.9434953534644419)*hy;
			y[43]= ya+(j*2+1+0.9434953534644419)*hy;
			y[44]= ya+(j*2+1-0.9613096946231363)*hy;
			y[45]= ya+(j*2+1+0.9613096946231363)*hy;
			y[46]= ya+(j*2+1-0.9758102337149845)*hy;
			y[47]= ya+(j*2+1+0.9758102337149845)*hy;
			y[48]= ya+(j*2+1-0.9869470350233716)*hy;
			y[49]= ya+(j*2+1+0.9869470350233716)*hy;
			y[50]= ya+(j*2+1-0.9946819193080071)*hy;
			y[51]= ya+(j*2+1+0.9946819193080071)*hy;
			y[52]= ya+(j*2+1-0.9989899477763282)*hy;
			y[53]= ya+(j*2+1+0.9989899477763282)*hy;
	
			for(int p=1; p<54; p++)
				for(int q=1; q<54; q++)
					total+=weight[p]*weight[q]*fun(qq,x[p],y[q]);
		}
	}
	total*=hx*hy;
	return total;
}

double percent=0;
double ddd=0;

double adaptive2d_polar( double (*fun)(double,double,double)
		,double qq,double xa, double xb,double ya,double yb)
{
	double torr = 1e-17;
	int th=32;
	double total=0;
	double hx = (xb-xa)/th;
	double hy = (yb-ya)/th;
	percent=0;
	ddd=0;
	#pragma omp parallel for reduction (+:total)
	for(int index=0; index<th*th; index++)
	{
		int i = index%th;
		int j = index/th;
		double initx = xa + hx*i;
		double endx = xa+hx*(i+1);
		double inity = ya + hy*j;
		double endy = ya+hy*(j+1);
		total+=adaptive2d_polar_internal(fun,qq,initx,endx,inity,endy,0,torr);
	}
	return total;
}
double adaptive2d_polar_internal(double (*fun)(double,double,double)
		,double qq,double xa,double xb,double ya,double yb,int step,double torr)
{
	double res1 = gq2d_polar(fun,qq,xa,xb,ya,yb,4);
	double res2 = gq2d_polar(fun,qq,xa,xb,ya,yb,1);
	double midx = 0.5*(xa+xb);
	double midy = 0.5*(ya+yb);
//	printf("%f~%f\t\t%f~%f\t%d\n",xa,xb,ya,yb,step);
	if( fabs(res1-res2)>torr  && step<10 )
		return adaptive2d_polar_internal(fun,qq,xa,midx,ya,midy,step+1,torr)
			+adaptive2d_polar_internal(fun,qq,midx,xb,ya,midy,step+1,torr)
			+adaptive2d_polar_internal(fun,qq,xa,midx,midy,yb,step+1,torr)
			+adaptive2d_polar_internal(fun,qq,midx,xb,midy,yb,step+1,torr);
	else
	{
		double dstep = 1.0/32.0/32.0/pow(4,step);
		#pragma omp critical
		{
			percent+= dstep;
			if(percent>ddd)
			{
			//	printf("%f\n",percent);
				ddd+=0.1;
			}

		}
		//	printf("%f\n",percent);
		return res1;
	}
}
