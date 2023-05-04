////equilibrium.cpp---------check if the system is in equilibrium 
#include <stdio.h>
#include <math.h>
#include <stdlib.h> //rand() and functions of this kind, mainly used here
#include <time.h>

////C=20mM,sigma=2.2A,using dimensionless length
#define pi (4*atan(1)) 
#define e 2.71828182845904523536
#define kB 1.380649e-23
#define R   (pow(2.0,1.0/6.0))
#define eps (5.0/6.0) //in E Etemp corresponding to R
#define lb 1.7

////system parameters
#define sigma 1  //only F uses
#define lx 51.2
#define ly 51.2 
#define lz 51.2 //lx,ly,lz side length of box£¬10 Debye length
#define L0 1.1 //distance between 2 monomers in one chain 
#define N 25 //41
#define n1 170 //202//+1
#define n2 120 //120//-1
#define n3 30 //+3
#define n4 30 //-3
#define r0 0.5 //radius of monomers
#define ri 0.5 //radius of ions 
#define nc 31 //crowding number
#define rc 7 //radius of crowdings
#define nion (n1+n2+n3+n4)
#define nall (n1+n2+n3+n4+nc)

////modeling parameters
#define step0 10
#define stepnc 1
//#define allstep ((n1+n2)*100.0)
//#define bin 50
//#define D 8
//#define sam 100.0  //sample step


////global variable
double temp=0.0;
//double U[(long long)allstep]={0.0};
int countaccept=0;  //count accept
int countnc=0;

////function prototype
double distance(double x1,double y1,double z1,double x2,double y2,double z2);
//void judgement(int i,double *x1,double *y1,double *z1,double *x2,double *y2,double *z2); 
//void movei(double Ei,int i,const double *x1,const double *y1,const double *z1,const double *q1,double *x2,double *y2,double *z2,const double *q2);
int move(int i,double *x1,double *y1, double *z1,double *q1,double *x2,double *y2,double *z2,double *q2);
double E(int i,double *x1,double *y1, double *z1,double *q1,double *x2,double *y2,double *z2,double *q2);
double Etemp(int i,double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2,double xtemp,double ytemp,double ztemp);
double pmf(double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2);
void print_location(FILE *fx,FILE *fy,FILE *fz,double *x2,double *y2,double *z2);
void print_number(FILE *fp,double *x2,double *y2,double D);
//void print_curve(FILE *fp,const double *x2,const double *y2);
double power6(double x);
double power7(double x);
double nearest(double ds, double init, double ls);
double corr(double ds, double init, double ls);



////main function
int main(int argc, char *argv[])
{
	double sam;
	double allstep;
	double D;
	D= strtod (argv[1],NULL);  //argv[1] is the distance between 2 rods, also the number behind the filename
	
	FILE *fp, *fq, *fx ,*fy, *fz, *fxnew, *fynew, *fznew;
	fp=fopen("eqstep.txt","rb");  //read in the equilibrium time
	int Dd;
        Dd= (int)(10*D);
	
        char npmf[20];
        char xlocal[20];
        char ylocal[20];
        char zlocal[20];
//-----------------------------
	char xnew[20];
	char ynew[20];
	char znew[20];
//-----------------------------
        sprintf(npmf,"npmf%d.dat",Dd);
        sprintf(xlocal,"x%d.dat",Dd);
        sprintf(ylocal,"y%d.dat",Dd);
        sprintf(zlocal,"z%d.dat",Dd);
//-------------------------------------
	sprintf(xnew,"xnew%d.dat",Dd);
	sprintf(ynew,"ynew%d.dat",Dd);
	sprintf(znew,"znew%d.dat",Dd);
//-------------------------------------
	fq=fopen(npmf,"wb+"); //record force
	fx=fopen(xlocal,"rb");//open the last location
	fy=fopen(ylocal,"rb");
	fz=fopen(zlocal,"rb");
//------------------------------
	fxnew=fopen(xnew,"wb+");
	fynew=fopen(ynew,"wb+");
	fznew=fopen(znew,"wb+");
	//fmean=fopen("fmean.dat","a+");
	
	////if we cannot open the file
	if (fx==NULL||fy==NULL||fz==NULL||fp==NULL){
		printf("File open error!");
		return -1;
	}
	

	////pick out the corresponding tau
	int list=0;
 	double rank[70],tau[70];
	while(!feof(fp))
        {
                fscanf(fp,"%lf %lf",&rank[list],&tau[list]);
                list++;
        }
	for(list=0;list<40;list++)
        {
                if ((int)(10*D)==(int)(rank[list])){sam=2*tau[list];break;}
        }	
	
	allstep=8000*sam;
	//allstep=1.0;
	

	////1 are RNA monomers, 2 are ions, 3 are crowdings
	double x1[N+N]={0.0};
	double x2[nall]={0.0};
	//double x3[nc]={0.0}; 
	double y1[N+N]={0.0};
	double y2[nall]={0.0};
	//double y3[nc]={0.0};
	double z1[N+N]={0.0};
	double z2[nall]={0.0};
	//double z3[nc]={0.0};
	double q1[N+N]={0.0};
	double q2[nall]={0.0};
	//double q3[nc]={0.0};  
	
	////Assign values to monomers£¬RNA kept still 
	int i;
	for (i=0;i<N;i++){
		q1[i]=-1.0;
		x1[i]=-D/2.0;
		y1[i]=0.0;
		z1[i]=((i+i)-N+1)*L0/2.0;//correct, monomer at the bottom is 0, then central monomer is (N+1)/2-1=(N-1)/2
		q1[i+N]=-1.0;
		x1[i+N]=D/2.0;
		y1[i+N]=0.0;
		z1[i+N]=((i+i)-N+1)*L0/2.0;
	}

	for (i=0;i<n1;i++){q2[i]=1.0;}
	for (i=n1;i<(n1+n2);i++){q2[i]=-1.0;}
	for (i=(n1+n2);i<(n1+n2+n3);i++){q2[i]=3.0;}
	for (i=(n1+n2+n3);i<(n1+n2+n3+n4);i++){q2[i]=-3.0;}
//	for (i=(n1+n2+n3+n4);i<nion;i++){q2[i]=0.0;}


//for (i=0;i<nion;i++){printf("%d\t%f\n",i,q2[i]);}

	
	////put ions randomly, i.e.ions' three-dimensional values are random
	srand((unsigned)time(NULL));
	
	double initialE;	
	////place ions

	for (i=0;i<nall;i++){
		//x[i]=i;
		fscanf(fx,"%lf\n",&x2[i]);
		fscanf(fy,"%lf\n",&y2[i]);
		fscanf(fz,"%lf\n",&z2[i]);
		//printf("%d\t%lf\n",i,x2[i]);
	}
	
	printf("Initial finished!\n");	
	
//**************************************************************************************************
	
	////move many steps to reach equilibrium
	srand((unsigned)time(NULL));
	double step1;
	double aveF=0.0;
	for (step1=1.0;step1<allstep+1;step1++)
	{
		
		////move ions
		i=(int) (rand()%(nall));		
		move(i,x1,y1,z1,q1,x2,y2,z2,q2);  //function call		

		////sample and
		//aveF+=pmf(x1,y1,z1,q1,x2,y2,z2,q2);
		if (0==(long long)step1%(long)sam) { 
			double aveF=0.0;
			aveF=pmf(x1,y1,z1,q1,x2,y2,z2,q2);
			fprintf(fq,"%f\n",aveF);
			fflush(stdout);
			//aveF=0;
			printf("%f\n",step1);
		}   
	
	
		////print number in a certain area
		//print_number(fq,x2,y2,D);  //function call
	
	//	if (0==fmod(step1,1000.0))	{printf("%f\n",step1);}
					
	}
	
	double meanF;
	//printf("allf=%f\n",aveF);
	//printf("divide=%d\n",(int)(allstep/sam));
	//meanF=aveF/((int)(allstep/sam));
	
	//printf("bin=%d\n",(int)(allstep/sam));
	//printf("mean_force=%f\n",meanF);
	
	//fprintf(fmean,"%f\t %f\n",D,meanF);
	//fflush(stdout);
	
	
	printf("Sample finished!\n");
	printf("couacc=%d\n",countaccept);
	
	print_location(fxnew,fynew,fznew,x2,y2,z2);

	fclose(fp);	
	fclose(fq);
	fclose(fx);
	fclose(fy);
	fclose(fz);
	fclose(fxnew);
	fclose(fynew);
	fclose(fznew);
	//fclose(fmean);

	double acco=countaccept; 
	printf("countaccept=%f\n",acco);
	printf("allstep=%f\n",allstep);
	printf("acceptance rate=%f\n",acco/(double) allstep);
	printf("acceptance nc=%f\n",countnc/(double) allstep);

		 
	return 0;
	
	 
}




////--------------------------------------------------------------------------------------------------main function




////distance function, calculating center-to-center distance of two particles
double distance(double x1,double y1,double z1,double x2,double y2,double z2)   //function definition
{
	double d2;
	double d;
	d2=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
	//printf("%f\t",d2);
	d=sqrt(d2);
	return d; 
} 

//*****************************************************************************

 
 
////random-move function---every ion
int move(int i,double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2)  //function definition
{
	////move a certain step length, (rand()%100-49.5)/50.0 range from -1~1
	
	double Ei;
	Ei=E(i,x1,y1,z1,q1,x2,y2,z2,q2);

	//step0 is the step of ions 
	double tx=((rand()*1.0/RAND_MAX)-0.5)*step0;     
	double ty=((rand()*1.0/RAND_MAX)-0.5)*step0; 
	double tz=((rand()*1.0/RAND_MAX)-0.5)*step0; 
	
	//step1 is the step of crowders
	double Tx=((rand()*1.0/RAND_MAX)-0.5)*stepnc;     
	double Ty=((rand()*1.0/RAND_MAX)-0.5)*stepnc; 
	double Tz=((rand()*1.0/RAND_MAX)-0.5)*stepnc; 
	
	//temporary variable
	double xtemp;
	double ytemp;
	double ztemp;
	
	//printf("tx=%f\n",tx);
	if (i<nion){
	xtemp=(double) x2[i]+tx;
	ytemp=(double) y2[i]+ty;
	ztemp=(double) z2[i]+tz;   //pass test1, correct
	}
	else {
	xtemp=(double) x2[i]+Tx;
	ytemp=(double) y2[i]+Ty;
	ztemp=(double) z2[i]+Tz;
	}
	
	
	////periodic boundary condition
	if(xtemp>lx/2.0) {xtemp=xtemp-lx;}
	if(ytemp>ly/2.0) {ytemp=ytemp-ly;}
	if(ztemp>lz/2.0) {ztemp=ztemp-lz;}
	if(xtemp<-lx/2.0) {xtemp=xtemp+lx;}
	if(ytemp<-lx/2.0) {ytemp=ytemp+ly;}
	if(ztemp<-lx/2.0) {ztemp=ztemp+lz;} 
	
	
	//printf("x=%f\ty=%f\tz=%f\n",x2[i],y2[i],z2[i]);
	////the same with judgement
	////ion and monomer
	int j;
	double d_im;
	for (j=0;j<N+N;j++){
		d_im=distance(x1[j],y1[j],z1[j],xtemp,ytemp,ztemp);
		if (i>(n1+n2+n3+n4)){ if (d_im<(rc+r0)) return 11;}
		else {if (d_im<(ri+r0)) return 12;}
	}
	////between ions
	int k;
	double d_ii;
	for (k=0;k<nion;k++){
		if(k-i){
		d_ii=distance(x2[k],y2[k],z2[k],xtemp,ytemp,ztemp);
		if (k<(n1+n2+n3+n4)&&i<(n1+n2+n3+n4)){if(d_ii<(ri+ri)) return 21;}
		else if(k>=(n1+n2+n3+n4)&&i>=(n1+n2+n3+n4)){if(d_ii<(rc+rc)) return 22;}
		else {if(d_ii<(ri+rc)) return 23;}
		}
	}
	////metroplis
	////first calculate the final evergy Ef
	double Ef;
	Ef=Etemp(i,x1,y1,z1,q1,x2,y2,z2,q2,xtemp,ytemp,ztemp);
	double dE;
	dE=Ef-Ei;
	//printf("dE=%f\n",dE);
	//printf("Ei=%f\tEf=%f\n",Ei,Ef);
	if (dE>0){	
		double r;
		r=rand()*1.0/RAND_MAX;
//printf("r=%f\n",r);
		double b;
		b=exp(-dE);
		//printf("b=%f\n",b);
		if (r>b){
			//printf("x3=%f\n",x2[i]);
			return 3;   
		}		
	}
	
	
	////pass the Ef-move to Ei-main, pass the temp position to original array
	x2[i]=xtemp;
	y2[i]=ytemp;
	z2[i]=ztemp;
	
	if (i<nion){++countaccept;}
	else {++countnc;}	
	//fprintf(fp,"%.12f\n",Ef);
	return 4;
} 
 
 
//****************************************************************************


////Coulomb: nearest
double nearst(double ds, double init, double ls)
{	
	double sion=0.0;
	if (ds>ls/2.0){sion=init-ls;}
        else if(ds<-ls/2.0){sion=init+ls;}
	else {sion=init;}
	return sion;
}

////Coulomb: correlation (first correlation, second nearest)
double corr(double ds, double init, double ls)
{
	double correlation=0.0;
	if ((0.0<=ds)&&(ds<ls/2.0)){correlation=init-ls;}
	else if ((-ls/2.0<ds)&&(ds<0.0)){correlation=init+ls;}
	else {correlation=init;}
	return correlation;
}



////system-energy function
double E(int i,double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2)  //function definition
{
	double u=0.0;
	////the ion and monomers
	int j;
	//double d_im;
	//double d_im2;
	for (j=0;j<N+N;j++){
		double xion=0.0;
                double yion=0.0;
                double zion=0.0;
		double xcorr=0.0;
		double ycorr=0.0;
		double zcorr=0.0;
		double dx=x1[j]-x2[i];
		double dy=y1[j]-y2[i];
		double dz=z1[j]-z2[i];
                //if (dx>lx/2.0){xion=x1[j]-lx;}
		//else if ((0.0<dx)&&(dx<lx/2.0){xion=x1[j];}
		//else if ((-lx/2.0<dx)&&(dx<0.0)){xion=x1[j];}
                //else if (x1[j]-x2[i]<-lx/2.0){xion=x1[j]+lx;}
                //else {xion=x1[j];}
                //if (y1[j]-y2[i]>ly/2.0){yion=y1[j]-ly;}
                //else if (y1[j]-y2[i]<-ly/2.0){yion=y1[j]+ly;}
                //else {yion=y1[j];}
                //if (z1[j]-z2[i]>lz/2.0){zion=z1[j]-lz;}
                //else if (z1[j]-z2[i]<-lz/2.0){zion=z1[j]+lz;}
                //else {zion=z1[j];}
		xion=nearst(dx,x1[j],lx);
		yion=nearst(dy,y1[j],ly);
		zion=nearst(dz,z1[j],lz);
		//xcorr=corr(dx,x1[j],lx);
		//ycorr=corr(dy,y1[j],ly);
		//zcorr=corr(dz,z1[j],lz);
		double d_im;
		//double d_im2;
		d_im=distance(x1[j],y1[j],z1[j],x2[i],y2[i],z2[i]);
		//d_im2=distance(xcorr,ycorr,zcorr,x2[i],y2[i],z2[i]);
		//printf("dim=%f\n",d_im);
		//printf("q1=%f\n",q1[j]);
		if (d_im<R){
			u+=lb*q2[i]*q1[j]/d_im+10.0/3.0*power6(1.0/d_im)*power6(1.0/d_im)-10.0/3.0*power6(1.0/d_im)+eps;
			//printf("u1=%f\t",u);
		}
		else {
			u+=lb*q2[i]*q1[j]/d_im;
			//printf("u2=%f\t",u);
		}
	}
	////the ion and other ions
	int k;
	//double d_ii;
	for (k=0;k<nion;k++)
	{
		double xion=0.0;
                double yion=0.0;
                double zion=0.0;
		double xcorr=0.0;
                double ycorr=0.0;
                double zcorr=0.0;
                double dx=x2[k]-x2[i];
                double dy=y2[k]-y2[i];
                double dz=z2[k]-z2[i];
                //if (x2[k]-x2[i]>lx/2.0){xion=x2[k]-lx;}
                //else if (x2[k]-x2[i]<-lx/2.0){xion=x2[k]+lx;}
                //else {xion=x2[k];}
                //if (y2[k]-y2[i]>ly/2.0){yion=y2[k]-ly;}
                //else if (y2[k]-y2[i]<-ly/2.0){yion=y2[k]+ly;}
                //else {yion=y2[k];}
                //if (z2[k]-z2[i]>lz/2.0){zion=z2[k]-lz;}
                //else if (z2[k]-z2[i]<-lz/2.0){zion=z2[k]+lz;}
                //else {zion=z2[k];}
		xion=nearst(dx,x2[k],lx);
                yion=nearst(dy,y2[k],ly);
                zion=nearst(dz,z2[k],lz);
                //xcorr=corr(dx,x2[k],lx);
                //ycorr=corr(dy,y2[k],ly);
                //zcorr=corr(dz,z2[k],lz);
		double d_ii;
		double d_ii2;
		d_ii=distance(x2[k],y2[k],z2[k],x2[i],y2[i],z2[i]);
		//d_ii2=distance(xcorr,ycorr,zcorr,x2[i],y2[i],z2[i]);
		//printf("dii=%f\n",d_ii);
		if (i-k){
			if (d_ii<R){
				u+=lb*q2[i]*q2[k]/d_ii+10.0/3.0*power6(1.0/d_ii)*power6(1.0/d_ii)-10.0/3.0*power6(1.0/d_ii)+eps;
				//printf("u3=%f\t",u);
			}
			else {
				u+=lb*q2[i]*q2[k]/d_ii;
				//printf("u4=%f\t",u);
			}
		}
	}
	//printf("u=%f\n",u);
	return u;
}


//*****************************************************************************

////changed-energy function
double Etemp(int i,double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2,double xtemp,double ytemp,double ztemp)  //function definition
{
	double u=0.0;
	////the ion and monomers
	int j;
	//double d_im;
	for (j=0;j<N+N;j++){
		double xion=0.0;
		double yion=0.0;
		double zion=0.0;
		double xcorr=0.0;
                double ycorr=0.0;
                double zcorr=0.0;
                double dx=x1[j]-xtemp;
                double dy=y1[j]-ytemp;
                double dz=z2[j]-ztemp;
		//if (x1[j]-xtemp>lx/2.0){xion=x1[j]-lx;}
                //else if (x1[j]-xtemp<-lx/2.0){xion=x1[j]+lx;}
                //else {xion=x1[j];}
                //if (y1[j]-ytemp>ly/2.0){yion=y1[j]-ly;}
                //else if (y1[j]-ytemp<-ly/2.0){yion=y1[j]+ly;}
                //else {yion=y1[j];}
                //if (z1[j]-ztemp>lz/2.0){zion=z1[j]-lz;}
                //else if (z1[j]-ztemp<-lz/2.0){zion=z1[j]+lz;}
                //else {zion=z1[j];}		
		xion=nearst(dx,x1[j],lx);
                yion=nearst(dy,y1[j],ly);
                zion=nearst(dz,z1[j],lz);
                //xcorr=corr(dx,x1[j],lx);
                //ycorr=corr(dy,y1[j],ly);
                //zcorr=corr(dz,z1[j],lz);
		double d_im;
		double d_im2;
		d_im=distance(x1[j],y1[j],z1[j],xtemp,ytemp,ztemp);
		//d_im2=distance(xcorr,ycorr,zcorr,xtemp,ytemp,ztemp);
		//printf("dim=%f\n",d_im);
		//printf("q1=%f\n",q1[j]);
		if (d_im<R){
			u+=lb*q2[i]*q1[j]/d_im+10.0/3.0*power6(1.0/d_im)*power6(1.0/d_im)-10.0/3.0*power6(1.0/d_im)+eps;
			//printf("u1=%f\t",u);
		}
		else {
			u+=lb*q2[i]*q1[j]/d_im;
			//printf("u2=%f\t",u);
		}
	}
	////the ion and other ions
	int k;
	//double d_ii;
	for (k=0;k<nion;k++)
	{
		double xion=0.0;
		double yion=0.0;
		double zion=0.0;
		double xcorr=0.0;
                double ycorr=0.0;
                double zcorr=0.0;
                double dx=x2[k]-xtemp;
                double dy=y2[k]-ytemp;
                double dz=z2[k]-ztemp;
		//if (x2[k]-xtemp>lx/2.0){xion=x2[k]-lx;}
		//else if (x2[k]-xtemp<-lx/2.0){xion=x2[k]+lx;}
		//else {xion=x2[k];}
		//if (y2[k]-ytemp>ly/2.0){yion=y2[k]-ly;}
		//else if (y2[k]-ytemp<-ly/2.0){yion=y2[k]+ly;}
		//else {yion=y2[k];}
		//if (z2[k]-ztemp>lz/2.0){zion=z2[k]-lz;}
		//else if (z2[k]-ztemp<-lz/2.0){zion=z2[k]+lz;}
		//else {zion=z2[k];}
		xion=nearst(dx,x2[k],lx);
                yion=nearst(dy,y2[k],ly);
                zion=nearst(dz,z2[k],lz);
                //xcorr=corr(dx,x2[k],lx);
                //ycorr=corr(dy,y2[k],ly);
                //zcorr=corr(dz,z2[k],lz);
		double d_ii;
		//double d_ii2;
		d_ii=distance(x2[k],y2[k],z2[k],xtemp,ytemp,ztemp);
		//d_ii2=distance(xcorr,ycorr,zcorr,xtemp,ytemp,ztemp);
		//printf("dii=%f\n",d_ii);
		if (i-k){
			if (d_ii<R){
				u+=lb*q2[i]*q2[k]/d_ii+10.0/3.0*power6(1.0/d_ii)*power6(1.0/d_ii)-10.0/3.0*power6(1.0/d_ii)+eps;
				//printf("u3=%f\t",u);
			}
			else {
				u+=lb*q2[i]*q2[k]/d_ii;
				//printf("u4=%f\t",u);
			}
		}
	}
	//printf("utemp=%f\n",u);
	return u;
}



//***************************************************************************

double pmf(double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2)
{
	int i;
	int j;
	int k;
	//double d_mm;
	//double d_mm2;
	//double d_im;
	//double d_im2;
	double Fle=0;
	double Fre=0;
////away is positive, closer is negative
	for (i=0;i<N;i++){
		//printf("Fle=%f\n",Fle);
		for (j=N;j<N+N;j++){
			/*double xion=0.0;
			double yion=0.0;
			double zion=0.0;
			double xcorr=0.0;
			double ycorr=0.0;
			double zcorr=0.0;
			double dx=x1[j]-x1[i];
			double dy=y1[j]-y1[i];
			double dz=z1[j]-z1[i];
			if (z1[j]-z1[i]<-lz/2.0){zion=z1[j]+lz;}
			else if (z1[j]-z1[i]>lz/2.0){zion=z1[j]-lz;}
			else {zion=z1[j];}
			xion=nearst(dx,x1[j],lx);
			yion=nearst(dy,y1[j],ly);
			zion=nearst(dz,z1[j],lz);
			xcorr=corr(dx,x1[j],lx);
			ycorr=corr(dy,y1[j],ly);
			zcorr=corr(dz,z1[j],lz);*/
			double d_mm;
			//double d_mm2;
			d_mm=distance(x1[i],y1[i],z1[i],x1[j],y1[j],z1[j]);
			//d_mm2=distance(x1[i],y1[i],z1[i],xcorr,ycorr,zcorr);
			if(d_mm<R){
				Fle+=(x1[j]-x1[i])*lb/(d_mm*d_mm*d_mm*sigma)+(x1[j]-x1[i])*(40.0*power7(1.0/d_mm)*power6(1.0/d_mm)-20.0*power7(1.0/d_mm))/(d_mm*sigma);
				}
		//	else if((x1[j]-x1[i])>lx/2.0){
		//		Fle+=(x1[j]-lx-x1[i])*3.2/(d_mm*d_mm*d_mm*sigma);}
		//	else if((x1[j]-x1[i])<-lx/2.0){
		//		Fle+=(x1[j]+lx-x1[i])*3.2/(d_mm*d_mm*d_mm*sigma);}
			else{	
				Fle+=(x1[j]-x1[i])*lb/(d_mm*d_mm*d_mm*sigma);    //neglect q1[i] and q1[j] since they multiply as 1
			}
		}
		
		for (k=0;k<nion;k++){
			double deltax;
			double deltax2;
			double xion=0.0;
			double yion=0.0;
			double zion=0.0;
			double xcorr=0.0;
			double ycorr=0.0;
			double zcorr=0.0;
			double dx=x2[k]-x1[i];
			double dy=y2[k]-y1[i];
			double dz=z2[k]-z1[i];
			//if (x2[k]>x1[i]+lx/2.0){xion=x2[k]-lx;}
			//else {xion=x2[k];}
			//if (z2[k]-z1[i]<-lz/2.0){zion=z2[k]+lz;}
			//else if (z2[k]-z1[i]>lz/2.0){zion=z2[k]-lz;}
			//else {zion=z2[k];}
			xion=nearst(dx,x2[k],lx);
			yion=nearst(dy,y2[k],ly);
			zion=nearst(dz,z2[k],lz);
                	//xcorr=corr(dx,x2[k],lx);
			//ycorr=corr(dy,y2[k],ly);
                      	//zcorr=corr(dz,z2[k],lz);						
			double d_im;
			double d_im2;
			d_im=distance(x1[i],y1[i],z1[i],xion,y2[k],zion);
			//d_im2=distance(x1[i],y1[i],z1[i],xcorr,ycorr,zcorr);
			deltax=xion-x1[i];
			//deltax2=xcorr-x1[i];		
			//dx=x2[k]-x1[i];
			//if(dx>lx/2.0) dx=-lx+dx;
			//if(dx<-lx/2.0) dx=lx+dx;
			if(d_im<R){
				Fle+=lb*q1[i]*q2[k]/(d_im*d_im*d_im*sigma)*deltax+(40.0*power7(1.0/d_im)*power6(1.0/d_im)-20.0*power7(1.0/d_im))*deltax/(d_im*sigma);
				}
				else{
					Fle+=lb*q1[i]*q2[k]/(d_im*d_im*d_im*sigma)*deltax;
				}//printf("%f\n",Fle);
		}
	}
	/*
	for (i=N;i<N+N;i++){
		printf("Fre=%f\n",Fre);
		for (j=0;j<N;j++){
			d_mm=distance(x1[i],y1[i],z1[i],x1[j],y1[j],z1[j]);
			Fre+=3.2/(d_mm*d_mm*d_mm*sigma)*(x1[i]-x1[j]);    //neglect q1[i] and q1[j] since they multiply as 1
		}
		
		for (k=0;k<n1+n2;k++){
			double dx;
			d_im=distance(x1[i],y1[i],z1[i],x2[k],y2[k],z2[k]);
			dx=x2[i]-x1[k];
			if(dx>lx/2.0) dx=-lx+dx;
			if(dx<-lx/2.0) dx=lx+dx;
			if(d_im<R){
				Fre+=3.2*q1[i]*q2[k]/(d_im*d_im*d_im*sigma)*dx+(40.0*power7(1/(sigma*d_im))*power6(1/(sigma*d_im))-20.0*power7(1/(sigma*d_im)))*dx/d_im;
				}
				else{
					Fre+=3.2*q1[i]*q2[k]/(d_im*d_im*d_im*sigma)*dx;
				}
		}
	}
	*/	
	double F;
	F=Fle;
	//printf("F=%f\n",F);
	
	return F;
}


//***************************************************************************************

////print function
void print_number(FILE *fp,double *x2,double *y2,double D)   //function definition
{////print the number of ions in a certain area
	int i;
	int count=0;
	for (i=0;i<nion;i++)
	{
		//printf("%d\t",i);
		double area=(5.0*r0)*(5.0*r0);
		if (((x2[i]-D/2)*(x2[i]-D/2)+y2[i]*y2[i])<area){
			count++;
		}
	}
	fprintf(fp,"%d\n",count);
	fflush(stdout);
}


////print function
void print_location(FILE *fx,FILE *fy,FILE *fz,double *x2,double *y2,double *z2)   //function definition
{////print the number of ions in a certain area
	int i;
	for (i=0;i<nion;i++)
	{
		fprintf(fx,"%f\n",x2[i]);
		fprintf(fy,"%f\n",y2[i]);
		fprintf(fz,"%f\n",z2[i]);
		fflush(stdout);
		}
}





double power6(double x)
{
	double y;
	y=x*x*x*x*x*x;
	return y;
}

double power7(double x)
{
	double y;
	y=x*x*x*x*x*x*x;
	return y;
}


////---------------------------------------------------------------------------------------my_functions




/* 


////judgement-assign function, if overlapped-recursion, if not-assign value
void judgement(int i,double *x1,double *y1,double *z1,double *x2,double *y2,double *z2) //function definition
{
	////for each ion i, test the overlapping of ith-ion and j-monomers & ith-ion and i-1-ions
	////assignment is inside for recursion
	//srand((unsigned)time(NULL));
		x2[i]=(double) (rand()%1001)/1000.0*lx-lx/2.0;  //original (rand()%100)/100*lx-lx/2.0
		y2[i]=(double) (rand()%1001)/1000.0*ly-ly/2.0;
		z2[i]=(double) (rand()%1001)/1000.0*lz-lz/2.0;
		////ion and monomer
		int j;
		double d_im;
		for (j=0;j<N+N;j++){
			d_im=distance(x1[j],y1[j],z1[j],x2[i],y2[i],z2[i]);
			if (d_im<(ri+r0)){
				judgement(i,x1,y1,z1,x2,y2,z2);  //function call
			}
		}
		////between ions
		int k;
		double d_ii;
		for (k=0;k<i;k++){
			d_ii=distance(x2[k],y2[k],z2[k],x2[i],y2[i],z2[i]);
			if(d_ii<(ri+ri)){
				judgement(i,x1,y1,z1,x2,y2,z2);  //function call
			}
		}
		////update ith ion xyz
}





////random-move function---every ion
void movei(double Ei,int i,const double *x1,const double *y1,const double *z1,const double *q1,double *x2,double *y2,double *z2,const double *q2)  //function definition
{
	////move a certain step length, (rand()%100-49.5)/50.0 range from -1~1
	static long ountall;
	static long countaccep;
	countall++;
	double tx=((rand()%101)-50)/50.0*step0;     //original (rand()%100-49.5)/50.0*step0;
	double ty=((rand()%101)-50)/50.0*step0;
	double tz=((rand()%101)-50)/50.0*step0;
	x2[i]=(double) x2[i]+tx;
	y2[i]=(double) y2[i]+ty;
	z2[i]=(double) z2[i]+tz;   //pass test1, correct
	////the same with judgement
	////ion and monomer
	int j;
	double d_im;
	for (j=0;j<N+N;j++){
		d_im=distance(x1[j],y1[j],z1[j],x2[i],y2[i],z2[i]);
		if (d_im<(ri+r0)){
			movei(Ei,i,x1,y1,z1,q1,x2,y2,z2,q2);  //function call
		}
	}
	////between ions
	int k;
	double d_ii;
	for (k=0;k<i;k++){
		d_ii=distance(x2[k],y2[k],z2[k],x2[i],y2[i],z2[i]);
		if(d_ii<(ri+ri)){
			movei(Ei,i,x1,y1,z1,q1,x2,y2,z2,q2);  //function call
		}
	}
	////more than judgement
	////metroplis
	////first calculate the final evergy Ef
	double Ef;
	Ef=E(i,x1,y1,z1,q1,x2,y2,z2,q2);
	double dE;
	dE=Ef-Ei;
	if (dE>0){
		double r=(rand()%1001)/1000;
		double b=exp(-dE);
		if (r>b){
			movei(Ei,i,x1,y1,z1,q1,x2,y2,z2,q2);   //function call
		}		
	}
	////periodic boundary condition
	if(x2[i]>lx/2.0) {x2[i]=x2[i]-lx;}
	if(y2[i]>ly/2.0) {y2[i]=y2[i]-ly;}
	if(z2[i]>lz/2.0) {z2[i]=z2[i]-lz;}
	if(x2[i]<-lx/2.0) {x2[i]=x2[i]+lx;}
	if(y2[i]<-lx/2.0) {y2[i]=y2[i]+ly;}
	if(z2[i]<-lx/2.0) {z2[i]=z2[i]+lz;} 
	
	////pass the Ef-move to Ei-main
	countaccept++;
	temp=Ef;
}

*/









/*
print_curve(FILE *fp,const double *x2,const double *y2)  //function definition
{////print density of several areas
	int  i;
	double density[k][100]={0.0};
	for (i=0;)
	
}
*/


