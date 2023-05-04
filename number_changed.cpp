////equilibrium.cpp---------check if the system is in equilibrium 
#include <stdio.h>
#include <math.h>
#include <stdlib.h> //rand() and functions of this kind, mainly used here
#include <time.h>

////r0=2.1A£¬C=20mM£¬sigma=2.2A using dimensionless length
#define pi (4*atan(1)) 
#define e 2.71828182845904523536
#define kB 1.380649e-23
#define R (pow(2.0,1.0/6.0))
#define eps (5.0/6.0)
#define lb 1.7

////system parameters
#define lx 51.2 
#define ly 51.2 
#define lz 51.2 //lx,ly,lz side length of box£¬10 Debye length
#define L0 1.1 //distance between 2 monomers in one chain 
#define N  25 //41
#define n1  170 //202//+1
#define n2  120 //120//-1
#define n3 30 //+3
#define n4 30 //-3
#define r0 0.5 //radius of monomers
#define ri 0.5 //radius of ions 
#define nc 31 //crowding number
#define rc 7 //radius of crowdings
#define nion (n1+n2+n3+n4) //number of all ions
#define nall (n1+n2+n3+n4+nc) //number of crowders


////modeling parameters
#define step0 10
#define stepnc 0.01
//#define allstep ((n1+n2)*500.0)
//#define bin 50
//#define D 4


////global variable
double temp=0.0;
//double U[(long long)allstep]={0.0};
long long countaccept=0;  //count accept of ion moving
long long  countnc=0;  //count accept of crowder moving

////function prototype
double distance(double x1,double y1,double z1,double x2,double y2,double z2);
void judgement(int i,double *x1,double *y1,double *z1,double *x2,double *y2,double *z2); 
//void movei(double Ei,int i,const double *x1,const double *y1,const double *z1,const double *q1,double *x2,double *y2,double *z2,const double *q2);
int move(int i,double *x1,double *y1, double *z1,double *q1,double *x2,double *y2,double *z2,double *q2);
//int cmove(int i,double *x1,double *y1, double *z1,double *q1,double *x2,double *y2,double *z2,double *q2);
double E(int i,double *x1,double *y1, double *z1,double *q1,double *x2,double *y2,double *z2,double *q2);
double Etemp(int i,double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2,double xtemp,double ytemp,double ztemp);
void print_location(FILE *fx,FILE *fy,FILE *fz,double *x2,double *y2,double *z2);
void print_number(FILE *fp,double *x2,double *y2,double D);
//void print_curve(FILE *fp,const double *x2,const double *y2);
double power6(double x);




////main function
int main(int argc, char *argv[])
{
//printf("%d\n",nion);
	double D;
	D= strtod (argv[1],NULL);  //argv[1] is the distance between 2 rods, also the number behind the filename
	double allstep;
	allstep=nion*strtod (argv[2],NULL);  //argv[2] is the runnning allstep, it will be added if not enough for tau

	FILE *fq, *fx ,*fy, *fz;
	int Dd;
	Dd= (int)(10*D);
	char num[20];
	char xlocal[20];
	char ylocal[20];
	char zlocal[20];
	sprintf(num,"num%d.dat",Dd);
	sprintf(xlocal,"x%d.dat",Dd);
	sprintf(ylocal,"y%d.dat",Dd);
	sprintf(zlocal,"z%d.dat",Dd);	
	fq=fopen(num,"wb+"); //record number
	fx=fopen(xlocal,"wb+");//record the last location
	fy=fopen(ylocal,"wb+");
	fz=fopen(zlocal,"wb+");
	
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
//for (i=0;i<nion;i++){printf("%f\n",q2[i]);}
	////put ions randomly, i.e.ions' three-dimensional values are random
	srand((unsigned)time(NULL));
	
	double initialE;	
	////place ions
	for (i=0;i<nall;i++){

		judgement(i,x1,y1,z1,x2,y2,z2);  //function call
		
	}
 
	printf("Random finished!\n");	
	

	
	////move many steps to reach equilibrium
	srand((unsigned)time(NULL));
	double step1;
	double nion_step=0;
	for (step1=1.0;step1<allstep+1;step1++)
	{
		
		////move ions
		double sample_diff;
		sample_diff=rand()%10;
		if (sample_diff<2.5) {i=(int) (rand()%(nion)); nion_step++;}
		else {i=(int) (rand()%(nc)+nion);}
		//printf("No.%d\n",i); 
		int c;
		c=move(i,x1,y1,z1,q1,x2,y2,z2,q2);  //function call

		//print the result
		//every bin print once
		//if (4==c){
			print_number(fq,x2,y2,D);  //function call
			//print_curve(fp,x2,y2);   //function call
			//printf("%d\n",step1);
		//}
		
	if (0==fmod(step1,1000.0))	{printf("%f\n",step1);}
					
	}
	printf("Moving finished!\n");
	printf("couacc=%d\n",countaccept);
	printf("countnc=%d\n",countnc);
	
	print_location(fx,fy,fz,x2,y2,z2);
	
	fclose(fq);
	fclose(fx);
	fclose(fy);
	fclose(fz);

	double acco=countaccept; 
	//printf("countaccept=%d\n",acco);
	printf("acceptance ni=%f\n",acco/(double) nion_step);
	printf("acceptance nc=%f\n",countnc/(double) (allstep-nion_step));

		 
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
			if (i>=(n1+n2+n3+n4)){
				if (d_im<(r0+rc)){
					judgement(i,x1,y1,z1,x2,y2,z2);}
			}
			else {
		 		if (d_im<(ri+r0)){
					judgement(i,x1,y1,z1,x2,y2,z2);}  //function call
			}
		}
		////between ions
		int k;
		double d_ii;
		for (k=0;k<i;k++){
			d_ii=distance(x2[k],y2[k],z2[k],x2[i],y2[i],z2[i]);
			if (i<(n1+n2+n3+n4)){
				if (d_ii<(ri+ri)){
					judgement(i,x1,y1,z1,x2,y2,z2);}
			}
			else {
				if(k<(n1+n2+n3+n4)){if (d_ii<(ri+rc)) judgement(i,x1,y1,z1,x2,y2,z2);}
				else {if(d_ii<(rc+rc)) judgement(i,x1,y1,z1,x2,y2,z2);}  //function call
			}
		}
		////update ith ion xyz
}


//***************************************************************************

 
 
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
	for (k=0;k<nall;k++){
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

////system-energy function
double E(int i,double *x1,double *y1,double *z1,double *q1,double *x2,double *y2,double *z2,double *q2)  //function definition
{
	double u=0.0;
	////the ion and monomers
	int j;
	double d_im;
	for (j=0;j<N+N;j++){
		d_im=distance(x1[j],y1[j],z1[j],x2[i],y2[i],z2[i]);
		//printf("dim=%f\n",d_im);
		//printf("q1=%f\n",q1[j]);
		if (d_im<R){
			u+=lb*q2[i]*q1[j]/d_im+10.0/3.0*power6(1.0/d_im)*power6(1.0/d_im)-10.0/3.0*power6(1.0/d_im)+5/6;
			//printf("u1=%f\t",u);
		}
		else {
			u+=lb*q2[i]*q1[j]/d_im;
			//printf("u2=%f\t",u);
		}
	}
	////the ion and other ions
	int k;
	double d_ii;
	for (k=0;k<nion;k++)
	{
		d_ii=distance(x2[k],y2[k],z2[k],x2[i],y2[i],z2[i]);
		//printf("dii=%f\n",d_ii);
		if (i-k){
			if (d_ii<R){
				u+=lb*q2[i]*q2[k]/d_ii+10.0/3.0*power6(1.0/d_ii)*power6(1.0/d_ii)-10.0/3.0*power6(1.0/d_ii)+5/6;   //no need to be divided by 2 since only calculate 1 ion with others
				//printf("u3=%f\t",u);
			}
			else {
				u+=lb*q2[i]*q2[k]/d_ii;    //no need to be divided by 2 since only calculate 1 ion with others
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
	double d_im;
	for (j=0;j<N+N;j++){
		d_im=distance(x1[j],y1[j],z1[j],xtemp,ytemp,ztemp);
		//printf("dim=%f\n",d_im);
		//printf("q1=%f\n",q1[j]);
		if (d_im<R){
			u+=lb*q2[i]*q1[j]/d_im+10.0/3.0*power6(1.0/d_im)*power6(1.0/d_im)-10.0/3.0*power6(1.0/d_im)+5/6;
			//printf("u1=%f\t",u);
		}
		else {
			u+=lb*q2[i]*q1[j]/d_im;
			//printf("u2=%f\t",u);
		}
	}
	////the ion and other ions
	int k;
	double d_ii;
	for (k=0;k<nion;k++)
	{
		d_ii=distance(x2[k],y2[k],z2[k],xtemp,ytemp,ztemp);
		//printf("dii=%f\n",d_ii);
		if (i-k){
			if (d_ii<R){
				u+=lb*q2[i]*q2[k]/d_ii+10.0/3.0*power6(1.0/d_ii)*power6(1.0/d_ii)-10.0/3.0*power6(1.0/d_ii)+5/6;  //no need to be divided by 2 since only calculate 1 ion with others
				//printf("u3=%f\t",u);
			}
			else {
				u+=lb*q2[i]*q2[k]/d_ii;   //no need to be divided by 2 since only calculate 1 ion with others
				//printf("u4=%f\t",u);
			}
		}
	}
	//printf("utemp=%f\n",u);
	return u;
}

//*************************************************************************

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


////---------------------------------------------------------------------------------------my_functions










/* 
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


