#include "sadrib.h"
#include "random.h"
#include <stdio.h>
#include <stdlib.h>

Ribbon *ribbon;

Vector *tt;
Vector *nn;
Vector *bb;
Vector *rr;

Statistics stat;

long seed;

double temperature;
double force;
double torque;
double *f_tt;
double *f_bb;

double *cos_psi;
double *cos_tau;
double *sin_psi;
double *sin_tau;

double tm[3][3];

double get_zeta();
double get_energy();
double get_torsion(); 
double get_curvature();
double get_end_to_end();
double dot(Vector *, Vector *);
double get_link_n(); // linking //
double get_twist();
double get_writhe();
Vector cproduct(Vector *, Vector *);

Int_Type num_of_iteration;
Int_Type num_of_particles;
 
void end();
void run();
void init();
void debug();
void frenet();
void measure();
void init_spiral();
void init_random();
void write_shape();
void allocate_memory();
void read_conf(char *);
void write_conf(FILE *);
void one_step(Int_Type);
void write_hi(FILE *, Int_Type);
void get_correlation_functions(double *, double *);

/*******************************************************************/

main() 
{
	init();
	run();
	end();
	return EXIT_SUCCESS;
}

/*******************************************************************/
void init()
{
	char f_name[32];
	Int_Type cflag;
 
   FILE *fptr;

   if ((fptr = fopen("Input.txt","r")) == NULL){
       printf("Error! opening file");

       // Program exits if the file pointer returns NULL.
       exit(1);
   }

   fscanf(fptr,"%u%u%lg%lg%lg%ld%u", &num_of_particles, &num_of_iteration, &temperature, &force, &torque, &seed, &cflag);

   fclose(fptr); 

	stat.ar = 0;	
	stat.zeta = 0;
	stat.link = 0; // linking 
	stat.twist = 0; //twist
	stat.writhe = 0; // writhe
	stat.energy = 0;
	stat.torsion = 0;
	stat.curvature = 0;
	stat.end_to_end = 0;
	stat.link1 = 0; // linking mtd 2
	stat.twist1 = 0; // twist mtd 2
	stat.writhe1 = 0; // writhe mtd 2
	
	tm[0][0] = 0;
	tm[0][1] = 0;
	tm[0][2] = 0;
	tm[1][0] = 0;
	tm[1][1] = 0;
	tm[1][2] = 0;
	tm[2][0] = 0;
	tm[2][1] = 0;
	tm[2][2] = 0;

  	switch(cflag){
  	case 1:
    	allocate_memory();
    	init_spiral();
    	break;
  	case 2:
    	allocate_memory();
    	init_random();
    	break;
  	case 3:
    	read_conf(f_name);
    	break;
  	}
}

/*******************************************************************/

void allocate_memory()
{
  	Int_Type i;

	ribbon = (Ribbon *) malloc(num_of_particles*sizeof(Ribbon));

	cos_psi = calloc(num_of_particles,sizeof(double));
	cos_tau = calloc(num_of_particles,sizeof(double));
	sin_psi = calloc(num_of_particles,sizeof(double));
	sin_tau = calloc(num_of_particles,sizeof(double));
	
	f_tt = (double *) malloc(num_of_particles*sizeof(double));
	f_bb = (double *) malloc(num_of_particles*sizeof(double));

	tt = (Vector *) malloc(num_of_particles*sizeof(Vector));	
	nn = (Vector *) malloc(num_of_particles*sizeof(Vector));
	bb = (Vector *) malloc(num_of_particles*sizeof(Vector));

	rr = (Vector *) malloc(num_of_particles*sizeof(double));

	stat.tt = (Vector *) calloc(num_of_particles,sizeof(Vector));
	stat.nn = (Vector *) calloc(num_of_particles,sizeof(Vector));
	stat.bb = (Vector *) calloc(num_of_particles,sizeof(Vector));
	
	stat.f_tt = calloc(num_of_particles,sizeof(double));
	stat.f_bb = calloc(num_of_particles,sizeof(double));
}

/*******************************************************************/

void init_random()
{
  double psi, tau;
  Int_Type i;

	for (i = 1; i < num_of_particles-1; i++){
		psi = acos(1-2*ran2(&seed));
		tau = PI*(1-2*ran2(&seed))/2.0;
		ribbon[i-1].psi = psi;
		ribbon[i-1].tau = tau;
	}
}

void init_spiral()
{
	double x, y, psi, tau, curv, tors;
  	Int_Type i;
   
    FILE *fptr;

   if ((fptr = fopen("CT.txt","r")) == NULL){
       printf("Error! opening file");

       // Program exits if the file pointer returns NULL.
       exit(1);
   }

   fscanf(fptr,"%lg%lg", &curv, &tors);

   fclose(fptr); 

	psi = acos(1-0.5*curv*curv);
	tau = acos(1-0.5*tors*tors);

  	for (i=1; i<num_of_particles-1; i++){
    	ribbon[i-1].psi = psi;
    	ribbon[i-1].tau = tau;
 	}
}

void read_conf(char *f_name)
{
	Int_Type i, n = 0;    
	double psi, tau;
	
	FILE *f_in;
	
	if ((f_in = fopen(f_name,"r")) == NULL){
		printf("Error: file %s not found\n",f_name);
		exit(0);
	}

	printf("Reading %s...\n",f_name);
	
	while (fscanf(f_in,"%*lg%*lg")!=EOF) n++;
	
	fseek(f_in,0,SEEK_SET);
	num_of_particles = n;
	allocate_memory();

	for (i=0; i<num_of_particles-2; i++){
		fscanf(f_in,"%lg%lg",&psi,&tau);
		ribbon[i].psi = psi;
		ribbon[i].tau = tau;
	}
	fclose(f_in);
}

void frenet()
{
	double sp, st, cp, ct, psi, tau;
	Int_Type i;
	
	tt[0].x = 0.;	
	tt[0].y = 0.;
	tt[0].z = 1.;
	
	nn[0].x = 1.;
	nn[0].y = 0.;
	nn[0].z = 0.;
	
	bb[0].x = 0.;
	bb[0].y = 1.;
	bb[0].z = 0.;

	rr[0].x = 0.0;
	rr[0].y = 0.0;
	rr[0].z = 0.0;
	
	for (i=1; i<num_of_particles; i++){
		
		psi = ribbon[i-1].psi;
		tau = ribbon[i-1].tau;

		sp = sin(psi);
		st = sin(tau);
		cp = cos(psi);
		ct = cos(tau);
		
		tt[i].x = cp*tt[i-1].x + sp*ct*nn[i-1].x + sp*st*bb[i-1].x;
		tt[i].y = cp*tt[i-1].y + sp*ct*nn[i-1].y + sp*st*bb[i-1].y;
		tt[i].z = cp*tt[i-1].z + sp*ct*nn[i-1].z + sp*st*bb[i-1].z;
		
		nn[i].x =-sp*tt[i-1].x + cp*ct*nn[i-1].x + cp*st*bb[i-1].x;
		nn[i].y =-sp*tt[i-1].y + cp*ct*nn[i-1].y + cp*st*bb[i-1].y;
		nn[i].z =-sp*tt[i-1].z + cp*ct*nn[i-1].z + cp*st*bb[i-1].z;
		
		bb[i].x =-st*nn[i-1].x + ct*bb[i-1].x;
		bb[i].y =-st*nn[i-1].y + ct*bb[i-1].y;
		bb[i].z =-st*nn[i-1].z + ct*bb[i-1].z;
		
		rr[i].x  = rr[i-1].x+tt[i-1].x;
		rr[i].y  = rr[i-1].y+tt[i-1].y;
		rr[i].z  = rr[i-1].z+tt[i-1].z;
	}	
}

Vector cproduct(Vector *v1, Vector *v2){
	
	Vector v;
	v.x =v1->y*v2->z-v2->y*v1->z;
	v.y =-v1->x*v2->z+v2->x*v1->z;
	v.z =v1->x*v2->y-v2->x*v1->y;
	return v; 
}

double get_energy(){
	
	double cp, ct, en=0;
	Int_Type i;

	for (i=1; i<num_of_particles-1; i++){
		cp = cos(ribbon[i-1].psi);
		ct = cos(ribbon[i-1].tau);
		en+= (2-cp-ct)*(2-cp-ct)/(1-cp); // Sadowsky energy
	}
	
	if (force){
		en -= force*get_zeta();
	}
	
	// linking //
	if (torque){
		en -= 2*PI*torque*get_link_n();
	}
	
	return en;
}

double get_link_n(){
	Vector dummy,r1,r2,r3,r4,r13,r14,r23,r24,n1,n2,n3,n4;
	double Omega =0;
	double norm1, norm2, norm3, norm4, ti, lk = 0, Wr = 0, Tw = 0;
	Int_Type i,j;
	
	frenet();

	ti = -(bb[1].x-bb[0].x)*nn[1].x -(bb[1].y-bb[0].y)*nn[1].y -(bb[1].z-bb[0].z)*nn[1].z;
	Tw += ti;

	for (i=2; i<num_of_particles-1; i++){

		ti = -(bb[i].x-bb[i-1].x)*nn[i].x -(bb[i].y-bb[i-1].y)*nn[i].y -(bb[i].z-bb[i-1].z)*nn[i].z;
		Tw += ti;

		for (j=1; j<i;j++){
			if ((i-j)>1){
				r1=rr[i];
				r2=rr[i+1];
				r3=rr[j];
				r4=rr[j+1];

				dummy=cproduct(&tt[j],&tt[i]);

				r13.x=r3.x-r1.x;
				r13.y=r3.y-r1.y;
				r13.z=r3.z-r1.z;

				r14.x=r4.x-r1.x;
				r14.y=r4.y-r1.y;
				r14.z=r4.z-r1.z;

				r23.x=r3.x-r2.x;
				r23.y=r3.y-r2.y;
				r23.z=r3.z-r2.z;

				r24.x=r4.x-r2.x;
				r24.y=r4.y-r2.y;
				r24.z=r4.z-r2.z;

				n1=cproduct(&r13,&r14);
				n2=cproduct(&r14,&r24);
				n3=cproduct(&r24,&r23);
				n4=cproduct(&r23,&r13);
				
				norm1=sqrt(dot(&n1,&n1));
				norm2=sqrt(dot(&n2,&n2));
				norm3=sqrt(dot(&n3,&n3));
				norm4=sqrt(dot(&n4,&n4));
			
				n1.x=n1.x/norm1;
				n1.y=n1.y/norm1;
				n1.z=n1.z/norm1;

				n2.x=n2.x/norm2;
				n2.y=n2.y/norm2;
				n2.z=n2.z/norm2;

				n3.x=n3.x/norm3;
				n3.y=n3.y/norm3;
				n3.z=n3.z/norm3;

				n4.x=n4.x/norm4;
				n4.y=n4.y/norm4;
				n4.z=n4.z/norm4;
			
				Omega += copysign(1.0, dot(&dummy,&r13))*(asin(dot(&n1,&n2))+ asin(dot(&n2,&n3))+ asin(dot(&n3,&n4))+ asin(dot(&n4,&n1)));
				}
			}	
	}

	Wr = Omega;
	lk = Tw + Wr;
	lk /= (2.0*PI);
	
	return lk;	
}

double get_twist(){
	
	double ti;
	double Tw=0;
	Int_Type i;
	
	frenet();
	
	for (i=1; i<num_of_particles; i++){
		
		ti = -(bb[i].x-bb[i-1].x)*nn[i].x -(bb[i].y-bb[i-1].y)*nn[i].y -(bb[i].z-bb[i-1].z)*nn[i].z;
		Tw += ti;
	}
	
	Tw = Tw/(2*PI);
	
	return Tw;

}

double get_writhe(){
	
	double ki;
	double Wr=0;
	Int_Type i;
	
	frenet();
	
	for (i=1; i<num_of_particles; i++){
		
		ki =  (tt[i].x-tt[i-1].x)*nn[i].x +(tt[i].y-tt[i-1].y)*nn[i].y +(tt[i].z-tt[i-1].z)*nn[i].z;
		Wr += ki*bb[i].z/(1+ tt[i].z);
		
	}
	
	Wr /= (2*PI);
	
	return Wr;
	
}

double dot(Vector *v1, Vector *v2){
	
	return v1->x*v2->x+v1->y*v2->y+v1->z*v2->z; 
}

double get_end_to_end()
{
	Int_Type i;
	Vector pt;
	double r2;
	
	frenet();
	
	pt.x = 0;
	pt.y = 0;
	pt.z = 0;
	
	for (i=0; i<num_of_particles-1; i++){
		pt.x += tt[i].x;
		pt.y += tt[i].y;
		pt.z += tt[i].z;
	}
	
	r2 = pt.x*pt.x+pt.y*pt.y+pt.z*pt.z;
	return r2;
}

/*******************************************************************/

double get_zeta()
{
	double z=0;
	Int_Type i;
	
	frenet();
	
	for (i=0; i<num_of_particles-1; i++){
		z += tt[i].z;
	}
	
	return z;
}


double get_curvature()
{
	double curvature=0;
	double ki;
	Int_Type i;
	
	frenet();
	
	for (i=0; i<num_of_particles-2; i++){
		curvature += 2*(1-cos(ribbon[i].psi));
		
	}
	
	return curvature/(double) num_of_particles;
}

double get_torsion()
{
	double torsion=0;
	Int_Type i;
	
	for (i=0; i<num_of_particles-2; i++){
		torsion += 2*(1-cos(ribbon[i].tau));
	}
	
	//return torsion;
	return torsion/(double) num_of_particles;
}

void get_cf(double *f1, double *f2){
	
	int i, j, *nu;
	
	frenet();
	
	nu = calloc(num_of_particles, sizeof(int));

	for (i=0; i<num_of_particles-1; i++){
		f1[i] = 0.;
		f2[i] = 0.;
	}
	
	for (i=0; i<num_of_particles-1; i++){
		for(j=0; j<num_of_particles-1; j++){
			f1[abs(i-j)]+= dot(&tt[i],&tt[j]);
			f2[abs(i-j)]+= dot(&bb[i],&bb[j]);
			nu[abs(i-j)]++;
		}
	}
	
	for (i=0; i<num_of_particles-1; i++){
			f1[i] /= (double) nu[i];
			f2[i] /= (double) nu[i];
	}	
	free(nu);
}

void one_step(Int_Type i)
{
	double e1, e2, z1, z2, cp, ct, psi, tau, psi_old, tau_old, boltzmann, lk1, lk2;
		
	psi = acos(1-2*ran2(&seed));
	tau = PI*(1-2*ran2(&seed))/2.0;
	
	cp = cos(ribbon[i].psi);
	ct = cos(ribbon[i].tau);


	e1 = (2-cp-ct)*(2-cp-ct)/(1-cp);
	e2 = (2-cos(psi)-cos(tau))*(2-cos(psi)-cos(tau))/(1-cos(psi));

	if (force){
		psi_old = ribbon[i].psi;
		tau_old = ribbon[i].tau;
		z1 = get_zeta();
		ribbon[i].psi = psi;
		ribbon[i].tau = tau;
		z2 = get_zeta();
		e1 -= force*z1;
		e2 -= force*z2;
		ribbon[i].psi = psi_old;
		ribbon[i].tau = tau_old;
	}

	if (torque){
		psi_old = ribbon[i].psi;
		tau_old = ribbon[i].tau;
		lk1 = get_link_n();
		ribbon[i].psi = psi;
		ribbon[i].tau = tau;
		lk2 = get_link_n();
		e1 -= 2*PI*torque*lk1;
		e2 -= 2*PI*torque*lk2;
		ribbon[i].psi = psi_old;
		ribbon[i].tau = tau_old;
	}

	boltzmann = exp((e1-e2)/temperature);	
	
	if (e2<e1){
		ribbon[i].psi = psi;
		ribbon[i].tau = tau;
		stat.ar++;
		return;
	}
	else if (ran2(&seed)<boltzmann){
		ribbon[i].psi = psi;	
		ribbon[i].tau = tau;	
		stat.ar++;	
	}
}

void measure()
{
	double w, psi, tau;	
	Int_Type i;
	
	get_cf(f_tt,f_bb);
	
	w = 2./num_of_iteration;

	psi = ribbon[0].psi;
	tau = ribbon[0].tau;

	tm[0][0] += w*cos(psi);
	tm[0][1] +=-w*sin(psi);
	tm[0][2] += 0;
	tm[1][0] += w*sin(psi)*cos(tau);
	tm[1][1] += w*cos(psi)*cos(tau);
	tm[1][2] +=-w*sin(tau);
	tm[2][0] += w*sin(psi)*sin(tau);
	tm[2][1] += w*cos(psi)*sin(tau);
	tm[2][2] += w*cos(tau);
	
	for (i=0; i<num_of_particles; i++){
		stat.tt[i].x += w*tt[i].x;
		stat.tt[i].y += w*tt[i].y;
		stat.tt[i].z += w*tt[i].z;
		stat.f_tt[i] += w*f_tt[i];
		stat.f_bb[i] += w*f_bb[i];
		cos_psi[i] += w*cos(ribbon[i].psi);
		cos_tau[i] += w*cos(ribbon[i].tau);
		sin_psi[i] += w*sin(ribbon[i].psi);
		sin_tau[i] += w*sin(ribbon[i].tau);
	}
	
	stat.end_to_end += w*get_end_to_end();	
	stat.curvature += w*get_curvature();
	stat.torsion += w*get_torsion();
	stat.energy += w*get_energy();
	stat.zeta += w*get_zeta();
	stat.link += w*get_link_n();
	stat.twist += w*get_twist();
	stat.writhe += w*get_writhe();
    stat.link1 += w*get_link1();
	stat.twist1 += w*get_twist1();
	stat.writhe1 += w*get_writhe1();
}

/*******************************************************************/

void run()
{
	Int_Type t, i, j, period=1000;
	
	FILE *f_hi, *f_abi, *f_ano;
	f_hi = fopen("hi.dat","w");			

	for (t=0; t<num_of_iteration; t++){
		for (i=0; i<num_of_particles-2; i++){
			one_step(i);			
		}
		if (t>num_of_iteration/2) measure();
		if((t+1)%period==0){
	      printf("i%u ",t);
	      write_hi(f_hi,t);
	      fflush(stdout);
	    }
	    if((t+1)%(100*period)==0){
	      printf("\n");
	      fflush(stdout);
	    }
	}

	fclose(f_hi);
}

void write_hi(FILE *f_ou, Int_Type t)
{
  	double energy, torsion, extension_z, link_n1, twist, writhe;

	energy = get_energy();
	torsion = get_torsion();
	extension_z= get_zeta();
	link_n1=get_link_n();
	twist=get_twist();
	writhe=get_writhe();

  	fprintf(f_ou,"%u\t%g\t%g\t%g\t%g\t%g\t%g\n", t, extension_z, link_n1, energy, torsion, twist, writhe);
}

void write_conf(FILE *f_ou)
{
  	double psi, tau;
	Int_Type i;
	
	for (i=0; i<num_of_particles-2; i++){
   		psi = fmod(ribbon[i].psi,2*PI);
   		tau = fmod(ribbon[i].tau,2*PI);
	  	fprintf(f_ou,"%g\t%g\n",psi,tau);
	}
}

void end()
{
	char color[32];
	double *cf1;
	Vector pt;
  	FILE *f_ou;

	Int_Type i, j;
	
	frenet();
	
	printf("\n");
	printf("\tEnergy %g\n",stat.energy);
	printf("\tTorsion %g\n",stat.torsion);
	printf("\tCurvature %g\n",stat.curvature);
	printf("\tEnd to end distance %g\n",stat.end_to_end);
	printf("\tExtension %g\n",stat.zeta);
	printf("\tLinking %g\n",stat.link); /* linking number */
	printf("\tTwist %g\n",stat.twist);		
	printf("\tWrithe %g\n",stat.writhe);
	printf("\tAcceptance ratio %g\n", stat.ar/(double)(num_of_iteration*num_of_particles));
	printf("\n");

	f_ou = fopen("en.dat","w");
	fprintf(f_ou,"Number of particles %u\n",num_of_particles);
	fprintf(f_ou,"Number of iteration %u\n",num_of_iteration);
	fprintf(f_ou,"Temperature %g\n",temperature);
	fprintf(f_ou,"Force %g\n",force);
	fprintf(f_ou,"Torque %g\n",torque); /* torque */
	fprintf(f_ou,"Energy %g\n",stat.energy);
	fprintf(f_ou,"Torsion %g\n",stat.torsion);
	fprintf(f_ou,"Curvature %g\n",stat.curvature);
	fprintf(f_ou,"End to end distance %g\n",stat.end_to_end);
	fprintf(f_ou,"Extension %g\n",stat.zeta);
	fprintf(f_ou,"Linking number %g\n",stat.link); /* linking number */
	fprintf(f_ou,"Twist number %g\n",stat.twist); /* twist number */
	fprintf(f_ou,"Writhe number %g\n",stat.writhe); /* writhe number */
	fprintf(f_ou,"Acceptance ratio %g\n", stat.ar/(double)(num_of_iteration*num_of_particles));
	fclose(f_ou);

	f_ou = fopen("liveplot.m","w");	
   	fprintf(f_ou,"Graphics3D[{{GrayLevel[0],{Line[{{0,0,0}");
	
	for (i=1; i<num_of_particles; i++){
		pt.x = 0;
		pt.y = 0;
		pt.z = 0;
		for (j=0; j<i; j++){
			pt.x += tt[j].x;
			pt.y += tt[j].y;
			pt.z += tt[j].z;
		}
		fprintf(f_ou,",{%g,%g,%g}",pt.x,pt.y,pt.z);
	}
	fprintf(f_ou,"}],Point[{0,0,0}]");
	for (i=1; i<num_of_particles; i++){
		pt.x = 0;
		pt.y = 0;
		pt.z = 0;
		for (j=0; j<i; j++){
			pt.x += tt[j].x;
			pt.y += tt[j].y;
			pt.z += tt[j].z;
		}
		fprintf(f_ou,",Point[{%g,%g,%g}]",pt.x,pt.y,pt.z);
	}
	fprintf(f_ou,"}}},Boxed->False]");
	fclose(f_ou);

	f_ou = fopen("cf.dat","w");
	for (i=0; i<num_of_particles-1; i++){
		fprintf(f_ou,"%u\t%g\t%g\n",i,stat.f_tt[i],stat.f_bb[i]);
	}
	fclose(f_ou);

	f_ou = fopen("anglebi.dat","w");
	for (i=1; i<num_of_particles-1; i++){
		fprintf(f_ou,"%u\t%g\n", i, acos(dot(&bb[i], &bb[i-1])));
	}
	fclose(f_ou);

	f_ou = fopen("angleno.dat","w");
	for (i=1; i<num_of_particles-1; i++){
		fprintf(f_ou,"%u\t%g\n", i, acos(dot(&nn[i], &nn[i-1])));
	}
	fclose(f_ou);
	
	f_ou = fopen("ou.dat","w");
	write_conf(f_ou);
	fclose(f_ou);
	
	f_ou = fopen("tm.dat","w");
	fprintf(f_ou,"%lg\t%lg\t%lg\n",tm[0][0],tm[0][1],tm[0][2]);
	fprintf(f_ou,"%lg\t%lg\t%lg\n",tm[1][0],tm[1][1],tm[1][2]);
	fprintf(f_ou,"%lg\t%lg\t%lg\n",tm[2][0],tm[2][1],tm[2][2]);
	fclose(f_ou);

}
