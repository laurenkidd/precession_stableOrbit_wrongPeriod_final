#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#define pi  3.14159
#define me   5.9736e24          //[kg]
#define ms   1.9891e30         //[kg]
#define	a   149598261e3      //[m]
#define e   0.0167112303531389  //eccentricity
#define G   6.67428e-11         //[N m^2 kg^-2]
#define Tday  86164.1            //[s] -> SIDEREAL DAY
#define	I1   8.008e37           //[kg m^2]
#define	I3    8.034e37 	    //[kg m^2]
#define theta0  23.45*pi/180.      //[radians]*/
#define Ex0  149598261e3*(1.0+0.0167112303531389) //a*(1+e) init dist from sun
#define	Ey0  0.0
#define	Ez0  0.0	
//kinetic energy
#define KE  (G*ms*me)/a * (1/(1+e) - 0.5)
//earth velocity
#define v0  sqrt((2*KE)/(me*(1+me/ms)))	
#define vex0  0.0
#define	vey0  v0
#define	vez0  0.0
//euler angles
#define phi0  0.0		//phi_0
#define	psi0  0.0		//psi_0
#define	vphi0  0.0	//dot(phi_0)
#define	vtheta0  0.0	//dot(theta_0)
#define	vpsi0  2.0*pi/Tday //dot(psi_0)
//initial conditions for sun
#define	Sx0  0.0 //let sun start at origin
#define	Sy0  0.0
#define	Sz0  0.0
#define	vsx0  0.0
#define	vsy0  -(me/ms)*v0
#define	vsz0  0.0
#define c0  2.0*pi/Tday






double functionx(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);
double (*funcx)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);
double functiony(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);
double (*funcy)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);
double functiony(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);
double (*funcy)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);
double functionz(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);
double (*funcz)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t);

double functiontheta(double Ex, double Ey, double Ez, double Sx,double Sy,  double Sz, double theta, double vtheta, double phi, double vphi);
double (*funct)(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi);
double functionphi(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi);
double (*funcph)(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi);


//double functionS(double x, double v, double t);
//double (*funcS)(double x, double v, double t);

double rkn4(double (*funcx)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t),
double (*funcy)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t),
 double (*funcz)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t),
 double (*funct)(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi),
double (*funcph)(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi),
double Ex[], double vex[],double Ey[],  double vey[], double Ez[], double vez[], double Sx[],
double vsx[],  double Sy[], double vsy[], double Sz[], double vsz[], double phi[], double vphi[],  double theta[], double vtheta[],  double psi[], double t, double h, long steps);

double functionx(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t){
	//	theta = theta0;
	Ez = 0;
	Sz = 0;
	
	//x position in CM-centered frame
	
	double d = sqrt((Sx-Ex)*(Sx-Ex)+ (Sy-Ey)*(Sy-Ey)+(Sz-Ez)*(Sz-Ez));
	double dx, dy, dz, d3;
	dx = Sx - Ex;
	dy = Sy - Ey;
	dz = Sz - Ez;
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	//printf("%f  %f %f %f %f\n", phi, theta, dy, dz, d3);
	double fx;
	fx = (-(G*me*ms)/(d*d*d))*dx + ((3*G*ms)/pow(d,5))*(I1-I3)*(dx/2. - (2.5*d3*d3*dx)/(d*d) + d3*(sin(phi)*sin(theta)));
	return -fx/me;
        
}

double functiony(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t){
	//theta = theta0;
	Ez = 0;
	Sz = 0;

	
	//y position in CM-centered frame
	
	double d = sqrt((Sx-Ex)*(Sx-Ex)+ (Sy-Ey)*(Sy-Ey)+(Sz-Ez)*(Sz-Ez));
	double dx, dy, dz, d3;
	dx = Sx - Ex;
	dy = Sy - Ey;
	dz = Sz - Ez;
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	double fy;
	fy = (-(G*me*ms)/(d*d*d))*dy + ((3*G*ms)/pow(d,5))*(I1-I3)*(dy/2. - (2.5*d3*d3*dy)/(d*d) + d3*(-cos(phi)*sin(theta)));
	return -fy/me;
}

double functionz(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t){
	//theta = theta0;
	Ez = 0;
	Sz = 0;
	
	//z position in CM-centered frame
	
	double d = sqrt((Sx-Ex)*(Sx-Ex)+ (Sy-Ey)*(Sy-Ey)+(Sz-Ez)*(Sz-Ez));
	double dx, dy, dz, d3;
	dx = Sx - Ex;
	dy = Sy - Ey;
	dz = Sz - Ez;
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	double fz;
	fz = (-(G*me*ms)/(d*d*d))*dz + ((3.0*G*ms)/pow(d,5))*(I1-I3)*(dz/2. - (2.5*d3*d3*dz)/(d*d) + d3*(cos(theta)));
	return -fz/me;
}

double functionphi(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi){
	//theta = theta0;
	//vtheta = 0;
	Ez = 0;
	Sz = 0;
	
	double d = sqrt((Sx-Ex)*(Sx-Ex)+ (Sy-Ey)*(Sy-Ey)+(Sz-Ez)*(Sz-Ez));
	double dx, dy, dz, d3;
	dx = Sx - Ex;
	dy = Sy - Ey;
	dz = Sz - Ez;
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	//c0 =  2*pi/Tday;
	//return (1.0/sin(theta))*(-2.0*vphi*vtheta*cos(theta)+(I3/I1)*c0*vtheta+((3*G*ms)/pow(d,5))*((I1-I3)/I1)*d3*(dx*cos(phi)+dy*sin(phi)));
	//going to assume vtheta is very small and neglact any dependance on it... also to see if this makes the angles stop blowing up
	return (1.0/sin(theta))*(((3*G*ms)/pow(d,5))*((I1-I3)/I1)*d3*(dx*cos(phi)+dy*sin(phi)));

}

double functiontheta(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi){
	//theta = theta0;
	//vtheta=0;
	Ez = 0;
	Sz = 0;
	
	double d = sqrt((Sx-Ex)*(Sx-Ex)+ (Sy-Ey)*(Sy-Ey)+(Sz-Ez)*(Sz-Ez));
	double dx, dy, dz, d3;
	dx = Sx - Ex;
	dy = Sy - Ey;
	dz = Sz - Ez;
	d3 = (dx*sin(phi) - dy*cos(phi))*sin(theta) + dz*cos(theta);
	//c0 =  2*pi/Tday;
	return vphi*vphi*sin(theta)*cos(theta)-(I3/I1)*c0*sin(theta)+(3*G*ms)/pow(d,5)*((I1-I3)/I1)*d3*((dx*sin(phi)-dy*cos(phi))*cos(theta)-dz*sin(theta));
}


int main() {

	


	long steps = 100000; 		//~3 yrs right now		//Number of steps.
	double t=0,h=1000; 			//time, step size(sec).
	//define arrays to be used
	double *Ex = malloc(steps * sizeof(double)); 	//Position array, x, earth
        double *Sx = malloc(steps * sizeof(double)); 	//position array, x, sun
	double *vex = malloc(steps * sizeof(double)); 	//velocity array, x, earth
        double *vsx = malloc(steps * sizeof(double)); 	//Velocity array, x sun

	double *Ey = malloc(steps * sizeof(double)); 	//Position array, y, earth
        double *Sy = malloc(steps * sizeof(double)); 	//position array, y, sun
	double *vey = malloc(steps * sizeof(double)); 	//velocity array, y, earth
        double *vsy = malloc(steps * sizeof(double)); 	//Velocity array, y, sun
	
	double *Ez = malloc(steps * sizeof(double)); 	//Position array, z, earth
        double *Sz = malloc(steps * sizeof(double)); 	//position array, z, sun
	double *vez = malloc(steps * sizeof(double)); 	//velocity array, z, earth
        double *vsz = malloc(steps * sizeof(double)); 	//Velocity array, z sun
	
	double *theta = malloc(steps * sizeof(double)); 	
        double *phi = malloc(steps * sizeof(double)); 	
	double *psi = malloc(steps * sizeof(double)); 
        double *vtheta = malloc(steps * sizeof(double)); 
	double *vphi = malloc(steps * sizeof(double)); 	



// Set pointer to point at the function to integrate. 
	funcx = functionx;
	funcy = functiony;
	funcz = functionz;
	funct = functiontheta;
	funcph = functionphi;

// Do integration. 


	rkn4(funcx,funcy,funcz,funct,funcph, Ex,vex,Ey, vey, Ez,  vez,  Sx,  vsx, Sy,  vsy, Sz, vsz, phi, vphi,  theta,
		  vtheta,  psi, t,h,steps);

	
// Print results to STDOUT 
	
	long int i;
        for ( i=0; i<steps; ++i){
		t += h;
		double d = sqrt((Sx[i]-Ex[i])*(Sx[i]-Ex[i])+ (Sy[i]-Ey[i])*(Sy[i]-Ey[i])+(Sz[i]-Ez[i])*(Sz[i]-Ez[i]));
		double AU = 1.5e11;
		double dmod = d/AU;
		if (i%1000 ==0){printf(" %e %f %e %e %f %f %f\n",dmod, t,Ex[i],Ey[i], phi[i], theta[i], psi[i]);}
        }



	free(Ex);
	free(Sx);
	free(vex);
	free(vsx);
	free(Ey);
	free(Sy);
	free(vey);
	free(vsy);
	free(Ez);
	free(Sz);
	free(vez);
	free(vsz);
	free(theta);
	free(phi);
	free(psi);
	free(vtheta);
	free(vphi);

	
	return 0;
}



	//runge kutta nystrom 4

	double rkn4(double (*funcx)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t),
double (*funcy)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t),
 double (*funcz)(double Ex, double vex, double Ey, double vey, double Ez, double vez, double Sx, double Sy, double Sz, double phi, double theta, double t),
 double (*funct)(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi),
double (*funcph)(double Ex, double Ey, double Ez, double Sx, double Sy, double Sz, double theta, double vtheta, double phi, double vphi),
double Ex[], double vex[],double Ey[],  double vey[], double Ez[], double vez[], double Sx[],
double vsx[],  double Sy[], double vsy[], double Sz[], double vsz[], double phi[], double vphi[],  double theta[], double vtheta[],  double psi[], double t, double h, long steps){
	double k1x,k2x,k3x,k4x, k1y, k2y, k3y, k4y, k1z, k2z, k3z, k4z;
	double k1t, k2t, k3t, k4t, k1ph, k2ph, k3ph, k4ph, k1sx, k1sy, k1sz, k2sx, k2sy, k2sz, k3sx, k3sy, k3sz, k4sx, k4sy, k4sz;
	
	
	Ex[0] = Ex0; //Initial position
	vex[0] = vex0; //Initial velocity
	Ey[0] = Ey0; //Initial position
	vey[0] = vey0; //Initial velocity
	Ez[0] = Ez0;
	vez[0] = vez0;
	Sx[0] = Sx0;
	Sy[0] = Sy0;
	Sz[0] = Sz0;
	theta[0] = theta0;
	phi[0] = phi0;
	vtheta[0] = vtheta0;
	vphi[0] = vphi0;
	vsx[0] = -vex0*me/ms;
	vsy[0] = -vey0*me/ms;
	vsz[0] = -vez0*me/ms;


	
	long i;
	for ( i=1; i<steps; ++i){
		//double vsx, vsy, vsz;
		
		//1//

		k1x = funcx(Ex[i-1],vex[i-1],Ey[i-1], vey[i-1], Ez[i-1], vez[i-1], Sx[i-1], Sy[i-1], Sz[i-1], phi[i-1], theta[i-1], t);
		//printf("%f",k1x);
		k1sx = -k1x*me/ms;
		k1y = funcy(Ex[i-1],vex[i-1],Ey[i-1], vey[i-1], Ez[i-1], vez[i-1], Sx[i-1], Sy[i-1], Sz[i-1], phi[i-1], theta[i-1], t);
		k1sy = -k1y*me/ms;
		k1z = funcx(Ex[i-1],vex[i-1],Ey[i-1], vey[i-1], Ez[i-1], vez[i-1], Sx[i-1], Sy[i-1], Sz[i-1], phi[i-1], theta[i-1], t);
		k1sz = -k1z*me/ms;
		k1t = funct(Ex[i-1], Ey[i-1], Ez[i-1], Sx[i-1], Sy[i-1], Sz[i-1], theta[i-1], vtheta[i-1], phi[i-1], vphi[i-1] );
		k1ph = funcph(Ex[i-1], Ey[i-1], Ez[i-1], Sx[i-1], Sy[i-1], Sz[i-1], theta[i-1], vtheta[i-1], phi[i-1], vphi[i-1]);		

		//2//

		k2x = funcx(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k1z/2., 				Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,Sy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,theta[i-1]+
			h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k2sx = -k2x*me/ms;
		k2y = funcy(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k1z/2., 				Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,Sy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,theta[i-1]+
			h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k2sy = -k2y*me/ms;
		k2z = funcz(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k1x/2.,Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k1y/2.,Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k1z/2., 				Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,Sy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,theta[i-1]+
			h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k2sz = -k2z*me/ms;		
		k2t = funct(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						Sy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k1t/2.,  phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., 				vphi[i-1]+h*k1ph/2. );
		k2ph = funcph(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						Sy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k1t/2.,  phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., 				vphi[i-1]+h*k1ph/2. );

		//3//
	
                k3x = funcx(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k2z/2., 				Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,Sy[i-1]+h*vsy[i-1]/2.+h*h*k2sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,
			theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k3sx = -k3x*me/ms;
                k3y = funcy(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k2z/2., 				Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,Sy[i-1]+h*vsy[i-1]/2.+h*h*k2sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,
			theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k3sy = -k3y*me/ms;
		k3z = funcz(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8.,vex[i-1]+h*k2x/2.,Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8.,vey[i-1]+h*k2y/2.,Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8.,vez[i-1]+h*k2z/2., 				Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8.,Sy[i-1]+h*vsy[i-1]/2.+h*h*k2sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8., phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8.,
			theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., t+h/2.);
		k3sz = -k3z*me/ms;
		k3t = funct(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						Sy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k2t/2.,  
			phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., vphi[i-1]+h*k2ph/2. );
		k3ph = funcph(Ex[i-1]+h*vex[i-1]/2.+h*h*k1x/8., Ey[i-1]+h*vey[i-1]/2.+h*h*k1y/8., Ez[i-1]+h*vez[i-1]/2.+h*h*k1z/8., Sx[i-1]+h*vsx[i-1]/2.+h*h*k1sx/8., 						Sy[i-1]+h*vsy[i-1]/2.+h*h*k1sy/8., Sz[i-1]+h*vsz[i-1]/2.+h*h*k1sz/8.,theta[i-1]+h*vtheta[i-1]/2.+h*h*k1t/8., vtheta[i-1]+h*k2t/2.,  
			phi[i-1]+h*vphi[i-1]/2.+h*h*k1ph/8., vphi[i-1]+h*k2ph/2. );
	

		//4//  
		k4x = funcx(Ex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, Ey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y, Ez[i-1]+h*vez[i-1]+h*h*k3z/2.,vez[i-1]+h*k3z,
			Sx[i-1]+h*vsx[i-1]+h*h*k3sx/2., Sy[i-1]+h*vsy[i-1]+h*h*k3sy/2., Sz[i-1]+h*vsz[i-1]+h*h*k3sz/2., phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 
			theta[i-1]+h*vtheta[i-1]+h*h*k3t/2., t+h);
		k4sx = -k4x*me/ms;
		k4y = funcy(Ex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, Ey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y, Ez[i-1]+h*vez[i-1]+h*h*k3z/2.,vez[i-1]+h*k3z,
			Sx[i-1]+h*vsx[i-1]+h*h*k3sx/2., Sy[i-1]+h*vsy[i-1]+h*h*k3sy/2., Sz[i-1]+h*vsz[i-1]+h*h*k3sz/2., phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 
			theta[i-1]+h*vtheta[i-1]+h*h*k3t/2., t+h);
		k4sy = -k4y*me/ms;
		k4z = funcz(Ex[i-1]+h*vex[i-1]+h*h*k3x/2.,vex[i-1]+h*k3x, Ey[i-1]+h*vey[i-1]+h*h*k3y/2.,vey[i-1]+h*k3y, Ez[i-1]+h*vez[i-1]+h*h*k3z/2.,vez[i-1]+h*k3z,
			Sx[i-1]+h*vsx[i-1]+h*h*k3sx/2., Sy[i-1]+h*vsy[i-1]+h*h*k3sy/2., Sz[i-1]+h*vsz[i-1]+h*h*k3sz/2., phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 
			theta[i-1]+h*vtheta[i-1]+h*h*k3t/2., t+h);
		k4sz = -k4z*me/ms;
		k4t = funct(Ex[i-1]+h*vex[i-1]+h*h*k3x/2., Ey[i-1]+h*vey[i-1]+h*h*k3y/2., Ez[i-1]+h*vez[i-1]+h*h*k3z/2.,Sx[i-1]+h*vsx[i-1]+h*h*k3sx/2., Sy[i-1]+h*vsy[i-1]+h*h*k3sy/2., 			Sz[i-1]+h*vsz[i-1]+h*h*k3sz/2., theta[i-1]+h*vtheta[i-1]+h*h*k3t/2.,vtheta[i-1]+h*k3t, phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., vphi[i-1]+h*k3ph);
		k4ph = funcph(Ex[i-1]+h*vex[i-1]+h*h*k3x/2., Ey[i-1]+h*vey[i-1]+h*h*k3y/2., Ez[i-1]+h*vez[i-1]+h*h*k3z/2.,Sx[i-1]+h*vsx[i-1]+h*h*k3sx/2., 
			Sy[i-1]+h*vsy[i-1]+h*h*k3sy/2., Sz[i-1]+h*vsz[i-1]+h*h*k3sz/2., theta[i-1]+h*vtheta[i-1]+h*h*k3t/2.,vtheta[i-1]+h*k3t, phi[i-1]+h*vphi[i-1]+h*h*k3ph/2., 					vphi[i-1]+h*k3ph);


		//next velocity term
		vex[i] = vex[i-1] + h*(k1x + 2.*k2x + 2.*k3x + k4x)/6.;
		vsx[i] = vsx[i-1] + h*(k1sx + 2.*k2sx + 2.*k3sx + k4sx)/6.;

		vey[i] = vey[i-1] + h*(k1y + 2.*k2y + 2.*k3y + k4y)/6.;
		vsy[i] = vsy[i-1] + h*(k1sy + 2.*k2sy + 2.*k3sy + k4sy)/6.;

		vez[i] = 0;  //vez[i] = vez[i-1] + h*(k1z + 2.*k2z + 2.*k3z + k4z)/6.;
		vsz[i] = 0.0; //vsz[i-1] + h*(k1sz + 2.*k2sz + 2.*k3sz + k4sz)/6.;
		
		//again testing const theta
		///////////vtheta[i] = 0.;
		vtheta[i] = vtheta[i-1] + h*(k1t + 2.*k2t + 2.*k3t + k4t)/6.;
		//vtheta[i] = 0.0;
		vphi[i] = vphi[i-1] + h*(k1ph + 2.*k2ph + 2.*k3ph + k4ph)/6.;

		//next position term
		Ex[i] = Ex[i-1] + h*vex[i-1] + h*h*(k1x+k2x+k3x)/6.;
		Sx[i] = Sx[i-1] + h*vsx[i-1] + h*h*(k1sx+k2sx+k3sx)/6.;	
		//Sx[i] = -Ex[i]*me/ms;

		Ey[i] = Ey[i-1] + h*vey[i-1] + h*h*(k1y+k2y+k3y)/6.;
		Sy[i] = Sy[i-1] + h*vsy[i-1] + h*h*(k1sy+k2sy+k3sy)/6.;	
		//Sy[i] = -Ey[i]*me/ms;

		Ez[i] = 0.0;
		Sz[i] = 0.0;
		//Ez[i] = Ez[i-1] + h*vez[i-1] + h*h*(k1z+k2z+k3z)/6.;
		//Sz[i] = Sz[i-1] + h*vsz[i-1] + h*h*(k1sz+k2sz+k3sz)/6.;	
		//Sz[i] = -Ez[i]*me/ms;

		theta[i] = theta[i-1] + h*vtheta[i-1] + h*h*(k1t+k2t+k3t)/6.;
		//theta[i] = fabs(fmod((theta[i-1] + h*vtheta[i-1] + h*h*(k1t+k2t+k3t)/6.),(2*pi)));
		theta[i] = theta0;
		//test const theta, because theta was blowing up 
		//////theta[i] = 23.45*pi/180.;
		//phi[i] = fabs(fmod((phi[i-1] + h*vphi[i-1] + h*h*(k1ph+k2ph+k3ph)/6.),(2*pi)));
		phi[i] = phi[i-1] + h*vphi[i-1] + h*h*(k1ph+k2ph+k3ph)/6.;
	
		//psi integration 
		double k1, k2, k3, k4;
		k1 = c0 + vphi[i-1]*cos(theta[i-1]); 
       		k2 = c0 + (vphi[i-1]+h*k1ph/2.)*cos(theta[i-1]+h*k1t/2.); 
        	k3 = c0 + (vphi[i-1]+h*k2ph/2.)*cos(theta[i-1]+h*k2t/2.); 
        	k4 = c0 + (vphi[i-1]+h*k3ph)*cos(theta[i-1]+h*k3t); 
        	psi[i]= (h/6.)*(k1+2*k2+2*k3+k4);


		t+=h;
		}	
		return 0;
}


