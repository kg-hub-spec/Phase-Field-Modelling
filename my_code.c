//gcc my_code.c -lm -Iusr/local/include/ -L/usr/local/lib -lfftw3 -lfftw -lm -lfftw3_threads
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<fftw3.h>
#include"pointgeneration.h"
#include<complex.h>
int main()
{
	fftw_init_threads ;
	FILE *ff;
	double dt,K,totaltime,x,delkx,A,delky,delkx1,delky1,*c1,kx,ky,kxpower2,kypower2,kx1,ky1,kx1power2,ky1power2,k1power4,*M;
	double idx,idy,dx,dy;
	int timestep,i,j,k,nx,ny,t,eta_max,num;
	char filename[50] ;
	dt=0.01;		K=6;	totaltime=10;	x=+9;		A=5.5;		dx=1; 	dy=1; 	timestep=500;	 nx=512;	ny=512; eta_max=5;
	fftw_complex *eta[eta_max];
	c1=(double*)malloc(nx*ny*sizeof(double));
	
	double kappa=2,L=1,epsilon=2,value,denominator,ksq;
	char name[100];
	long seed = 786;
	double percentp;
	percentp=0.7;
	int r=12, nwant,count; 					
	nwant=ceil((percentp-0.15)*nx*ny/(0.85*M_PI*r*r));
	
	int Xmin, Xmax, Ymin, Ymax,Zmin,Zmax, mindst=2*r;
	Xmin=r+2;	Xmax=nx-(r+2);	Ymin=r+2;	Ymax=ny-(r+2);		Zmin=0;		Zmax=eta_max-1;
	for(i=0;i<nx;i++){
		for (j=0;j<ny;j++){
		c1[j+ny*i] =0.15;
		}
	}	
	for(num=0;num<eta_max;num++)
	{
		eta[num]=(fftw_complex*) fftw_malloc(nx*ny*sizeof(fftw_complex));
	}
	for (i=0;i<eta_max;i++){
		for(j=0;j<nx;j++){
			for(k=0;k<ny;k++){
			eta[i][k+ny*j]=0;}
		}
	}
	count = pointgeneration(c1,eta,nwant,Xmax,Xmin,Ymax,Ymin,Zmax,Zmin,mindst,nx,ny,r);
	ff=fopen("initial_conc.dat","w");
	for(i=0;i<nx;i++)
	{ for(j=0;j<ny;j++)
		{ fprintf(ff,"%lf\t",c1[j+ny*i]);
		} fprintf(ff,"\n");
		} 
	fclose(ff);
	ff=fopen("initial_etas.dat","w");
	for (k=0;k<eta_max;k++)
	{ for(i=0;i<nx;i++)
		{ for(j=0;j<ny;j++){
		fprintf(ff,"%lf\t",creal(eta[k][j+ny*i]));
		}fprintf(ff,"\n");
		} fprintf(ff,"\n");
		} 
	fclose(ff);
	delkx=2.0*M_PI/((double)nx*dx);
	delky=2.0*M_PI/((double)ny*dy);
	delkx1=2.0*M_PI/((double)nx*dx);
	delky1=2.0*M_PI/((double)ny*dy);
	M=(double*)malloc(nx*ny*sizeof(double));
	fftw_complex *c2, *g_c,*g_eta[eta_max],*gradMmewx,*gradMmewy,*gradmewx,*gradmewy,*mew,*gradMmew;
	c2=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	g_c=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	gradmewx=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	gradmewy=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	gradMmewx=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	gradMmewy=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	gradMmew=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	mew=(fftw_complex*)fftw_malloc(nx*ny*sizeof(fftw_complex));
	for(num=0;num<eta_max;num++){
		g_eta[num]=(fftw_complex*) fftw_malloc(nx*ny*sizeof(fftw_complex));
	}

	fftw_plan plan1,plan2;
	fftw_plan_with_nthreads(8) ;
	plan1=fftw_plan_dft_2d(nx,ny,c2,c2,FFTW_FORWARD,FFTW_ESTIMATE);
	plan2=fftw_plan_dft_2d(nx,ny,c2,c2,FFTW_BACKWARD,FFTW_ESTIMATE);
	int A1=1, A2=2, A3=-3, A4=2, A5=36, A6=4; // Changed this as per paper
	
	double eta_square_sum;
	double var1, var2;
	
	for (i=0;i<nx;i++){
		for(j=0;j<ny;j++){
			c2[j+ny*i]=c1[j+ny*i];
		}
	}
		
	for(t=0;t<=totaltime/dt;t++){
		for (i=0;i<nx;i++){
			for (j=0;j<ny;j++){
				eta_square_sum=0;
				for (num=0;num<eta_max;num++){
					eta_square_sum=eta_square_sum + pow(eta[num][j+ny*i],2);
				}
				g_c[j+ny*i] = 2*A1*c2[j+ny*i]-A2*(eta_square_sum);
				M[j+ny*i]=1+x*creal(c2[j+ny*i]); 
			}
		}
		for (num=0;num<eta_max;num++){
			for (i=0;i<nx;i++){
				for (j=0;j<ny;j++){
					var2=1;
					var1=0;
					for (k=0;k<eta_max;k++){
						if (k!=num){
							var2=var2*pow(eta[k][j+ny*i],2);
							var1=var1+pow(eta[k][j+ny*i],2);
						}
					}
					g_eta[num][j+ny*i] = 2*A2*(1-c2[j+ny*i])*eta[num][j+ny*i] + 4*A3*pow(eta[num][j+ny*i],3) + 6*A4*pow(eta[num][j+ny*i],5) + A5*2*eta[num][j+ny*i]*var1 + A6*2*eta[num][j+ny*i]*var2;
					}
			}
		}
		fftw_execute_dft(plan1,c2,c2);
		fftw_execute_dft(plan1,g_c,g_c);
		for(i=0;i<nx;i++){
			if(i<nx/2){
				kx=(double)(i*delkx);
			}
			else{
				kx=(double)((i-nx)*delkx);
			}
			kxpower2=kx*kx;
			for(j=0;j<ny;j++){		
				if(j<ny/2){
					ky=(double)(j*delky);
				}
				else{
					ky=(double)((j-ny)*delky);
				}
				kypower2=ky*ky;
				mew[j+ny*i]= g_c[j+ny*i] + 2.0*K*(kxpower2+kypower2)*c2[j+ny*i];  // chemical potential mew=g + 2K (kx^2+ky^2) c
				gradmewx[j+ny*i]=-(kx*mew[j+ny*i])*_Complex_I;              
				gradmewy[j+ny*i]=-(ky*mew[j+ny*i])*_Complex_I ;  
			}
		}
		fftw_execute_dft(plan2,gradmewx,gradmewx);
		fftw_execute_dft(plan2,gradmewy,gradmewy);
		fftw_execute_dft(plan2,c2,c2);
		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){   
				c2[j+ny*i]/=(double)(nx*ny);
				gradmewx[j+ny*i]/=(double)(nx*ny);
				gradmewy[j+ny*i]/=(double)(nx*ny);
				gradMmewx[j+ny*i]=M[j+ny*i]*gradmewx[j+ny*i] ;
				gradMmewy[j+ny*i]=M[j+ny*i]*gradmewy[j+ny*i] ;
			}
		}
        fftw_execute_dft(plan1,gradMmewx,gradMmewx);
		fftw_execute_dft(plan1,gradMmewy,gradMmewy);
		fftw_execute_dft(plan1,c2,c2);
		for(i=0;i<nx;i++){
			if(i<nx/2){  
				kx1=(double)(i*delkx1);
			}
			else{
				kx1=(double)((i-nx)*delkx1);
			}
			kx1power2=kx1*kx1;
			for(j=0;j<ny;j++){
				if(j<ny/2){
					ky1=(double)(j*delky1);
				}
				else{
					ky1=(double)((j-ny)*delky1);
				}
				ky1power2=ky1*ky1;
				k1power4=pow(kx1power2+ky1power2,2);
				gradMmew[j+ny*i]=kx1*gradMmewx[j+ny*i]+ky1*gradMmewy[j+ny*i] ;
				c2[j+ny*i]= (c2[j+ny*i]*(1.0+2.0*A*dt*K*k1power4)-dt*gradMmew[j+ny*i]*_Complex_I)/(1.0+2.0*A*dt*K*k1power4);
			} 
		}		
		
		// Evolving Ettas
		for(num=0;num<eta_max;num++)
		{
		fftw_execute_dft(plan1,eta[num],eta[num]);
		fftw_execute_dft(plan1,g_eta[num],g_eta[num]);
		}
		for(i=0;i<nx;i++){
			if(i<nx/2){
				kx=2.0*M_PI*(double)i/((double)nx*dx);
			}
			else{
				kx=2.0*M_PI*(double)(i-nx)/((double)nx*dx);
			}
			for(j=0;j<ny;j++){
				if(j<ny/2){
					ky=2.0*M_PI*(double)j/((double)ny*dy);
				}
				else{
					ky=2.0*M_PI*(double)(j-ny)/((double)ny*dy);
				}
				ksq = kx*kx+ky*ky;
				denominator = 1.0+2.0*L*kappa*dt*ksq;
				for(num=0;num<eta_max;num++){
					eta[num][j+ny*i] = (eta[num][j+ny*i] - L*dt*g_eta[num][j+ny*i])/denominator;
				}
			}	
		}
		
		
		for(num=0;num<eta_max;num++){
			fftw_execute_dft(plan2,eta[num],eta[num]);
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
				eta[num][j+ny*i]/=(double)(nx*ny);
				}
			}
		}
		
		fftw_execute_dft(plan2,c2,c2);
		for(i=0;i<nx;i++){
			for(j=0;j<ny;j++){
				c2[j+ny*i]/=(double)(nx*ny);
			}
		}
		if(t%timestep==0){
			sprintf(filename,"%fconc_out.dat",(int)t*dt) ;
			ff=fopen(filename,"w");
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					fprintf(ff,"%f\t%f\t%lf\n",(float)i*dx,(float)j*dy,(creal(c2[j+ny*i])));
				}
			//fprintf(ff,"\n");	
			}
			fclose(ff);
		}
		//if(t%timestep==0){
			//sprintf(filename,"%deta_output22.dat",t) ;
			//ff=fopen(filename,"w");
			//for(k=0;k<eta_max;k++){
				//for(i=0;i<nx;i++){
					//for(j=0;j<ny;j++)
					//{
						//fprintf(ff,"%lf\t",(creal(eta[k][j+ny*i])));
					//}
				//fprintf(ff,"\n");	
				//}
			//fprintf(ff,"\n");	
			//}
			//fclose(ff);
		//}
	}
	fftw_cleanup_threads();
	fftw_cleanup();
	fftw_free(c2);
	fftw_free(g_c);
	fftw_free(mew);
	fftw_free(gradmewx);
	fftw_free(gradmewy);
	fftw_free(gradMmewx);
	fftw_free(gradMmewy);
	fftw_free(gradMmew);
	free(M);
	fftw_destroy_plan(plan1);
	fftw_destroy_plan(plan2);
	for(num=0;num<eta_max;num++){
		fftw_free(eta[num]);
		fftw_free(g_eta[num]);
	}
	
}
