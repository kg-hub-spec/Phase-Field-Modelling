
int pointgeneration(double *conc, fftw_complex **eta, int nwant,int Xmax,int Xmin,int Ymax,int Ymin,int Zmax,int Zmin, int mindst, int Nx, int Ny, int r)
{	int *X, *Y, *Z;
	X=(int*)malloc(sizeof(int)*nwant);
	Y=(int*)malloc(sizeof(int)*nwant);
	Z=(int*)malloc(sizeof(int)*nwant);
	int ntrials=100000000;
	int i,j,k,count=0,flag; 				//here Z is etta//
	int x,y,z;
	//Initialising Zero Values to Arrays//
	X[0]=rand()%(Xmax-Xmin+1)+Xmin;
	Y[0]=rand()%(Ymax-Ymin+1)+Ymin;
	Z[0]=rand()%(Zmax-Zmin+1)+Zmin;
	for (i=0;i<ntrials & count<nwant;i++){
		flag=0;
		x=rand()%(Xmax-Xmin+1)+Xmin; // coordinates generated will be inclusive of min and max values//
		y=rand()%(Ymax-Ymin+1)+Ymin;
		z=rand()%(Zmax-Zmin+1)+Zmin;
		for(j=0;j<=count;j++){
			if (((X[j]-x)*(X[j]-x)+(Y[j]-y)*(Y[j]-y))<mindst*mindst){
				flag=1;
			}
		}
		if (flag==0){
			count=count+1;
			X[count]=x;
			Y[count]=y;
			Z[count]=z;
		}
		
	}
	
	for (i=0;i<count;i++){
		for (j=X[i]-r;j<=X[i]+r;j++){
				for (k=Y[i]-r;k<=Y[i]+r;k++){
					if ((j-X[i])*(j-X[i])+(k-Y[i])*(k-Y[i])<=r*r){
						conc[((k+Ny)%Ny)+Ny*((j+Nx)%Nx)]=1;
						eta[Z[i]][k+Ny*j]=1;
					}
				}
			}
	}
free(X);
free(Y);
free(Z);
return count; 
}