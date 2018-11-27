/*---------------------------------------------------*/
/* module  : FonctionDemo.c                          */
/* auteur  : Mignotte Max                            */
/* date    : 19/11/99                                */              
/* langage : C                                       */
/* labo    : DIRO                                    */
/*---------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "FonctionDemo.h"

/*------------------------------------------------*/
/* FONCTIONS -------------------------------------*/
/*------------------------------------------------*/
/*---------------------------------------------------------*/
/*  Alloue de la memoire pour une matrice 1d de float      */
/*---------------------------------------------------------*/
float* fmatrix_allocate_1d(int hsize)
 {
  float* matrix;

  matrix=(float*)malloc(sizeof(float)*hsize); 
  if (matrix==NULL) printf("probleme d'allocation memoire");

  return matrix; 
 }

/*----------------------------------------------------------*/
/*  Alloue de la memoire pour une matrice 2d de float       */
/*----------------------------------------------------------*/
float** fmatrix_allocate_2d(int vsize,int hsize)
 {
  int i;
  float** matrix;
  float *imptr;

  matrix=(float**)malloc(sizeof(float*)*vsize);
  if (matrix==NULL) printf("probleme d'allocation memoire");

  imptr=(float*)malloc(sizeof(float)*hsize*vsize);
  if (imptr==NULL) printf("probleme d'allocation memoire");
 
  for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }

/*----------------------------------------------------------*/
/* Libere la memoire de la matrice 1d de float              */
/*----------------------------------------------------------*/
void free_fmatrix_1d(float* pmat)
 { 
  free(pmat); 
 }

//----------------------------------------------------------*/
/* Libere la memoire de la matrice 2d de float              */
/*----------------------------------------------------------*/
void free_fmatrix_2d(float** pmat)
 { 
  free(pmat[0]);
  free(pmat);
 }

/*----------------------------------------------------------*/
/* Chargement de l'image de nom <name> (en pgm)             */
/*----------------------------------------------------------*/
void LoadImagePgm(char* name,float** mat,int length,int width)
 {
  int i,j,k;
  unsigned char var;
  char buff[NBCHAR];

  char stringTmp1[NBCHAR],stringTmp2[NBCHAR],stringTmp3[NBCHAR];
 
  int ta1,ta2,ta3;
  FILE *fic;

  /*-----nom du fichier pgm-----*/
  strcpy(buff,name);
  strcat(buff,".pgm");
  printf("---> Ouverture de %s",buff);

  /*----ouverture du fichier----*/
  fic=fopen(buff,"r");
  if (fic==NULL)
    { printf("\n- Grave erreur a l'ouverture de %s  -\n",buff);
      exit(-1); }

  /*--recuperation de l'entete--*/
  fgets(stringTmp1,100,fic);
  fgets(stringTmp2,100,fic);
  fscanf(fic,"%d %d",&ta1,&ta2);
  fscanf(fic,"%d\n",&ta3);

  /*--affichage de l'entete--*/
  printf("\n\n--Entete--");
  printf("\n----------");
  printf("\n%s%s%d %d \n%d\n",stringTmp1,stringTmp2,ta1,ta2,ta3);
   
  /*--chargement dans la matrice--*/
     for(i=0;i<length;i++)
      for(j=0;j<width;j++)  
        { fread(&var,1,1,fic);
          mat[i][j]=var; }

   /*---fermeture du fichier---*/
  fclose(fic);
 }


/*----------------------------------------------------------*/
/* Sauvegarde de l'image de nom <name> au format pgm        */
/*----------------------------------------------------------*/
void SaveImagePgm(char* name,float** mat,int length,int width)
 {
  int i,j,k;
  char buff[NBCHAR];
  FILE* fic;
  time_t tm;

  /*--extension--*/
  strcpy(buff,name);
  strcat(buff,".pgm");

  /*--ouverture fichier--*/
  fic=fopen(buff,"w");
    if (fic==NULL) 
        { printf(" Probleme dans la sauvegarde de %s",buff); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format pgm\n",name);

  /*--sauvegarde de l'entete--*/
  fprintf(fic,"P5");
  fprintf(fic,"\n# IMG Module, %s",ctime(&tm));
  fprintf(fic,"%d %d",width,length);
  fprintf(fic,"\n255\n");

  /*--enregistrement--*/
     for(i=0;i<length;i++)
      for(j=0;j<width;j++) 
        fprintf(fic,"%c",(char)mat[i][j]);
   
  /*--fermeture fichier--*/
   fclose(fic); 
 } 

/*-----------------*/
/* FOURIER 2D -----*/
/*-----------------*/
/*------------------------------------------------*/
/*  FOURN ----------------------------------------*/
/*------------------------------------------------*/
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}


/*----------------------------------------------------------*/
/* FFTDD                                                    */
/*----------------------------------------------------------*/
void FFTDD(float** mtxR,float** mtxI,int lgth, int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1)); 

 /*Remplissage de nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*Remplissage de data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,FFT);

 /*Remplissage*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2];  }

 /*Conversion en Matrice*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

    mtxR[posy][posx]=ImgFreqR[i]/(wdth*lgth);  
    mtxI[posy][posx]=ImgFreqI[i]/(wdth*lgth); }

 /*Liberation memoire*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}


/*----------------------------------------------------------*/
/* IFFTDD                                                   */
/*----------------------------------------------------------*/
void IFFTDD(float** mtxR,float**  mtxI,int lgth,int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1));

 /*Remplissage de nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*Remplissage de data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,IFFT);

 /*Remplissage*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2]; }

 /*Conversion en Matrice*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

   mtxR[posy][posx]=ImgFreqR[i];  
   mtxI[posy][posx]=ImgFreqI[i]; }

 /*Liberation memoire*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}

/*-----------------*/
/* FOURIER 1D -----*/
/*-----------------*/
/*------------------------------------------------*/
/*  FOUR1 ----------------------------------------*/
/*------------------------------------------------*/
void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

/*----------------------------------------------------------*/
/* FFT1D                                                    */
/*----------------------------------------------------------*/
void FFT1D(float* vctR,float* vctI,int wdth)
{
 int i,j;

 float* data;
 unsigned long nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*(2*wdth)+1);

 /*Remplissage de nn*/
 nn=wdth;

 /*Remplissage de data*/
 for(i=0;i<wdth;i++) 
   { data[(2*i)+1]=vctR[i];
     data[(2*i)+2]=vctI[i]; }

 /*FFT1D*/
 four1(data,nn,FFT);

 /*Recuperation des donnes*/
 for(i=0;i<wdth;i++)
  { vctR[i]=data[(2*i)+1];
    vctI[i]=data[(2*i)+2];  }

 /*Recalage*/
 for(i=0;i<wdth;i++)
  { vctR[i]/=(wdth);
    vctI[i]/=(wdth);  }

 /*Liberation memoire*/
 free(data);
}

/*----------------------------------------------------------*/
/* IFFT1D                                                    */
/*----------------------------------------------------------*/
void IFFT1D(float* vctR,float* vctI,int wdth)
{
 int i,j;

 float* data;
 unsigned long nn;

 /*allocation memoire*/
 data=(float*)malloc(sizeof(float)*(2*wdth)+1);

 /*Remplissage de nn*/
 nn=wdth;

 /*Remplissage de data*/
 for(i=0;i<wdth;i++) 
   { data[(2*i)+1]=vctR[i];
     data[(2*i)+2]=vctI[i]; }

 /*FFT1D*/
 four1(data,nn,IFFT);

 /*Recuperation des donnes*/
 for(i=0;i<wdth;i++)
  { vctR[i]=data[(2*i)+1];
    vctI[i]=data[(2*i)+2];  }

 /*Liberation memoire*/
 free(data);
}

/*----------------------------------------------------------*/
/* ReMkeVct                                                 */
/*----------------------------------------------------------*/
void ReMkeVct(float* Vct,int wdth)
{
 int j;
 int cj;
 float*  Vct_tmp;

 /*Initialisation*/
 cj=(int)(wdth/2);

 /*Allocation memoire*/
 Vct_tmp=fmatrix_allocate_1d(wdth);

 /*Recadrage*/
 for(j=0;j<cj;j++) Vct_tmp[cj+j]=Vct[j];

 for(j=cj;j<wdth;j++) Vct_tmp[j-cj]=Vct[j];

 /*Transfert*/
 for(j=0;j<wdth;j++) Vct[j]=Vct_tmp[j];

 /*desallocation memoire*/
 free(Vct_tmp);
}

/*----------------------------------------------------------*/
/* Module d'un vecteur                                      */
/*----------------------------------------------------------*/
void ModVct(float* VctM,float* VctR,float* VctI,int wdth)
{
 int j;

 /*Calcul du module*/
 for(j=0;j<wdth;j++)
 VctM[j]=sqrt((VctR[j]*VctR[j])+(VctI[j]*VctI[j]));
}


/*----------------------------------------------------------*/
/* Module de la matrice                                     */
/*----------------------------------------------------------*/
void Mod(float** matM,float** matR,float** matI,int lgth,int wdth)
{
 int i,j;

 /*Calcul du module*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 matM[i][j]=sqrt((matR[i][j]*matR[i][j])+(matI[i][j]*matI[i][j]));
}

/*----------------------------------------------------------*/
/* ReMkeImg                                                 */
/*----------------------------------------------------------*/
void ReMkeImg(float** mat,int lgth,int wdth)
{
 int i,j;
 int ci,cj;
 float** mattmp;

 /*Initialisation*/
 ci=(int)(lgth/2);
 cj=(int)(wdth/2);

 /*Allocation memoire*/
 mattmp=fmatrix_allocate_2d(lgth,wdth);

 /*Recadrage*/
 for(i=0;i<ci;i++) for(j=0;j<cj;j++)
 mattmp[ci+i][cj+j]=mat[i][j];

 for(i=ci;i<lgth;i++) for(j=cj;j<wdth;j++)
 mattmp[i-ci][j-cj]=mat[i][j];

 for(i=0;i<ci;i++) for(j=cj;j<wdth;j++)
 mattmp[ci+i][j-cj]=mat[i][j];

 for(i=ci;i<lgth;i++) for(j=0;j<cj;j++)
 mattmp[i-ci][cj+j]=mat[i][j];

 /*Transfert*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
  mat[i][j]=mattmp[i][j];

 /*desallocation memoire*/
 free_fmatrix_2d(mattmp);
}

/*----------------------------------------------------------*/
/* Inverse le sens d'une matrice                            */
/*----------------------------------------------------------*/
void InverseSense(float** mat,int lgth,int wdth)
{
 int i,j;
 int ci,cj;
 float** mattmp;

 /*Allocation memoire*/
 mattmp=fmatrix_allocate_2d(lgth,wdth);

 /*Transfert*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mattmp[lgth-i-1][wdth-j-1]=mat[i][j];

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][wdth-j-1]=mattmp[i][j];

 /*desallocation memoire*/
 free_fmatrix_2d(mattmp);
}

/*----------------------------------------------------------*/
/* Mult                                                     */
/*----------------------------------------------------------*/
void Mult(float** mat,float coef,int lgth,int wdth)
{
 int i,j;

 /*multiplication*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    { mat[i][j]*=coef;
      if (mat[i][j]>GREY_LEVEL) mat[i][j]=GREY_LEVEL; }
}

/*----------------------------------------------------------*/
/* Recal                                                    */
/*----------------------------------------------------------*/
void Recal(float** mat,int lgth,int wdth)
{
 int i,j;
 float max,min;

 /*Initialisation*/
 max=0.0;
 min=100000000;

 /*Recherche du min*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if (mat[i][j]<min) min=mat[i][j];

 /*plus min*/
   for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][j]-=min;

 /*Recherche du max*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    if (mat[i][j]>max) max=mat[i][j];

 /*Recalibre la matrice*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   mat[i][j]*=(GREY_LEVEL/max);      
}

/*----------------------------------------------------------*/
/* Mult 2 matrices complexes                                */
/*----------------------------------------------------------*/
void MultMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=mat1Rin[i][j]*mat2Rin[i][j]-mat1Iin[i][j]*mat2Iin[i][j];
     matIout[i][j]=mat1Rin[i][j]*mat2Iin[i][j]+mat2Rin[i][j]*mat1Iin[i][j]; }
}

/*----------------------------------------------------------*/
/* Matrice complexe au carre                                */
/*----------------------------------------------------------*/
void SquareMatrix(float** matRout,float** matIout,float** matRin,float** matIin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=CARRE(matRin[i][j])-CARRE(matIin[i][j]);
     matIout[i][j]=2*matRin[i][j]*matIin[i][j]; }
}


/*--------------*/
/* BRUIT -------*/
/*--------------*/
/*************************************************/
/* retourne un nombre aleatoire entre zero et un */
/*************************************************/
float randomize(void)
{ return ((float)rand()/RAND_MAX); }

/*****************************************************/
/* bruit gaussien                                    */
/*****************************************************/
float gaussian_noise(float var,float mean)
{
 float noise,theta;

 /*generation du bruit*/
 noise=sqrt(-2*var*log(1.0-((float)rand()/RAND_MAX)));
 theta=(float)rand()*1.9175345E-4-PI;
 noise=noise*cos(theta);
 noise+=mean;
 if (noise>GREY_LEVEL) noise=GREY_LEVEL;
 if (noise<0) noise=0;
 return noise;
}

/*****************************************************/
/* bruit gaussien sur la matrice                     */
/*****************************************************/
void add_gaussian_noise(float** mat,int lgth,int wdth,float var)
{
 int i,j;

 /*Boucle sur l'image*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 if (var!=0.0) mat[i][j]=gaussian_noise(var,mat[i][j]);
}


/*---------------*/
/* MESURE -------*/
/*---------------*/



/*---------------*/
/* HISTO --------*/
/*---------------*/
/*****************************************************/
/* calcul de l'histogramme                           */
/*****************************************************/
void compute_histo(float** mat,int lgth,int wdth,float* hist)
{
 int i,j;

  for(i=0;i<GREY_LEVEL;i++) hist[i]=0.0;

  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if ((mat[i][j]>=0)&&(mat[i][j]<GREY_LEVEL))
       hist[(int)(mat[i][j])]++;
  
  for(i=0;i<GREY_LEVEL;i++)  hist[i]/=(wdth*lgth);
}

/*******************************************************/
/* sauvegarde de l'histogramme sous forme d'image PGM  */
/*******************************************************/   
void SaveHistoPgm(char* name,float* hist)
{
 int i,j;
 float min,max;
 float** Mat;
 float* histtmp;

 /*Allocation memoire*/
 Mat=fmatrix_allocate_2d(256,256);
 histtmp=fmatrix_allocate_1d(256);
 
 /*Initialisation*/
 for(i=0;i<256;i++) for(j=0;j<256;j++) Mat[i][j]=BLACK; 

 /*Recalage*/
 for(max=0.0,min=1000000.0,i=0;i<GREY_LEVEL;i++)
   { if (hist[i]>max) max=hist[i];
     if (hist[i]<min) min=hist[i]; }

 for(i=0;i<GREY_LEVEL;i++) histtmp[i]=((hist[i]-min)/max)*255.0; 

 /*Affichage*/
 for(i=0;i<GREY_LEVEL;i++) for(j=0;j<histtmp[i];j++) Mat[GREY_LEVEL-j][i]=WHITE;

 /*Sauvegarde*/
 SaveImagePgm(name,Mat,256,256);

 /*Liberation memoire*/
 free_fmatrix_2d(Mat); 
 free_fmatrix_1d(histtmp);  
}

/*---------------*/
/* MARKOV -------*/
/*---------------*/
/******************************************************/
/* retourne la probabilite d'un pixel y d'appartenir  */
/* a la classe correspondante aux parametres          */
/* coef, moy et var d une gaussienne                  */
/******************************************************/
float funcgauss(float pix,float moy,float var)
{
return (float)(1/(sqrt(2*PI*var)))*exp(-1*(CARRE((float)pix-moy))/(2*var));
}

/*---------------*/
/* GRAPH --------*/
/*---------------*/
/****************************************************************/
/* amplitude en x et y d'une gaussienne 2D centree au           */
/*    coord (xmoy,ymoy) et de variance var                      */
/****************************************************************/
float funcgauss2D(int row,int col,int row_moy,int col_moy,float var)
{
 float dist2;
 float amp;

 dist2=CARRE((float)row-(float)row_moy)+CARRE((float)col-(float)col_moy);
 amp=(1.0/(2.0*PI*var))*exp(-1.0*(dist2/(2*var)));

 return amp;
}

/*---------------*/
/* CONVOLUTION --*/
/*---------------*/


/* ---- Convol_U2DBlur ---- */
/* ------------------------ */ 
void Convol_U2DBlur(float** mat,int lgth, int wdth,int size)
 {
  int i,j,k,l,m;
  float** mattmp;
  float tmp;
  float nb;

  /*Allocation memoire*/
  mattmp=fmatrix_allocate_2d(lgth,wdth);

  /*Uniform 2D Blur*/
  for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++) 
    {
     nb=0.0; tmp=0.0;
     for(l=-(int)(size/2);l<=(int)(size/2);l++)
     for(m=-(int)(size/2);m<=(int)(size/2);m++)
       {
	if (((i+l)>0)&&((i+l)<lgth)&&((j+m)>0)&&((j+m)<wdth))
          {  nb++; tmp+=mat[i+l][j+m]; }
       
        else continue;
       }    
     mattmp[i][j]=(tmp/nb);
    }

   /*Recopiage matrice*/
  for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++) 
   mat[i][j]=mattmp[i][j];

  /*Liberation memoire*/
  free_fmatrix_2d(mattmp);
 }

/*----------------------------------------------------------*/
/* Recal                                                    */
/*----------------------------------------------------------*/
void RecalMoy(float** matin,float** matout,int lgth,int wdth)
{
 int i,j;
 float moyin,moyout;

 /*Recherche des deux moyennes*/
 moyin=moyout=0.0;
 for(i=0;i<lgth;i++)  for(j=0;j<wdth;j++)
   { moyin+=matin[i][j];
     moyout+=matout[i][j];}

 /*Recalage*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    matin[i][j]*=(moyout/moyin);

  /*Verification*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   {
    if (matin[i][j]<0.0)   matin[i][j]=0.0;
    if (matin[i][j]>255.0) matin[i][j]=255.0;
   }
}
