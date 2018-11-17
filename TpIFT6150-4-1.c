/*------------------------------------------------------*/
/* Prog    :                                           */
/* Auteur  :                                            */
/* Date    :                                            */
/* version :                                            */ 
/* langage : C                                          */
/* labo    : DIRO                                       */
/*------------------------------------------------------*/

/*------------------------------------------------*/
/* FICHIERS INCLUS -------------------------------*/
/*------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "FonctionDemo.h"

/*------------------------------------------------*/
/* DEFINITIONS -----------------------------------*/
/*------------------------------------------------*/
#define NAME_IMG_IN  "lenna"

#define NAME_IMG_OUT0 "ImgOut0"
#define NAME_IMG_OUT1 "ImgOut1" 
#define NAME_IMG_OUT2 "ImgOut2" 
#define NAME_IMG_OUT3 "ImgOut3"  
#define NAME_IMG_OUT4 "ImgOut4"
#define NAME_IMG_OUT5 "ImgOut5" 

#define NB_PROJECTIONS 180
#define LENGTH 128
#define WIDTH  128

#define LENGTH_RADON NB_PROJECTIONS
#define WIDTH_RADON  WIDTH

/*------------------------------------------------*/
/* FONCTIONS -------------------------------------*/                     
/*------------------------------------------------*/


/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL ---------------------------*/                     
/*------------------------------------------------*/
int main()
 {
  int i,j;
  float rotat;
  
  float** MatriceImgG;
  float** MatriceRadon;
  float** MatriceRadonRFFT;
  float** MatriceRadonIFFT;
  float** MatriceRadonMFFT;
  float** MatRFFT;
  float** MatIFFT;
  float** MatMFFT;
  float** Mat1;
  float** Mat2;
  float*  VctR;
  float*  VctI;
  
  /*Allocation memoire des matrices*/
  MatriceImgG=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatriceRadon=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonRFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonIFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatriceRadonMFFT=fmatrix_allocate_2d(LENGTH_RADON,WIDTH_RADON);
  MatRFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatIFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  MatMFFT=fmatrix_allocate_2d(LENGTH,WIDTH);
  Mat1=fmatrix_allocate_2d(LENGTH,WIDTH);
  Mat2=fmatrix_allocate_2d(LENGTH,WIDTH);

  /*Allocation memoire des vecteurs*/
  VctR=fmatrix_allocate_1d(WIDTH);
  VctI=fmatrix_allocate_1d(WIDTH);

  /*Initialisation a zero de toutes les matrices*/
  for(i=0;i<LENGTH;i++) for(j=0;j<WIDTH;j++) 
   { MatriceImgG[i][j]=0.0;
     MatRFFT[i][j]=0.0;
     MatIFFT[i][j]=0.0;
     MatMFFT[i][j]=0.0;
     Mat1[i][j]=0.0;
     Mat2[i][j]=0.0; }

  for(i=0;i<LENGTH_RADON;i++) for(j=0;j<WIDTH_RADON;j++)
     { MatriceRadon[i][j]=0.0;
       MatriceRadonRFFT[i][j]=0.0;
       MatriceRadonIFFT[i][j]=0.0;
       MatriceRadonMFFT[i][j]=0.0; }

  /*Initialisation a zero de tous les vecteurs*/
   for(i=0;i<WIDTH;i++) 
     { VctR[i]=0.0;
       VctI[i]=0.0; }
    
   //On charge l'image dans MatriceImgG
   //----------------------------------
   LoadImagePgm(NAME_IMG_IN,MatriceImgG,LENGTH,WIDTH);
   /* Interpolation bilineaire = deux directions a interpoler*/

  
 
  //-----------------------
  //Nouvelles Fonctions ---
  //-----------------------
  /*----------------------------------------------------------------------*/
  /* Transforme de Fourier monodimensionnelle:                            */
  /* ----------------------------------------                             */
  /* FFT1D(VctR,VctI,WIDTH)                                               */
  /*                                                                      */
  /* VctR: vecteur associe au valeurs reelles                             */   
  /* VctI: vecteur associe au valeurs imaginaires                         */
  /* WIDTH      : Largeur des deux vecteurs                               */
  /* ------                                                               */
  /* Resultat de cette FFT:                                               */
  /* VctR: Partie reelle de la FFT                                        */
  /* VctI: Partie imaginaire de la FFT                                    */
  /*----------------------------------------------------------------------*/    
  /*----------------------------------------------------------------------*/
  /* Transforme de Fourier monodimensionnelle inverse:                    */
  /* ------------------------------------------------                     */
  /* IFFTDD(VctR,VctI,WIDTH)                                              */
  /*                                                                      */ 
  /* VctR: vecteur associe au valeurs reelles                             */   
  /* VctI: vecteur associe au valeurs imaginaires                         */
  /* WIDTH      : Largeur des deux vecteurs                               */
  /* ------                                                               */
  /* Resultat de cette FFT inverse:                                       */
  /* VctR: Partie reelle de la FFT inverse                                */
  /* VctI: Partie imaginaire de la FFT inverse                            */
  /*----------------------------------------------------------------------*/   
  /*----------------------------------------------------------------------*/
  /* ReMkeVct(Vct,WIDTH)                                                  */
  /* --------------------------                                           */
  /* Recadre le Vecteur Vct de largeur WIDTH                              */
  /* selon les 2 cadrants                                                 */
  /*----------------------------------------------------------------------*/  
  /*----------------------------------------------------------------------*/
  /* Module de la FFT1D                                                   */
  /* ------------------                                                   */
  /* ModVct(VctM,VctR,VctI,WIDTH)                                         */ 
  /*                                                                      */
  /* VctR: partie reelle du vecteur                                       */
  /* VctI: partie imaginaire du vecteur                                   */
  /* ------                                                               */
  /* Resultat:                                                            */
  /* VctM: module du vecteur                                              */
  /*----------------------------------------------------------------------*/


  /*-------- FIN ---------------------------------------------*/
  /*----------------------------------------------------------*/
  /*Sauvegarde des matrices sous forme d'image pgms*/
  SaveImagePgm(NAME_IMG_OUT0,MatriceImgG,LENGTH,WIDTH);
  
  /*Liberation memoire pour les matrices*/
  free_fmatrix_2d(MatriceImgG); 
  free_fmatrix_2d(MatriceRadon);
  free_fmatrix_2d(MatriceRadonRFFT);
  free_fmatrix_2d(MatriceRadonIFFT);
  free_fmatrix_2d(MatriceRadonMFFT);
  free_fmatrix_2d(MatRFFT);
  free_fmatrix_2d(MatIFFT);
  free_fmatrix_2d(MatMFFT);
  free_fmatrix_2d(Mat1); 
  free_fmatrix_2d(Mat2); 

  /*Liberation memoire pour les vecteurs*/
  free(VctR);
  free(VctI);   

  /*Commande systeme: visualisation de Ingout.pgm*/
  system("display ImgOut0.pgm&");
  
  /*retour sans probleme*/ 
  printf("\n C'est fini ... \n\n\n");
  return 0; 	 
}
