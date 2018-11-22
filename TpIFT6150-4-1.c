/*------------------------------------------------------*/
/* Prog    : TpIFT6150                                  */
/* Auteur  : Ahmed Amine DAOU & Warshe Wrushabh         */
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
#define xb WIDTH/2
#define yb LENGTH/2
#define xc WIDTH
#define yc LENGTH/2



#define LENGTH_RADON NB_PROJECTIONS
#define WIDTH_RADON  WIDTH

/*------------------------------------------------*/
/* FONCTIONS -------------------------------------*/
/*------------------------------------------------*/






float distancebc(){

       return sqrt(CARRE(xc-xb)+CARRE (yc-yb));
}
float distanceb(int xa,int ya){
   return  sqrt(CARRE(xa-xb)+CARRE(ya-yb));

}

float distancec(int xa,int ya){
   return  sqrt(CARRE(xa-xc)+CARRE(ya-yc));

}

float theta(int xa,int ya)
{
    float numerateur, denumerateur, angle;
    //initialisation
     numerateur=denumerateur=angle=0.0;

    //produit scalaire des vecteur BC et BA
    numerateur=(xa-xb)*(xc-xb) +(ya-yb)*(yc-yb);


    //produit des normes
    denumerateur=distanceb(xa,ya)*distancebc(); //||AB||*||BC||

    if (denumerateur != 0.0) angle=acos(numerateur/denumerateur);
    else angle=0.0;

    if (ya>yc) angle=PI-angle;

    return (angle*(180.0/PI));
}

void reconstituerspectre(float** MatRFFT, float** MatIFFT, float** MatriceRadonRFFT, float** MatriceRadonIFFT) {
    /*int i, j, x, y, angle1, rayon1;

    float angle, rayon, f1, f2;
    int angle,distance;
    for(i=0; i<LENGTH/2; i++)
        for(j=0; j<WIDTH; j++) {
            //partie reelle
            MatRFFT[angle][distance]=

            //partie imaginaire

        }*/
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL ---------------------------*/
/*------------------------------------------------*/
int main(int argc, char** argv)
{
    int i, j, k;
    float **MatriceImgG;
    float **MatriceRadon;
    float **MatriceRadonRFFT;
    float **MatriceRadonIFFT;
    float **MatriceRadonMFFT;
    float **MatRFFT;
    float **MatIFFT;
    float **MatMFFT;
    float **Mat1;
    float **Mat2;
    float *VctR;
    float *VctI;

    /*Allocation memoire des matrices */
    MatriceImgG = fmatrix_allocate_2d(LENGTH, WIDTH);
    MatriceRadon = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatriceRadonRFFT = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatriceRadonIFFT = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatriceRadonMFFT = fmatrix_allocate_2d(LENGTH_RADON, WIDTH_RADON);
    MatRFFT = fmatrix_allocate_2d(LENGTH, WIDTH);
    MatIFFT = fmatrix_allocate_2d(LENGTH, WIDTH);
    MatMFFT = fmatrix_allocate_2d(LENGTH, WIDTH);
    Mat1 = fmatrix_allocate_2d(LENGTH, WIDTH);
    Mat2 = fmatrix_allocate_2d(LENGTH, WIDTH);

    /*Allocation memoire des vecteurs */
    VctR = fmatrix_allocate_1d(WIDTH);
    VctI = fmatrix_allocate_1d(WIDTH);

    /*Initialisation a zero de toutes les matrices */
    for (i = 0; i < LENGTH; i++)
        for (j = 0; j < WIDTH; j++) {
            MatriceImgG[i][j] = 0.0;
            MatRFFT[i][j] = 0.0;
            MatIFFT[i][j] = 0.0;
            MatMFFT[i][j] = 0.0;
            Mat1[i][j] = 0.0;
            Mat2[i][j] = 0.0;
        }

    for (i = 0; i < LENGTH_RADON; i++)
        for (j = 0; j < WIDTH_RADON; j++) {
            MatriceRadon[i][j] = 0.0;
            MatriceRadonRFFT[i][j] = 0.0;
            MatriceRadonIFFT[i][j] = 0.0;
            MatriceRadonMFFT[i][j] = 0.0;
        }

    /*Initialisation a zero de tous les vecteurs */
    for (i = 0; i < WIDTH; i++) {
        VctR[i] = 0.0;
        VctI[i] = 0.0;
    }

    //On charge l'image dans MatriceImgG
    //----------------------------------
    LoadImagePgm(NAME_IMG_IN, MatriceImgG, LENGTH, WIDTH);

    for(k=0; k<NB_PROJECTIONS; k++) {
        //rotation bilineaire d'angle k
        bilineaire(MatriceImgG, Mat1, k );

        //Construction de matrice Radon
        for(i=0; i<LENGTH; i++)
            for(j=0; j<WIDTH; j++) {
                //somme de niveaux de gris
                MatriceRadon[k][j] += Mat1[i][j];
            }
        //construction des VctR a partir de la matrice de Radon
        for(i=0; i<WIDTH; i++) {
            VctR[i] = MatriceRadon[k][i];
            VctI[i] = 0.0;
        }
        //Premier Remake vector de VctR et VctI
        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);
        //transformee de fourier pour chaque ligne(vecteur)
        FFT1D(VctR, VctI, WIDTH);
        //deuxieme remake vector de la partie reelle et imaginaire obtenue par FFT
        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);


       //Maintenant MatriceRadonRFFT & MatriceRadonIFFT constituent la transformée de Radon de Lenna
        for(j=0; j<WIDTH; j++) {
            MatriceRadonRFFT[k][j] = VctR[j];
            MatriceRadonIFFT[k][j] = VctI[j];
        }
    }


    /*----------------------------------------------------------*/
    /*Sauvegarde de Lenna sous forme d'image pgms */
    SaveImagePgm(NAME_IMG_OUT0, MatriceImgG, LENGTH, WIDTH);
    SaveImagePgm("Radon(c)", MatriceRadonRFFT, LENGTH_RADON, WIDTH_RADON);



    //calcul du module de la transformée de radon
    Mod(MatriceRadonMFFT, MatriceRadonRFFT, MatriceRadonIFFT, LENGTH_RADON, WIDTH_RADON);

    SaveImagePgm(NAME_IMG_OUT1, MatriceRadonMFFT, LENGTH_RADON, WIDTH_RADON);


    //reconstruction
    reconstituerspectre(MatRFFT, MatIFFT, MatriceRadonRFFT, MatriceRadonIFFT);


    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    ReMkeImg(MatIFFT, LENGTH, WIDTH);



    IFFTDD(MatRFFT, MatIFFT, LENGTH, WIDTH);
    //mettre l'image au centre
    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    //pour visualisation
    Recal(MatRFFT, LENGTH, WIDTH);
    //RecalMoy(MatRFFT, MatriceImgG, LENGTH, WIDTH);
    SaveImagePgm(NAME_IMG_OUT2, MatRFFT, LENGTH, WIDTH);


    /*Liberation memoire pour les matrices */
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

    /*Liberation memoire pour les vecteurs */
    free(VctR);
    free(VctI);

    /*Commande systeme: visualisation de Imgout.pgm */
    // system("display ImgOut0.pgm&");

    /*retour sans probleme */
    printf("\n C'est fini ... \n\n\n");
    //printf("%f",distanceb(64,0));
    printf("%f",theta(10,120));
    return 0;
}

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
