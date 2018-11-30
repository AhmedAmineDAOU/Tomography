/*------------------------------------------------------*/
/* Prog    : TpIFT6150-4-3                              */
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
#define LENGTH 128
#define WIDTH  128

#define xb WIDTH/2
#define yb LENGTH/2
#define xc WIDTH
#define yc LENGTH/2
#define LENGTH_RADON NB_PROJECTIONS
#define WIDTH_RADON  WIDTH




#define NB_PROJECTIONS 90

/*------------------------------------------------*/
/* FONCTIONS -------------------------------------*/
/*------------------------------------------------*/

float f(float** Mat, int j, int i) {

    if(i < 0 || i >= LENGTH) {i = (i + LENGTH) % LENGTH;}

    if(j < 0 || j >= WIDTH){ j = (j + WIDTH) % WIDTH; }

    return Mat[i][j];
}



void bilineaire(float** source, float** dest, float theta) {
    int x1,y1,j,i;
    float xp,yp,fxpy,fxpy1;

    for(int i=0; i<LENGTH; i++) {
        for(int j=0; j<WIDTH; j++) {

            xp =  (j - WIDTH/2) * cos(-theta) + (i - LENGTH/2) * sin(-theta) + WIDTH/2;
            yp = -(j - WIDTH/2) * sin(-theta) + (i - LENGTH/2) * cos(-theta) + LENGTH/2;

            x1 = floor(xp);
            y1 = floor(yp);

            fxpy = f(source, x1, y1) + (xp - x1)*(f(source, x1 + 1, y1) - f(source, x1, y1));
            fxpy1 = f(source, x1, y1 + 1)+ (xp - x1) * (f(source, x1 + 1, y1 + 1) - f(source, x1, y1 + 1));

            dest[i][j] = fxpy + (yp - y1) * (fxpy1 - fxpy);
        }
    }
}


/*----------------------------------------------------*/
/* theta en degre entre BA et BC                      */
/* A(ptar,ptac) - B(ptbr,ptbc) - C(ptcr,ptcc)         */
/*----------------------------------------------------*/
float AngleDeg(int ptar,int ptac,int ptbr,int ptbc,int ptcr,int ptcc)
{
    float num,den;
    float theta;

    //initialisation
    num=den=theta=0.0;

    //calcul potentiel
    num+=(float)(ptac-ptbc)*(ptcc-ptbc);
    num+=(float)(ptar-ptbr)*(ptcr-ptbr);
    den+=(float)sqrt(CARRE(ptac-ptbc)+CARRE(ptar-ptbr));
    den*=(float)sqrt(CARRE(ptcc-ptbc)+CARRE(ptcr-ptbr));

    if (den!=0.0) theta=acos(num/den);
    else theta=0.0;

    if (ptar>ptcr) theta=2*PI-theta;

    return (theta*(180.0/PI));
}
/*Reconstitue le spectre à partir de la transformée de Radon
   param: FR spectre Fourier partie reelle
        : FI spectre Fourier partie imaginaire
        : RR radon partie reelle
        : RI radon partie imaginaire
*/
void reconstituerspectre(float** FR, float** FI, float** RR, float** RI) {
    int  x, y, thetap, rhop;

    float theta, rho, fxpy, fxpy1;
    //parcourir moitié haute des matrice
    //suffisante car les matrices de radon et de Fourier sont (reelle->paire) (imaginaire->impaire)

    for(int i=0; i<LENGTH/2; i++)
        for(int j=0; j<WIDTH; j++) {

            x = j - WIDTH/2;
            y = i - WIDTH/2;

            rho = fmod(round(sqrt(CARRE(x) + CARRE(y)) + WIDTH/2), WIDTH);

            theta = fmod(180 - round(AngleDeg(x, y, 0, 0, xc, 0)), 180) / (180.0 / NB_PROJECTIONS);


            thetap = floor(theta);
            rhop = floor(rho);

            fxpy = RR[thetap][rhop] + (theta - thetap)*(RR[(thetap + 1) % NB_PROJECTIONS][rhop] - RR[thetap][rhop]);
            fxpy1 = RR[thetap][(rhop + 1) % WIDTH_RADON]
                + (theta - thetap) * (RR[(thetap + 1) % NB_PROJECTIONS][(rhop + 1) % WIDTH_RADON] - RR[thetap][(rhop + 1) % WIDTH_RADON]);

            //paire===> f(-x)= f(x)

            FR[i][j] = fxpy + (rho - rhop) * (fxpy1 - fxpy);
            FR[LENGTH - 1 - i][(WIDTH - j) % WIDTH] = fxpy + (rho - rhop) * (fxpy1 - fxpy);

            // Partie Imaginaire spectre
            fxpy = RI[thetap][rhop] + (theta - thetap)*(RI[(thetap + 1) % NB_PROJECTIONS][rhop] - RI[thetap][rhop]);
            fxpy1 = RI[thetap][(rhop + 1) % WIDTH_RADON]
                + (theta - thetap) * (RI[(thetap + 1) % NB_PROJECTIONS][(rhop + 1) % WIDTH_RADON] - RI[thetap][(rhop + 1) % WIDTH_RADON]);
            //Impaire===> f(-x)= - f(x)
            FI[i][j] = -(fxpy + (rho - rhop) * (fxpy1 - fxpy));
            FI[LENGTH - 1 - i][(WIDTH - j) % WIDTH] = fxpy + (rho - rhop) * (fxpy1 - fxpy);
        }
}

/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL ---------------------------*/
/*------------------------------------------------*/
int main(int argc, char** argv)
{
    int i, j;
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


    LoadImagePgm(NAME_IMG_IN, MatriceImgG, LENGTH, WIDTH);

    for(int p=0; p < NB_PROJECTIONS ; p++) {

        bilineaire(MatriceImgG, Mat1, (p * 180.0 / NB_PROJECTIONS) / 360.0 * 2 * PI);

        for(i=0; i<LENGTH; i++)
            for(j=0; j<WIDTH; j++) {
                MatriceRadon[p][j] += Mat1[i][j];
            }

        for(i=0; i<WIDTH; i++) {
            VctR[i] = MatriceRadon[p][i];
            VctI[i] = 0.0;
        }

        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);

        FFT1D(VctR, VctI, WIDTH);
        ReMkeVct(VctR, WIDTH);
        ReMkeVct(VctI, WIDTH);

        for(j=0; j<WIDTH; j++) {
            MatriceRadonRFFT[p][j] = VctR[j];
            MatriceRadonIFFT[p][j] = VctI[j];
        }
    }

    /*----------------------------------------------------------*/
    /*Sauvegarde des matrices sous forme d'image pgms */

    //module de la transformee de radon
    Mod(MatriceRadonMFFT, MatriceRadonRFFT, MatriceRadonIFFT, LENGTH_RADON, WIDTH_RADON);
    
    //reconstruction
    reconstituerspectre(MatRFFT, MatIFFT, MatriceRadonRFFT, MatriceRadonIFFT);

    Mod(MatMFFT, MatRFFT, MatIFFT, LENGTH, WIDTH);
    //Remake image partie reelle & imaginaire du spectre
    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    ReMkeImg(MatIFFT, LENGTH, WIDTH);
    //transformee de fourier inverse
    IFFTDD(MatRFFT, MatIFFT, LENGTH, WIDTH);
    //remake image de la reconstruction
    ReMkeImg(MatRFFT, LENGTH, WIDTH);

    //pour visualisation
    Recal(MatRFFT, LENGTH, WIDTH);

    //changer le nom de l'image suivant l'angle choisi
    if (NB_PROJECTIONS == 90) SaveImagePgm(NAME_IMG_OUT4, MatRFFT, LENGTH, WIDTH);
    else if (NB_PROJECTIONS ==45 ) SaveImagePgm(NAME_IMG_OUT5, MatRFFT, LENGTH, WIDTH);


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

    
    /*retour sans probleme */
    printf("\n C'est fini ... \n\n\n");
    return 0;
}

