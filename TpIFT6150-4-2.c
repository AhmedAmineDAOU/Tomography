/*------------------------------------------------------*/
/* Prog    : TpIFT6150-4-2                              */
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

float f(float** Mat, int j, int i) {

    if(i < 0 || i >= LENGTH) {i = (i + LENGTH) % LENGTH;}

    if(j < 0 || j >= WIDTH) {j = (j + WIDTH) % WIDTH;}

    return Mat[i][j];
}

/* Rotation d'un angle theta avec l'interpolation plus proche voisin
   param: FR spectre Fourier partie reelle
        : FI spectre Fourier partie imaginaire
        : RR radon partie reelle
        : RI radon partie imaginaire
*/

void plusProcheVoisin(float** source, float** dest, float angle) {

    float xp,yp;
    angle=angle*PI/180;

    for(int i=0; i<LENGTH; i++) {
        for(int j=0; j<WIDTH; j++) {

            xp =  (j - WIDTH/2) * cos(-angle) + (i - LENGTH/2) * sin(-angle) + WIDTH/2;
            yp = -(j - WIDTH/2) * sin(-angle) + (i - LENGTH/2) * cos(-angle) + LENGTH/2;

            xp = floor(xp);
            yp = floor(yp);

            dest[i][j] = f(source, xp, yp);
        }
    }
}
/*----------------------------------------------------*/
/* angle en degre entre BA et BC                      */
/* A(ptar,ptac) - B(ptbr,ptbc) - C(ptcr,ptcc)         */
/*----------------------------------------------------*/
float AngleDeg(int ptar,int ptac,int ptbr,int ptbc,int ptcr,int ptcc)
{
 float num,den;
 float angle;

 //initialisation
 num=den=angle=0.0;

 //calcul potentiel
 num+=(float)(ptac-ptbc)*(ptcc-ptbc);
 num+=(float)(ptar-ptbr)*(ptcr-ptbr);
 den+=(float)sqrt(CARRE(ptac-ptbc)+CARRE(ptar-ptbr));
 den*=(float)sqrt(CARRE(ptcc-ptbc)+CARRE(ptcr-ptbr));

 if (den!=0.0) angle=acos(num/den);
 else angle=0.0;

 if (ptar>ptcr) angle=2*PI-angle;

 return (angle*(180.0/PI));
}


/*Reconstitue le spectre à partir de la transformée de Radon
   param: FR spectre Fourier partie reelle
        : FI spectre Fourier partie imaginaire
        : RR radon partie reelle
        : RI radon partie imaginaire
*/
void reconstituerspectre(float** FR, float** FI, float** RR, float** RI) {
    int  x, y;
    float theta, rho;
    //parcourir moitié haute des matrice
    //suffisante car les matrices de radon et de Fourier sont (reelle->paire) (imaginaire->impaire)
    for(int i=0; i<LENGTH/2; i++)
        for(int j=0; j<WIDTH; j++) {

            x = j - WIDTH/2;
            y = i - LENGTH/2;

            rho = fmod(floor(sqrt(CARRE(x) + CARRE(y)) + WIDTH/2), WIDTH);

            theta = fmod(NB_PROJECTIONS - floor(AngleDeg(x, y, 0, 0, xc, 0)), NB_PROJECTIONS);

            FR[i][j] = RR[(int) theta][(int) rho];
            FI[i][j] = -RI[(int) theta][(int) rho];

            FR[LENGTH - 1 - i][(WIDTH - j) % WIDTH] = RR[(int) theta][(int) rho];
            FI[LENGTH - 1 - i][(WIDTH - j) % WIDTH] = RI[(int) theta][(int) rho];
        }
}


/*------------------------------------------------*/
/* PROGRAMME PRINCIPAL ---------------------------*/
/*------------------------------------------------*/
int main(int argc, char** argv)
{
    int i, j, p;
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

    for(p=0; p < NB_PROJECTIONS ; p++) {
        //rotation bilineaire d'angle p
        plusProcheVoisin(MatriceImgG, Mat1, p );

        //Construction de matrice Radon
        for(i=0; i<LENGTH; i++)
            for(j=0; j<WIDTH; j++) {
                //somme de niveaux de gris
                MatriceRadon[p][j] += Mat1[i][j];
            }
        //construction des VctR a partir de la matrice de Radon
        for(i=0; i<WIDTH; i++) {
            VctR[i] = MatriceRadon[p][i];
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
            MatriceRadonRFFT[p][j] = VctR[j];
            MatriceRadonIFFT[p][j] = VctI[j];
        }
    }


    /*----------------------------------------------------------*/
    /*Sauvegarde de Lenna sous forme d'image pgms */
    SaveImagePgm(NAME_IMG_OUT0, MatriceImgG, LENGTH, WIDTH);
    SaveImagePgm("Radon(c)", MatriceRadonRFFT, LENGTH_RADON, WIDTH_RADON);


    //reconstruction
    reconstituerspectre(MatRFFT, MatIFFT, MatriceRadonRFFT, MatriceRadonIFFT);

    //Remake image partie reelle & imaginaire du spectre
    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    ReMkeImg(MatIFFT, LENGTH, WIDTH);



    IFFTDD(MatRFFT, MatIFFT, LENGTH, WIDTH);
    //mettre l'image au centre
    ReMkeImg(MatRFFT, LENGTH, WIDTH);
    //pour visualisation
    Recal(MatRFFT, LENGTH, WIDTH);
    //RecalMoy(MatRFFT, MatriceImgG, LENGTH, WIDTH);
    SaveImagePgm(NAME_IMG_OUT3, MatRFFT, LENGTH, WIDTH);


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
    //printf("%f",distanceb(64,0));

    return 0;



}


