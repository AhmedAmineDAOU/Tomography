/*------------------------------------------------------*/
/* Prog    : TpIFT6150-4-1                              */
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

/* Rotation d'un angle theta avec l'interpolation bilineaire
   param: FR spectre Fourier partie reelle
        : FI spectre Fourier partie imaginaire
        : RR radon partie reelle
        : RI radon partie imaginaire
*/

void bilineaire(float** source, float** dest, float angle) {
    int x,y;
    float xp,yp,fxpy,fxpy1;
    angle=angle*PI/180;

    for(int i=0; i<LENGTH; i++) {
        for(int j=0; j<WIDTH; j++) {
            //x'
            xp =  (j - WIDTH/2) * cos(-angle) + (i - LENGTH/2) * sin(-angle) + WIDTH/2;
            //y'
            yp = -(j - WIDTH/2) * sin(-angle) + (i - LENGTH/2) * cos(-angle) + LENGTH/2;

            x = floor(xp);
            y = floor(yp);
            //f(x,y')
            fxpy = f(source, x, y) + (xp - x)*(f(source, x + 1, y) - f(source, x, y));
            //f(x',y+1)
            fxpy1 = f(source, x, y + 1)+ (xp - x) * (f(source, x + 1, y + 1) - f(source, x, y + 1));
            //f(x',y')
            dest[i][j] = fxpy + (yp - y) * (fxpy1 - fxpy);
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
    int  x, y, thetap, rhop;
    float theta, rho, fxpy, fxpy1;
    //parcourir moitié haute des matrice
    //suffisante car les matrices de radon et de Fourier sont (reelle->paire) (imaginaire->impaire)
    for(int i=0; i<LENGTH/2; i++)
        for(int j=0; j<WIDTH; j++) {

            x = j - WIDTH/2;
            y = i - LENGTH/2;

            rho = fmod(floor(sqrt(CARRE(x) + CARRE(y)) + WIDTH/2), WIDTH);

            theta = fmod(NB_PROJECTIONS - floor(AngleDeg(x, y, 0, 0, xc, 0)), NB_PROJECTIONS);

            thetap = floor(theta);
            rhop = floor(rho);

            // Partie reelle du spectre
            fxpy  = RR[thetap][rhop] + (theta - thetap) * (RR[thetap][rhop] - RR[thetap][rhop]);
            fxpy1 = RR[thetap][rhop]+ (theta - thetap) * (RR[thetap][rhop] - RR[thetap][rhop]);

            //paire===> f(-x)= f(x)
            FR[i][j] = fxpy + (rho - rhop) * (fxpy1 - fxpy);
            FR[LENGTH - 1 - i][WIDTH - j] = fxpy + (rho - rhop) * (fxpy1 - fxpy);

            // Partie Imaginaire spectre
            fxpy = RI[thetap][rhop] + (theta - thetap)*(RI[thetap][rhop] - RI[thetap][rhop]);
            fxpy1 = RI[thetap][rhop]+ (theta - thetap) * (RI[thetap][rhop] - RI[thetap][rhop]);
            //Impaire===> f(-x)= - f(x)
            FI[i][j] = -(fxpy + (rho - rhop) * (fxpy1 - fxpy));
            FI[LENGTH - 1 - i][WIDTH - j] = fxpy + (rho - rhop) * (fxpy1 - fxpy);

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
        bilineaire(MatriceImgG, Mat1, p );

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



    //calcul du module de la transformée de radon
    Mod(MatriceRadonMFFT, MatriceRadonRFFT, MatriceRadonIFFT, LENGTH_RADON, WIDTH_RADON);

    SaveImagePgm(NAME_IMG_OUT1, MatriceRadonMFFT, LENGTH_RADON, WIDTH_RADON);


    //reconstruction
    reconstituerspectre(MatRFFT, MatIFFT, MatriceRadonRFFT, MatriceRadonIFFT);

    //Remake vector partie reelle & imaginaire du spectre
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


    /*retour sans probleme */
    printf("\n C'est fini ... \n\n\n");
    //printf("%f",distanceb(64,0));

    return 0;

//NON UTIILSEES

/*distancebc: distance entre B et C*/
float distancebc(){ return sqrt(CARRE(xc-xb)+CARRE (yc-yb)); }

/*distanceb: distance entre A et B
  param: coordonnées du point A
*/
float distanceb(int xa,int ya){ return  sqrt(CARRE(xa-xb)+CARRE(ya-yb)); }

/*distancec: distance entre a c
*/
float distancec(int xa,int ya){ return  sqrt(CARRE(xa-xc)+CARRE(ya-yc)); }

/*rayon : calcule l'angle entre BC et BA
  calcule le produit scalaire BC.BA (numerateur du cos)
  calcule le produit des normes des vecteurs BC ET BA (denumerateur du cos)
  puis calcule l'arctangente de la fraction =>arcos(cos(rayon))
  param: coordonnées du point A
 *return: angle rayon
*/
float theta2(int xa,int ya)
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


    return (angle*(180.0/PI));
}



}

