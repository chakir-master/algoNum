#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <locale.h>

void gaussPivot();

///
void saisieMatrice();
void permute_ligne(int ligne1, int ligne2);
void calcul_ligne(int ligne1, int ligne2);
void affiche_matrice();

///Variables globales
int nLigne, ctrlLigne, ctrlVect, ctrlMat;
double A[20][20], B[20], X[20];

int main()
{
    saisieMatrice();
    gaussPivot();
    return 0;
}

void saisieMatrice()
{
    int cpt, cpt2;
    system("cls");
    do
    {
        printf("\nSaisir le nombre de lignes de la matrice (<=20) : ");
        ctrlLigne = scanf("%d", &nLigne);

        while(ctrlLigne == 0)
        {
            printf("\nVeuillez saisir une valeur réelle ");
            printf("\nRessaisir le nombre de ligne de la matrice : ");
            scanf("%d", &nLigne);
            fflush(stdin);
        }

    }while(nLigne>=20);

    printf("\nVeuillez saisir les éléments de la matrice carrée ligne par ligne : \n");
    for(cpt=1; cpt<=nLigne; cpt++)
    {
        for(cpt2=1; cpt2<=nLigne; cpt2++)
        {
            do
            {
                printf("A[%d][%d] = ", cpt, cpt2);
                ctrlMat = scanf("%lf", &A[cpt][cpt2]);
                fflush(stdin);

                while(ctrlMat == 0)
                {
                    printf("\nLa valeur saisi doit être un nombre réèl...");
                    printf("\nVeuillez ressaisir la valeur : ");
                    ctrlMat = scanf("%lf", &A[cpt][cpt2]);
                    fflush(stdin);
                }
            }while(ctrlMat == 0);
        }
    }

    printf("\nVeuillez saisir les éléments du vecteur B : \n");
    for(cpt=1; cpt<=nLigne; cpt++)
    {
        do
        {
            printf("B[%d] = ", cpt);
            ctrlVect = scanf("%lf", &B[cpt]);
            fflush(stdin);

            while(ctrlVect == 0)
            {
                printf("La valeur saisi doit être un nombre réèl.");
                printf("Veuillez ressaisir la valeur : ");
                ctrlVect = scanf("%lf", &B[cpt]);
            }
        }while(ctrlVect == 0);
    }
}

void permute_ligne(int ligne1, int ligne2)
{
    double tmp;
    int j;
    for(j=0; j<nLigne; j++)
    {
        tmp = A[ligne1][j];
        A[ligne1][j] =A[ligne1][j];
        A[ligne1][j] = tmp;
    }
}

void calcul_ligne(int ligne1, int ligne2)
{
    int j;
    for(j=ligne1+1; j<nLigne; j++)
    {
        A[ligne2][nLigne] = A[ligne2][j] - A[ligne1][j]/A[ligne1][ligne1];
    }
    A[ligne2][nLigne] = A[ligne2][nLigne] - (A[ligne1][nLigne]*A[ligne2][ligne1])/A[ligne1][ligne1];
    A[ligne2][ligne1] = 0;
}

void affiche_matrice()
{
    int i, j;
    for(i=0; i<nLigne; i++)
    {
        printf("|\t");
        for(j=0; j<nLigne; j++)
        {
            printf("%.1f\t", A[i][j]);
        }
        printf("|\n\n");
    }
}

void calcul_matrice_triangulaire()
{
    int i, j;
    double s;
    B[nLigne-1] = A[nLigne-1][nLigne] / A[nLigne-1][nLigne-1];
    for(i=nLigne-2; i>=0; i--)
    {
        s=0;
        for(j=nLigne-1; j>1; j--)
        {
            s = s + A[i][j]*B[j];
        }
        B[i] = (A[i][nLigne]-s) / A[i][i];
    }
}

void gaussPivot()
{
    int i, j;
    for(i=0; i<nLigne; j++)
    {
        for(j=i+1; j<nLigne; j++)
        {
            if(A[i][0]<A[j][0])
            {
                permute_ligne(i, j);
            }
        }
        for(j=i+1; j<nLigne; j++)
        {
            calcul_ligne(i, j);
        }
        ///nLigne += 1;
        affiche_matrice();
        ///calcul_matrice_triangulaire();
        ///nLigne += 1;
        ///affiche_matrice();
    }
}
