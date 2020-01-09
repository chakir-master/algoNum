                                                                        /*Reformattage du code */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <locale.h>
///Les procedures  d'affichage
void affiche_matrice(float *t,int valmax, int n);
void affiche_solution(float * t1,int n,int valmax);
///les procedures de calcul

void saisie_matrice_triangulaire(float *t, int n, int valmax);
void saisie_matrice(float *t, int n, int valmax);
void saisie_matrice_dim1(float *t,int n1, int n, int valmax);

void calcul_matrice_unaire(float *t,float *t1, int n, int valmax);
void calcul_matrice_triangulaire(float *t,float *t1, int n, int valmax);

void operation(float *t,float *t1, int n, int valmax);


void Gauss_sans_pivot(float *t,float * t1,int n,int valmax);
void Gauss_pivot_partielle(float *t,float * t1,int n,int valmax);
void Gauss_pivot_total(float *t,float * t1,int n,int valmax);
void Crout(float *t,float * t1,int n,int valmax);
void Cholesky(float *t, float *t1, int n, int valmax);
void jacobie(float *t, float *t1, int n, int valmax);
void Gauss_seidel(float *t, float *t1, int n, int valmax);


void calcul_ligne(float *t,float * t1,int n,int valmax,int l1,int l2);
void permute_ligne(float *t,float * t1,int n,int valmax,int l1,int l2);
void produit_matrice(float *ti,float *t, float *t1, int n, int valmax);
void transpose_matrice(float *ti,float *t,int n, int valmax);

///procedures speciales
void choix(float *t,float *t1,int n, int valmax);
void clear_matrice(float *t,int valmax, int n);

//les fonctions
int saisie_controle(int max);
int validite(float *t,int valmax, int n);
float difference(float *t, int c1, int c2,int n, int valmax);
float distance_absolue(float *t, int c, int n, int valmax);

///Declaration des fonctions
void affiche_matrice(float *t,int valmax, int n)
{
    int i,j;
    for(i=0; i<n; i++)
    {
        printf("|\t");
        for(j=0; j<n; j++)
        {
            printf("%.1f\t",*(t+i*valmax+j));
        }
        printf("|\n\n");
    }
}///Affiche les elements du tableau4
void affiche_solution(float * t1,int n,int valmax)
{
    int i;
    for(i=0; i<n; i++)
    {
        printf("\nX%d=%.2f",i+1,*(t1+i*valmax+0));
    }
}///Affiche la solution sous format de resolution d'equation
void clear_matrice(float *t, int n, int valmax)
{
    int i,j;
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            *(t+i*valmax+j)=0;
        }
    }
}
int saisie_controle(int max)
{
    int val,test;
    do
    {
        printf("Saisir votre choix. (inferieur a %d) :\n",max);
        test=scanf("%d",&val);
    }
    while(test==0|| val<0 || val>max);
    return val;
}/// Saisie controlee de  valeur entier



void saisie_matrice(float *t, int n, int valmax)
{
    int test,i,j;
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            do
            {
                printf("\nA[%d][%d]:",i+1,j+1);
                test=scanf("%f",t+i*valmax+j);
                fflush(stdin);
            }
            while(test==0);
        }
    }
}

void saisie_matrice_dim1(float *t,int n1, int n, int valmax)
{
    int test,i;
    for(i=0; i<n1; i++)
    {
        do
        {
            printf("\nB[%d]:",i+1);
            test=scanf("%f",t+i*valmax+n);
            fflush(stdin);
        }
        while(test==0);
    }

}

void saisie_matrice_triangulaire(float *t, int n, int valmax)
{
    int test,i,j;
    for(i=0; i<n; i++)
    {
        for(j=i; j<n; j++)
        {
            do
            {
                printf("\nA[%d][%d]:",i+1,j+1);
                test=scanf("%f",t+i*valmax+j);
                fflush(stdin);
            }
            while(test==0);
        }
    }
}
/**************************Test sur les matrices***********************************/
int validite(float *t,int valmax, int n)
{
    int i,j,rep=1,sum;
    for(j=0; j<n; j++)
    {
        sum=0;
        for(i=0; i<n; i++)
        {
            if(*(t+i*valmax+j)==0)
            {
                sum++;
            }
            if(sum>=n)
            {
                return 0;
            }
        }
    }
    return rep;
}
/**************************Les calculs sur les matrices***************************/
float distance_absolue(float *t, int c, int n, int valmax)
{
    int i;
    float sum=0;
    for(i=0; i<n; i++)
    {
        sum=sum+pow((*(t+i*valmax+c)),2);
    }
    sum=sqrt(sum);
    return sum;
}/// Valuier absolue c 'est a dire la distance graphique entre deux points
float difference(float *t, int c1, int c2,int n, int valmax)
{
    int i;
    float solution;
    for(i=0; i<n; i++)
    {
        (*(t+i*valmax+c1))=(*(t+i*valmax+c1))-(*(t+i*valmax+c2));
    }
    solution=distance_absolue((float *)t, c1,n,valmax );
    return solution;
}///Difference entre deux coonnes d'une matrice

void permute_ligne(float *t,float * t1,int n,int valmax,int l1,int l2)
{
    float tmp;
    int j;

     for(j=0; j<n; j++)
    {
        tmp=*(t+l1*valmax+j);
        *(t+l1*valmax+j)=*(t+l2*valmax+j);
        printf("\n l = %d", l1);
        *(t+l2*valmax+j)=tmp;

    }
}///Cette fonction permute deux lignes
void calcul_matrice_unaire(float *t,float *t1, int n, int valmax)
{
    int i;
    for(i=0; i<n; i++)
    {
        *(t1+i*valmax+0)=(*(t+i*valmax+n))/(*(t+i*valmax+i));
    }
}
void calcul_matrice_triangulaire(float *t,float *t1, int n, int valmax)
{
    int i,j;
    float s;
    (*(t1+(n-1)*valmax+0))=(*(t+(n-1)*valmax+(n)))/(*(t+(n-1)*valmax+(n-1)));
    for(i=n-2; i>=0; i--)
    {
        s=0;
        for(j=n-1; j>i ; j--)
        {
            s=s+(*(t+i*valmax+j))*(*(t1+(j)*valmax+0));
        }
        (*(t1+i*valmax+0))=((*(t+i*valmax+n))-s)/(*(t+i*valmax+i));
    }
}
void transpose_matrice(float *ti,float *t,int n, int valmax)
{
    int i,j;
    float temp;
    clear_matrice((float *)ti,valmax,n);
    ///printf("\nMatrice ti apres affectation:\n");
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            (*(ti+i*valmax+j))=(*(t+i*valmax+j));
        }
    }
    //affiche_matrice((float *)ti,valmax,n);
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            temp=*(t+i*valmax+j);
            *(ti+i*valmax+j)=*(t+j*valmax+i);
            *(ti+j*valmax+i)=temp;
        }
    }
    ///printf("\nMatrice transopse ti apres affectation:\n");
    affiche_matrice((float *)ti,valmax,n);
}///Transpose de matrice
void produit_matrice(float *ti,float *t, float *t1, int n, int valmax)
{
    int i,j,k;
    clear_matrice((float *)ti,valmax,n);
    for(i=0; i<n; i++)
    {
        for(j=0; j<n; j++)
        {
            *(ti+i*valmax+j)=0;
            for(k=0; k<n; k++)
            {
                *(ti+i*valmax+j)=*(ti+i*valmax+j)+ (*(t+i*valmax+k))*(*(t1+k*valmax+j));
            }
        }
    }
}
/*******************************Les operations directes****************************/

void calcul_ligne(float *t,float * t1,int n,int valmax,int l1,int l2)
{
    int j;
    for(j=l1+1; j<n; j++)
    {
        *(t+(l2)*valmax+j)=(*(t+(l2)*valmax+j))-((*(t+l1*valmax+j))*(*(t+l2*valmax+l1)))/(*(t+l1*valmax+l1));
    }
    *(t+(l2)*valmax+n)=(*(t+(l2)*valmax+n))-((*(t+l1*valmax+n))*(*(t+l2*valmax+l1)))/(*(t+l1*valmax+l1));
    *(t+(l2)*valmax+l1)=0;
}///Il effectue des calculs entre deux lignes


void Gauss_sans_pivot(float *t,float * t1,int n,int valmax)
{
    int i,j;
    for(i=0; i<n; i++)
    {
        for(j=i+1; j<n; j++)
        {
            calcul_ligne((float *)t,(float * )t1,n,valmax,i,j);
        }
    }
    affiche_matrice((float *)t,valmax,n+1);
    calcul_matrice_triangulaire((float *)t,(float *)t1,n,valmax);
    affiche_matrice((float *)t,valmax,n+1);
}

void Gauss_pivot_partielle(float *t,float * t1,int n,int valmax)
{
    int i,j;
    for(i=0; i<n; i++)
    {
        for(j=i+1; j<n; j++)
        {
            if(*(t+i*valmax+0)<*(t+j*valmax+0))
            {
                permute_ligne((float *)t,(float * )t1,n,valmax,i,j);
            }
        }
        for(j=i+1; j<n; j++)
        {
            calcul_ligne((float *)t,(float * )t1,n,valmax,i,j);
        }
    }
    affiche_matrice((float *)t,valmax,n+1);
    calcul_matrice_triangulaire((float *)t,(float *)t1,n,valmax);
    affiche_matrice((float *)t,valmax,n+1);
    system("PAUSE");
}

void Gauss_jordan(float *t,float * t1,int n,int valmax)
{
    int i,j;
    for(i=0; i<n; i++)
    {
        for(j=0; j<i; j++)
        {
            *(t+i*valmax+j)=(*(t+i*valmax+j))/(*(t+i*valmax+i));
        }
        for(j=i+1; j<n+1; j++)
        {
            *(t+i*valmax+j)=(*(t+i*valmax+j))/(*(t+i*valmax+i));
        }
        *(t+i*valmax+i)=1;

        for(j=0; j<i; j++)
        {
            calcul_ligne((float *)t,(float * )t1,n,valmax,i,j);
        }
        for(j=i+1; j<n; j++)
        {
            calcul_ligne((float *)t,(float * )t1,n,valmax,i,j);
        }
    }
    affiche_matrice((float *)t,valmax,n+1);
    calcul_matrice_triangulaire((float *)t,(float *)t1,n,valmax);
}

void Crout(float *t,float * t1,int n,int valmax)
{
    float R[valmax][valmax],L[valmax][valmax],Rt[valmax][valmax],Lt[valmax][valmax];
    float sum;
    float Y[valmax][valmax],B[valmax][valmax];
    int i,j,k;
    clear_matrice((float *)R,valmax,n);
    clear_matrice((float *)L,valmax,n);
    clear_matrice((float *)B,valmax,n);
    for(i=0; i<n; i++)
    {
        R[0][i]=(*(t+0*valmax+i));
        ///printf("\nR[0][%d]=%f",i,*(t+0*valmax+i));
    }
    for(i=1; i<n; i++)
    {
        L[i][0]=(*(t+i*valmax+0))/(*(t+0*valmax+0));
        ///printf("\nL[%d][0]=%f",i,(*(t+i*valmax+0))/(*(t+0*valmax+0)));
    }

    for(i=1; i<n; i++)
    {
        for(j=i; j<n; j++)
        {
            sum=0;
            for(k=0; k<i; k++)
            {
                sum=sum+ L[i][k]*R[k][j];
                ///printf("\n%f=%f+%f*%f",sum,sum, L[i][k],R[k][j]);
            }
            R[i][j]=(*(t+i*valmax+j))-sum;
            ///printf("\nR[%d][%d]:%f=%f-%f",i+1,j+1,R[i][j],*(t+i*valmax+j),sum);
        }
        for(k=i+1; k<n; k++)
        {
            sum=0;
            for(j=0; j<i; j++)
            {
                sum=sum+ L[k][j]*R[j][i];
                ///printf("\n%f=%f+%f*%f",sum,sum, L[k][j],R[j][i]);
            }
            L[k][i]=((*(t+k*valmax+i))-sum)/R[i][i];
            /// printf("\nL[%d][%d]:%f=(%f-%f)/%f",i+1,j+1, L[k][i],*(t+k*valmax+i),sum,R[i][i]);
        }
    }
    for(i=0; i<n; i++)
    {
        L[i][i]=1;
    }
    affiche_matrice((float *)B,valmax,n);
    for(i=0; i<n; i++)
    {
        B[i][0]=*(t+i*valmax+n);
    }
    printf("\nMatrice R\n");
    affiche_matrice((float *)R,valmax,n);
    printf("\nMatrice L\n");
    affiche_matrice((float *)L,valmax,n);
    printf("\nAffiche transpose de L\n");
    transpose_matrice((float *)Lt,(float *)L,n,valmax);
    printf("\nAffiche transpose de R\n");
    transpose_matrice((float *)Rt,(float *)R,n,valmax);
    printf("\nMatrice d'origine\n");
    affiche_matrice((float *)B,valmax,n);
    produit_matrice((float*)Y,(float *)Lt,(float *)B,n,valmax);
    printf("\nMatrice Y=B*L^(-1)\n");
    affiche_matrice((float *)Y,valmax,n);
    produit_matrice((float*)t1,(float *)Rt,(float *)Y,n,valmax);
    printf("\nMatrice X=Y*R^(-1)");
    affiche_matrice((float *)t1,valmax,n);


}

void Cholesky(float *t, float *t1, int n, int valmax)
{
    float L[valmax][valmax],Lt[valmax][valmax];
    float sum;
    float Y[valmax][valmax],B[valmax][valmax];
    int i,j,k;
    clear_matrice((float *)L,valmax,n);
    clear_matrice((float *)B,valmax,n);
    clear_matrice((float *)t,valmax,n);
    /*-------------------------------Saisie de la matrice sous forme symetrique------------------------*/
    printf("La matrice que vous avez saisie ne peut correspondre à la méthode choisit.");
    printf("\nCette fois-ci vous n aurez qu\'a saisir une matrice triangulaire.\n");
    printf("\nRessaisisez la matrice A");
    saisie_matrice_triangulaire((float *)t,n,valmax);
    affiche_matrice((float *)t,valmax,n);
    saisie_matrice_dim1((float *)t,n, n, valmax);
    affiche_matrice((float *)t,valmax,n+1);
    for(i=1; i<n; i++)
    {
        for(j=0; j<i; j++)
        {
            *(t+i*valmax+j)=(*(t+j*valmax+i));
        }
    }
    /*-------------------------------Division de la matrice principale en deux matrices---------------*/
    L[0][0]=sqrt(*(t+0*valmax+0));
    for(i=1; i<n; i++)
    {
        L[i][1]=(*(t+i*valmax+1))/L[0][0];
    }
    for(j=1; j<n; j++)
    {
        for(k=0; k<j; k++)
        {
            sum=sum+ L[j][k]*L[j][k];
        }
        L[j][j]=sqrt((*(t+j*valmax+j))-sum);
    }
    for(i=0; i<n; i++)
    {
        B[i][0]=*(t+i*valmax+n);
    }
    printf("\nMatrice d'origine\n");
    affiche_matrice((float *)B,valmax,n);
    printf("\nMatrice L\n");
    affiche_matrice((float *)L,valmax,n);
    printf("\nAffiche transpose de L\n");
    transpose_matrice((float *)Lt,(float *)L,n,valmax);
    produit_matrice((float*)Y,(float *)Lt,(float *)B,n,valmax);
    printf("\nMatrice Y=B*L^(-1)\n");
    affiche_matrice((float *)Y,valmax,n);
    produit_matrice((float*)t1,(float *)L,(float *)Y,n,valmax);
    printf("\nMatrice X=Y*R^(-1)");
    affiche_matrice((float *)t1,valmax,n);
}

/********************************Les methodes iteratives**************************/
void Gauss_seidel(float *t, float *t1, int n, int valmax)
{
    float sol[valmax][valmax];
    float sum=0;
    int i,j,value=0;
    float X0,X1;
    clear_matrice((float *)sol,valmax,n);
    printf("\nSaisir la matrice X0 solution de de l'equation\n");
    saisie_matrice_dim1((float *)sol,n, 0, valmax);

    do
    {
        for(i=0; i<n; i++)
        {
            sum=(*(t+i*valmax+n));
            for(j=0; j<i; j++)
            {
                sum=sum-((*(t+i*valmax+j))*sol[j][value]);
            }
            for(j=i+1; j<n; j++)
            {
                sum=sum-((*(t+i*valmax+j))*sol[j][value+1]);
            }
            sol[i][value+1]=(sum)/(*(t+i*valmax+i));
        }
        X1=difference((float *)t, value, value+1,n,  valmax);
        X0=distance_absolue((float *)t, value-1,n, valmax);
        ///affiche_matrice((float *)sol,valmax,value);
        printf("\nXn-1:%.2f\t |Xn+1-Xn|:%.2f",X0,X1);
        value++;
    }
    while(X1!=0);
    for(i=0; i<n; i++)
    {
        *(t1+i*valmax+0)=sol[i][value-2];
    }
    ///affiche_matrice((float *)t1,valmax,value);

}

void jacobie(float *t, float *t1, int n, int valmax)
{
    float sol[valmax][valmax];
    float sum=0;
    int i,j,value=0;
    float X0,X1;
    clear_matrice((float *)sol,valmax,n);
    printf("\nSaisir la matrice X0 solution de de l'equation\n");
    saisie_matrice_dim1((float *)sol,n, 0, valmax);

    do
    {
        for(i=0; i<n; i++)
        {
            sum=(*(t+i*valmax+n));
            for(j=0; j<n; j++)
            {
                if(i!=j)
                {
                    sum=sum-((*(t+i*valmax+j))*sol[j][value]);
                }
            }
            sol[i][value+1]=(sum)/(*(t+i*valmax+i));
        }
        X1=difference((float *)t, value, value+1,n,  valmax);
        X0=distance_absolue((float *)t, value-1,n, valmax);
        //affiche_matrice((float *)sol,valmax,value);
        printf("\nXn-1:%.2f\t |Xn+1-Xn|:%.2f",X0,X1);
        value++;
    }
    while(X1!=0);
    for(i=0; i<n; i++)
    {
        *(t1+i*valmax+0)=sol[i][value-2];
    }
    ///affiche_matrice((float *)t1,valmax,value);
}


void operation(float *t,float *t1, int n, int valmax)
{
    int methode;
    system("cls");
    printf("\t\t*       LES METHODES DE RESOLUTION DES SYSTEMES D'EQUATIONS LINEAIRES      *\n");
    printf("\n1-Gauss sans pivot");
    printf("\n2-Gauss pivot partiel");
    printf("\n3-Gauss jordan");
    printf("\n4-Crout");
    printf("\n5-Doolittle");
    printf("\n6-Cholesky");
    printf("\n7-jacobie");
    printf("\n8-Gauss-seidel");
    printf("\nSaisir votre choix :\n\n");
    methode=saisie_controle(8);
    switch (methode)
    {
        case 1:
            Gauss_sans_pivot((float *)t,(float *) t1,n,valmax);
            break;
        case 2:
            Gauss_pivot_partielle((float *)t,(float *) t1,n,valmax);
            break;
        case 3:
            Gauss_jordan((float *)t,(float *) t1,n,valmax);
            break;
        case 4:
            Crout((float *)t,(float *) t1,n,valmax);
            break;
        case 5:
            Crout((float *)t,(float *) t1,n,valmax);
            break;
        case 6:
            Cholesky((float *)t,(float *) t1,n,valmax);
            break;
        case 7:
            jacobie((float *)t,(float *) t1,n,valmax);
            break;
        case 8:
            Gauss_seidel((float *)t,(float *)t1,n,valmax);
            break;
    }
}

/********************************************************************************/
void choix(float *t,float *t1, int n, int valmax)
{
    int rep;
    printf("\nVeuiller saisir la matrice A");
    saisie_matrice((float *)t,n,valmax);
    affiche_matrice((float *)t, valmax, n);
    rep=validite((float *)t, valmax, n);
    if(rep==0)
    {
        printf("La matrice saisit donne des solutions non finies.\nVeuilez recommencer!");
    }
    else
    {
        printf("\nEntrer la matrice B\n");
        saisie_matrice_dim1((float *)t,n, n, valmax);
        affiche_matrice((float *)t, valmax, n+1);
        operation((float *)t,(float *) t1,n,valmax);
        affiche_solution((float*)t1, n, valmax);
    }

}

int main()
{
    system("cls");
    system("color 0A");
    setlocale(LC_CTYPE,"");
    char reponse;
    do
    {
        int valmax=20,n;
        float A[valmax][valmax],B[valmax][valmax];
        printf("\nRESOLUTION DE SYSTEMES D'EQUATIIONS NON LINEAIRES ");
        printf("\nAx = B (A etant une matrice carre)");
        printf("\nSaisir dimension de A :\n");
        n=saisie_controle(valmax);
        clear_matrice((float *)A,valmax,n);
        clear_matrice((float *)B,valmax,2);
        choix((float *)A,(float *)B,n,valmax);
        printf("\nVoulez-vous continuer? \nAppuyer N ou n pour quitter:");
        scanf("%c",&reponse);
        fflush(stdin);
        //system("cls");
    }
    while(toupper(reponse)!='N');
    return 0;
}
