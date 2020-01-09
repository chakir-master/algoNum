#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <locale.h>

/// PROTOTYPE DES FONCTIONS

///Equations non-lin�iares
//Generaux...
void equation_lineaire();
void systeme_equation_lineaire();
float saisirEntier(char *message);
int saisirEntiern(char *message);
//fonctions
void dichotomie();
void lagrange();
void point_fixe();
void secante();
void newton();
void corde1();
void corde2();
//particuliers
double f(double x);
double df(double x);
float phi(float x);
float derivee_f(float x);

///Système d'équation linéaire
//fonctions des methodes
void Gauss_sans_pivot(float *t,float * t1,int n,int valmax);
void Gauss_pivot_partielle(float *t,float * t1,int n,int valmax);
void Gauss_pivot_total(float *t,float * t1,int n,int valmax);
void Crout(float *t,float * t1,int n,int valmax);
void Cholesky(float *t, float *t1, int n, int valmax);
void jacobie(float *t, float *t1, int n, int valmax);
void Gauss_seidel(float *t, float *t1, int n, int valmax);
//fonctions particulières
    //affichage
void affiche_matrice(float *t,int valmax, int n);
void affiche_solution(float * t1,int n,int valmax);
    //saisie
void saisie_matrice_triangulaire(float *t, int n, int valmax);
void saisie_matrice(float *t, int n, int valmax);
void saisie_matrice_dim1(float *t,int n1, int n, int valmax);
    //calculs
void calcul_matrice_unaire(float *t,float *t1, int n, int valmax);
void calcul_matrice_triangulaire(float *t,float *t1, int n, int valmax);
void operation(float *t,float *t1, int n, int valmax);
void calcul_ligne(float *t,float * t1,int n,int valmax,int l1,int l2);
void permute_ligne(float *t,float * t1,int n,int valmax,int l1,int l2);
void produit_matrice(float *ti,float *t, float *t1, int n, int valmax);
void transpose_matrice(float *ti,float *t,int n, int valmax);
    //controles
int saisie_controle(int max);
int validite(float *t,int valmax, int n);
void choix(float *t,float *t1,int n, int valmax);
void clear_matrice(float *t,int valmax, int n);
float difference(float *t, int c1, int c2,int n, int valmax);
float distance_absolue(float *t, int c, int n, int valmax);



///Definition de fonctions d'equations non-lineaire

float phi(float x)
{
    float ans;
    //ans = 2/(x-1);
    ans = 2*x - 1;
    //ans = sqrt(x+2);
    //ans = -sqrt(x+2;
    //ans =
    return ans;
}
double f(double x) //image de la fonction
{
    double ans;
    //ans = x*x - x -2;
    ans = (x-1)*(x-1);
    return ans;
}

double df(double x) //inage de la derive
{
    double ans;
    //ans = 3*pow(x,2) - 18*x + 26;
    ans = 2*x -2;
    return ans;
}

float derivee_f(float x)
{
    float ans;
    ans = 2*x -2;
    return ans;
}

///Dichotomie
void dichotomie()
{
    char rep;

    do
    {
        int iteration=0,k=0,trouve=0;
        double a=0,b=0,tolerence=0,m=0;
        system("cls");
        printf("\t\t*             METHODE DE Dichotomie           *\n");
        do
        {
            a = saisirEntier("\t\tSaisir a  : ");
            b = saisirEntier("\t\tSaisir b  : ");
            tolerence= saisirEntier("\t\tSaisir la tolerance : 10^-");
            iteration = saisirEntiern("\t\tSaisir le nombre d'iteration  : ");
            tolerence= fabs(tolerence);
            tolerence=1/pow(10,tolerence);
            if (a==b) printf("\n\t\tSaisir a et b tel que (a<b).... Recommencez");
        }
        while(a>b ||a==b);
        printf("\n\t\tIntervalle i = [%lf ; %lf]\n", a, b);
        printf("\t\tf(%lf) = %lf et f(%lf) = %lf\n\n", a, f(a), b, f(b));
        if(fabs(f(b))<= tolerence && fabs(f(a))<= tolerence)
        {
            printf("\n\t\t On a 2 solutions : %.4f et %.4f",a,b);
        }
        else if(fabs(f(a))<= tolerence)
        {
            printf("\n\t\tLa solution est x= %.4f",a);
        }
        else if(fabs(f(b))<= tolerence)
        {
            printf("\n\t\tLa solution est x= %.4f",b);
        }
        else if(f(a)*f(b)>0)
        {
            printf("\n\t\tcette equation admet un nombre paire de solutions");
            exit(0);
        }
        else
        {
            m = (b-a)/2;
            do

            {
                k++;
                if(f(m) == 0)
                {
                    trouve=1;
                    break;
                }
                else if(f(m)*f(b)<0)
                {
                    a = m;
                }
                else
                {
                    b=m;
                }
                m = (b-a)/2;
            }
            while(((fabs(b-a)>=tolerence)||(fabs(f(m))>=tolerence))&&(k<iteration)&&trouve==0);
            if(!trouve)
            {
                printf("\n\t\tLa solution x appartient a l'intervalle [%.4f ,%.4f] ",a,b);
            }
            else
            {
                printf("\n\t\tLa solution est x = %.4f ",m);
            }

        }
        do
        {
            printf("\n\n\t\tvoulez vous recommencer(O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tveuillez saisir soit O (Oui) soit N (Non) : ");
                printf("\t\tvoulez vous refaire une autre partie(O/N) : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }
        while(rep != 'O' && rep != 'N');
    }
    while(rep=='O');
}

///Lagrange
void lagrange()
{
    char rep;

    do
    {
        int iteration=0,trouve=0;
        double a=0,b=0,tolerence=0,m=0,inf=0,sup=0;
        system("cls");
        printf("\n\n\t\t*****************************************************************\n");
        printf("\t\t*             RESOLUTION PAR LA METHODE DE LAGRANGE           *\n");
        printf("\t\t*****************************************************************\n\n");
        do
        {
            a = saisirEntier("\t\tSaisir a  : ");
            b = saisirEntier("\t\tSaisir b  : ");
            tolerence= saisirEntier("\t\tSaisir la tolerance : 10^-");
            iteration = saisirEntiern("\t\tSaisir le nombre d'iteration  : ");
            tolerence= fabs(tolerence);
            tolerence=1/pow(10,tolerence);
            if (a==b) printf("\n\t\ta == b;.... Resaisissez a et b");
        }
        while(a>b ||a==b);
        printf("\n\t\tIntervalle  i = [%lf ; %lf]\n", a, b);
        printf("\t\tf(%lf) = %lf et f(%lf) = %lf\n\n", a, f(a), b, f(b));
        if(fabs(f(b))<= tolerence && fabs(f(a))<= tolerence)
        {
            printf("\n\t\tLes Zeros: %.4f et %.4f",a,b);
        }
        else if(fabs(f(a))<= tolerence)
        {
            printf("\n\t\tLa solution est x= %.4f",a);
        }
        else if(fabs(f(b))<= tolerence)
        {
            printf("\n\t\tLa solution est x= %.4f",b);
        }
        else if(f(a)*f(b)>0)
        {
            printf("\n\t\tL'equation admet un nombre paire de solution");
            exit(0);
        }
        else
        {
            m= a-(b-a)*f(a)/(f(b)-f(a));
            inf = a;
            sup=b;
            while(trouve== 0 && iteration>0 && tolerence< (sup-inf))
            {
                if (f(inf)*f(m)<0)
                {
                    sup=m;

                }
                else if (f(inf)*f(m)>0)
                {
                    inf= m;

                }
                else
                {
                    if (!f(inf))
                        printf("\n\t\tLa solution est x= %f",inf);
                    if(!f(m))
                        printf("\n\t\tLa solution est x= %f",m);
                    trouve=1;
                }
                m=inf-(sup-inf)*f(inf)/(f(sup)-f(inf));
                iteration--;
            }
            if(!trouve)
                printf("\n\t\tLa solution est %f < x < %f a %f pres",inf,sup,tolerence);

        }
        do
        {
            printf("\n\n\t\tVoulez vous-recommencer (O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tveuillez saisir soit O (Oui) soit N (Non) : ");
                printf("\t\tvoulez vous refaire une autre partie(O/N) : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }
        while(rep != 'O' && rep != 'N');
    }
    while(rep=='O');
}

///Point fixe
void point_fixe()
{

    char rep;
    do
    {
        double a=0,tolerence=0,m=0;
        system("cls");
        printf("\t\t*              METHODE DU POINT FIXE            *\n");

        a = saisirEntier("\n\tSaisir x0 : ");
        tolerence= saisirEntier("\n\tEntrer la tolerance : 10^-");
        tolerence= fabs(tolerence);
        tolerence=1/pow(10,tolerence);
        if (f(a)==0)
        {
            printf("\n\t\tLa solution est x = %.4f",a);
        }
        else
        {
            m = a;
            int i =0;
            while(fabs(phi(m)-m)>tolerence)
            {
                i++;
                m = phi(m);
                printf("\n\t\t X%d = %.4f",i,m);
            }

            printf("\n\n\t\tLa solution est x = %.4f",m);
        }
        do
        {
            printf("\n\n\t\tvoulez vous recommencer (O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tveuillez saisir soit O (Oui) soit N (Non) : ");
                printf("\t\tvoulez vous refaire une autre partie(O/N) : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }
        while(rep != 'O' && rep != 'N');
    }
    while(rep=='O');
}

///Secante
void secante()
{

    char rep;
    do
    {
        double a,b,tolerence,m=0;
        system("cls");
        printf("\t\t*              RESOLUTION PAR LA METHODE DE LA SECANTE            *\n");
        do
        {
            a = saisirEntier("\t\tEntrez la valeur de X0 : ");
            b = saisirEntier("\t\tEntrez la valeur de X1 : ");
            tolerence= saisirEntier("\t\tEntrer le crit�re d'arr�t : 10^-");
            tolerence= fabs(tolerence);
            tolerence=1/pow(10,tolerence);
            if (a==b) printf("\n\t\tLes valeures initiales sont egales;.... Resaisissez\n");
        }
        while(a>b ||a==b);
        if(fabs(f(b))<= tolerence && fabs(f(a))<= tolerence)
        {
            printf("\n\t\t On a deux solutions qui sont : %.4f et %.4f",a,b);
        }
        else if(fabs(f(a))<= tolerence)
        {
            printf("\n\t\tLa solution est x= %.4f",a);
        }
        else if(fabs(f(b))<= tolerence)
        {
            printf("\n\t\tLa solution est x= %.4f",b);
        }
        else if(f(a)*f(b)>0)
        {
            printf("\n\t\tcette equation a probablement un nombre paire de solution");
            exit(0);
        }
        else
        {
            do
            {
                m = b-(b-a)*f(b)/(f(b)-f(a));
                a = b;
                b = m;
            }
            while((fabs(b-a)>tolerence)&&(fabs(f(m))>tolerence));
            printf("\n\t\tLa solution est x= %.4f",m);
        }
        do
        {
            printf("\n\n\t\tvoulez vous refaire une autre partie(O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tveuillez saisir soit O (Oui) soit N (Non) : ");
                printf("\t\tvoulez vous refaire une autre partie(O/N) : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }
        while(rep != 'O' && rep != 'N');
    }
    while(rep=='O');
}

///Newton
void newton()
{
    double erreur, x, x_prec;
    int n, iteration = 0;
    const int iter_max = 100;
    system("cls");
    printf("\nEquations non-lineaires : methode de Newton");
    printf("\n----------------------------------------------");


    do
    {
        printf("\nSaisir la tolerance : ");
        printf("\ntolerance = 10^-");
        scanf("%d", &n);
        fflush(stdin);
        erreur = 1 / pow(10, n);



    }while(erreur < 0 || erreur > 0.1);


        printf("\nSaisir x0 : ");
        scanf("%lf", &x_prec);
        fflush(stdin);



    for (int cpt = 0; cpt<iter_max; cpt++)
    {
        iteration++;
        x = x_prec - ( (f(x_prec)) / (df(x_prec)) );

        double test_conv = fabs(x-x_prec)/fabs(x);

        if(test_conv < erreur)
        {
            printf("\nLa convergence est atteinte en X = %lf", x);
            printf("\nLe zero est : %lf", x);
            exit(10);
        }
        x_prec = x;
    }

    if(iteration == iter_max)
    {
        printf("La convergence n'est pas atteinte apr�s %d iterations", iteration);
        exit(10);
    }
}

///Corde1
void corde1()
{
    double a, b, erreur, x, x_prec;
    int n, iteration = 0;
    const int iter_max = 100;
    system("cls");
    printf("\nEquations non-lineaires : methode de la corde1");
    printf("\n----------------------------------------------");
    printf("\nSoit [a, b], l'intervalle de definition de la fonction :\n");
    do
    {
        printf("\nSaisir a : \n");
        scanf("%lf", &a);
        fflush(stdin);
        printf("\nSaisir b : ");
        scanf("%lf", &b);
        fflush(stdin);

        if(a>b)
        {
            printf("\nVeuillez saisir a et b tel que \"a<b\"");
            printf("\nRessaisir a : \n");
            scanf("%lf", &a);
            fflush(stdin);
            printf("\nRessaisir b : ");
            scanf("%lf", &b);
            fflush(stdin);
        }

    }while(a>b);

    do
    {
        printf("\nSaisir la tolerance : ");
        printf("\ntolerance = 10^-");
        scanf("%d", &n);
        fflush(stdin);
        erreur = 1 / pow(10, n);

        printf("\nErreur = %lf", erreur);

    }while(erreur < 0 || erreur > 0.1);

    do
    {
        printf("\nSaisir x0 : ");
        scanf("%lf", &x_prec);
        fflush(stdin);
        if((x_prec < a) || (x_prec > b))
        {
            printf("\nLa valeur initiale doit etre comprise dans l'intervalle. ");
        }

    }while((x_prec < a) || (x_prec > b));

    for (int cpt = 0; cpt<iter_max; cpt++)
    {
        iteration++;
        x = x_prec - ( ((b-a)*f(x_prec)) / (f(b)-f(a)) );

        double test_conv = fabs(x-x_prec)/fabs(x);

        if(test_conv < erreur)
        {
            printf("\nLa convergence est atteinte en X = %lf", x);
            printf("\nf(%lf) = %lf", x, f(x));
            printf("\nAlors le zero est : %lf", x);
            exit(10);
        }
        x_prec = x;
    }

    if(iteration == iter_max)
    {
        printf("La convergence n'est pas atteinte apr�s %d iterations", iteration);
        exit(10);
    }
}


///Corde2

void corde2()
{
    double a, b, erreur, x, x_prec, x_init;
    int n, iteration = 0;
    const int iter_max = 1000;
    system("cls");
    printf("\nEquations non-lineaires : methode de la corde2");
    printf("\n----------------------------------------------");
    printf("\nSoit [a, b], l'intervalle de definition de la fonction :\n");
    do
    {
        printf("\nSaisir a : \n");
        scanf("%lf", &a);
        fflush(stdin);
        printf("\nSaisir b : ");
        scanf("%lf", &b);
        fflush(stdin);

        if(a>b)
        {
            printf("\nVeuillez saisir a et b tel que \"a<b\"");
            printf("\nRessaisir a : \n");
            scanf("%lf", &a);
            fflush(stdin);
            printf("\nRessaisir b : ");
            scanf("%lf", &b);
            fflush(stdin);
        }

    }while(a>b);

    do
    {
        printf("\nSaisir la tolerance : ");
        printf("\ntolerance = 10^-");
        scanf("%d", &n);
        fflush(stdin);
        erreur = 1 / pow(10, n);

    }while(erreur < 0 || erreur > 0.1);

    do
    {
        printf("\nSaisir x0 : ");
        scanf("%lf", &x_init);
        fflush(stdin);
        if((x_init < a) || (x_init > b))
        {
            printf("\nLa valeur initiale doit etre comprise dans l'intervalle [ %lf, %lf] . ", a, b);
        }

    }while((x_init < a) || (x_init > b));

    for (int cpt = 0; cpt<iter_max; cpt++)
    {
        iteration++;
        x = x_prec - ( (f(x_prec)) / (df(x_init)) );

        double test_conv = fabs(x-x_prec)/fabs(x);

        if(test_conv < erreur)
        {
            printf("\nLa convergence est atteinte en X = %lf", x);
            printf("\nf(%lf) = %lf", x, f(x));
            printf("\nAlors le zero est : %lf\n", x);
            exit(10);
        }
        x_prec = x;
    }

    if(iteration == iter_max)
    {
        printf("La convergence n'est pas atteinte apr�s %d iterations", iteration);
        exit(10);
    }
}

///FIN E-N-L
///****************************************************************************************************************
///****************************************************************************************************************

///Definition de fonctions de systeme d'equations non-lineaire
//Gauss sans pivot
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

//Gauss avec pivot
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

//Gauss Jordan
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

//Crout
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
    }
    for(i=1; i<n; i++)
    {
        L[i][0]=(*(t+i*valmax+0))/(*(t+0*valmax+0));
    }

    for(i=1; i<n; i++)
    {
        for(j=i; j<n; j++)
        {
            sum=0;
            for(k=0; k<i; k++)
            {
                sum=sum+ L[i][k]*R[k][j];
            }
            R[i][j]=(*(t+i*valmax+j))-sum;
        }
        for(k=i+1; k<n; k++)
        {
            sum=0;
            for(j=0; j<i; j++)
            {
                sum=sum+ L[k][j]*R[j][i];
            }
            L[k][i]=((*(t+k*valmax+i))-sum)/R[i][i];

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
    printf("\nTranspose de L\n");
    transpose_matrice((float *)Lt,(float *)L,n,valmax);
    printf("\nTranspose de R\n");
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

//Doolittle
void doolittle()
{
    printf("En cours...");
}

//Cholesky
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

//Jacobi
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
        printf("\nXn-1:%.2f\t |Xn+1-Xn|:%.2f",X0,X1);
        value++;
    }
    while(X1!=0);
    for(i=0; i<n; i++)
    {
        *(t1+i*valmax+0)=sol[i][value-2];
    }
}

//Gauss Seidel
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
        printf("\nXn-1:%.2f\t |Xn+1-Xn|:%.2f",X0,X1);
        value++;
    }
    while(X1!=0);
    for(i=0; i<n; i++)
    {
        *(t1+i*valmax+0)=sol[i][value-2];
    }
}

//defintions des fonctions particulieres
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
}

void affiche_solution(float * t1,int n,int valmax)
{
    int i;
    for(i=0; i<n; i++)
    {
        printf("\nX%d=%.2f",i+1,*(t1+i*valmax+0));
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

void calcul_ligne(float *t,float * t1,int n,int valmax,int l1,int l2)
{
    int j;
    for(j=l1+1; j<n; j++)
    {
        *(t+(l2)*valmax+j)=(*(t+(l2)*valmax+j))-((*(t+l1*valmax+j))*(*(t+l2*valmax+l1)))/(*(t+l1*valmax+l1));
    }
    *(t+(l2)*valmax+n)=(*(t+(l2)*valmax+n))-((*(t+l1*valmax+n))*(*(t+l2*valmax+l1)))/(*(t+l1*valmax+l1));
    *(t+(l2)*valmax+l1)=0;
}

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
}

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
}

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
}

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
}

///FIN S-E-L


///**********main() - equation non lineaire
void equation_lineaire()
{

    int choix_met;
    char rep;
    do{
             system("cls");

    printf("\t\t*       LES METHODES DE RESOLUTION DES EQUATIONS NON LINEAIRES      *\n");

    printf("\n\t\t\t1- Dichotomie");
    printf("\n\t\t\t2- Lagrange");
    printf("\n\t\t\t3- Point fixe");
    printf("\n\t\t\t4- Secante");
    printf("\n\t\t\t5- Newton");
    printf("\n\t\t\t6- Corde 1");
    printf("\n\t\t\t7- Corde 2");

    do
    {
        printf("\n\n\t\tVotre choix : ");
        scanf("%d", &choix_met);
        fflush(stdin);
    }
    while( choix_met < 1 || choix_met > 7);

    switch(choix_met)
    {
    case 1 :
        dichotomie();
        break;
    case 2 :
        lagrange();
        break;
    case 3 :
        point_fixe();
        break;
    case 4 :
        secante();
        break;
    case 5 :
        newton();
        break;
    case 6 :
        corde1();
        break;
    case 7 :
        corde2();
        break;
    }
    do
        {
            printf("\n\n\t\tRevenir au Menu des methodes non lineaire(O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tveuillez saisir soit O (Oui) soit N (Non) : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }while(rep != 'O' && rep != 'N');
    }
    while(rep=='O');

}


///**********main() - systeme d'equa lineaire
void systeme_equation_lineaire()
{
    char reponse;
    char rep;
    do
    {

        system("cls");
        setlocale(LC_CTYPE,"");
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
        }while(toupper(reponse)!='N');

        //*****************************************************************************//
        do
        {
            printf("\n\n\t\tRevenir au Menu des methodes non lineaire(O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tveuillez saisir soit O (Oui) soit N (Non) : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }while(rep != 'O' && rep != 'N');
    }while(rep=='O');
}


float saisirEntier(char *message)
{
    float n;
    int retour;
    // controle de la saisie
    do
    {
        printf("\n %s",message);
        retour= scanf(" %f",&n);
        fflush(stdin);
    }
    while(!retour);
    return n;
}

int saisirEntiern(char *message)
{
    int n;
    int retour;
    // controle de la saisie
    do
    {
        printf("\n %s",message);
        retour= scanf(" %d",&n);
        fflush(stdin);
    }
    while(!retour);
    return n;
}





int main()
{
    setlocale(LC_CTYPE,"");
    //Variables
    char choix_ini;
    char rep;
    do{
             system("cls");

  ;
    printf("\n\t\t*                LES METHODES NUMERIQUES                  *\n");

    printf("\n\t\t\tA- Equation non lineaire ");
    printf("\n\t\t\tB- Systeme d'equation lineaire");

    do
    {
        printf("\n\n\t\tVotre choix : ");
        scanf("%c", &choix_ini);
        fflush(stdin);
        choix_ini = toupper(choix_ini);

    }
    while(choix_ini != 'A' && choix_ini != 'B' );

    switch(choix_ini)
    {
    case 'A':
        equation_lineaire();
        break;
    case 'B':
        systeme_equation_lineaire();
        break;
    }
    do
        {
            printf("\n\n\t\tRevenir aux Methodes Numeriques(O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tveuillez saisir soit O (Oui) soit N (Non) : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }while(rep != 'O' && rep != 'N');
    }
    while(rep=='O');
    return 0;
}

