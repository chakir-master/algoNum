#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <locale.h>

/// PROTOTYPE DES FONCTIONS

///Equations non-lineiares
void equation_lineaire();
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

///Systeme d'equation lineaire
void systeme_equation_lineaire();
//fonctions principales
void gauss(float a[19][19],float b[19],int n);
void gaussPivot(float a[19][19],float b[19],int n);
void gaussJordan(float a[19][19],float b[19],int n);
void crout(float a[19][19],float b[19],int n);
void doolittle(float a[19][19],float b[19],int n);
void cholesky(float a[19][19],float b[19],int n);
void jacobie(float a[19][19],float b[19],int n);
void gaussSeidel(float a[19][19],float b[19],int n);
//fonctions particulieres
//Saisies
float saisirEntier(char *message);
int saisirEntiern(char *message);
///operations
float norme(float x[19],int n);
void saisirMatrice(float A[19][19],float B[19],int n);
void afficheSysteme(float A[19][19],float B[19],int n);
void afficheMatrice(float A[19][19],int n);
void rechercheZero(float A[19][19],float B[19],int n);
void coMatrices(float A[19][19],float c[19][19],int i,int j,int n);
float determinant(float A[19][19],int n);
///Variables globales
//systeme d equations lineaires
float A[19][19],B[19];
int n,i,j;
char valider;

///Interpolation lineaire
void interpolation();
//fonctions principales
void Newton(int n, float* x, float* y);
void Lagrange(int n,float* x, float* y);
void Moindre(int n, float* x, float* y);
///fonctions particulieres
float* AllocatVect(int n);
float** AllocatMat(int ligne, int colonne);
void SaisieDeDonnees(int n,float* x, float* y);
float* choleskyInter(float** Mat,float* b,int n);
int ctrsaisieint(void);
float ctrsaisiefloat(void);
float* CalculCoefDePhi(int n,float* x, float* y);
float* CalculCoefDep(int n, float* x, float* y,int p);
float** CalculDiffDiv(int n, float* x, float* y);
float* CalculCoefDeN(int n, float* x, float* y);

///Equations différentielles
void equaDif();
//fonctions principales
void TraitementEuler(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h);
void TraitementKunta(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h);
void TraitementEulerModifie(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h);
void TraitementKuntaOrdre4(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h);
//fonctions particulieres
//suites
float SuiteDeleur(float h, float xn, float yn);
float SuiteDeKunta(float h, float xn, float yn);
///operations
int Menu(void);
float* AllocatVect(int n);
float FonctionDerivee(float x, float y);
float ctrsaisiefloat(void);
int ctrsaisieint(void);

///*********************************************************************
///************************************************                     *********
///Definition de fonctions d'equations non-lineaire                     ************************
///************************************************                     *********
///*********************************************************************
float phi(float x)
{
    float ans;
    ///ans = -(x*x) +1;
    ///ans = sqrt(-x+1);
    ///ans = -(sqrt(-x+1));
    ans = 1/(x+1);
    return ans;
}
double f(double x) //image de la fonction
{
    double ans;
    ans = x*x + x -1;
    return ans;
}

double df(double x) //image de la derive
{
    double ans;
    ans = 2*x + 1;
    return ans;
}

float derivee_f(float x)
{
    float ans;
    ans = 2*x + 1;
    return ans;
}
//*************************
//Methodes
//*************************

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
        while(a>b || a==b);

        printf("\n\t\tIntervalle i = [%lf ; %lf]\n", a, b);
        printf("\t\tf(%lf) = %lf et f(%lf) = %lf\n\n", a, f(a), b, f(b));
        if(fabs(f(b))<= tolerence && fabs(f(a))<= tolerence)
        {
            printf("\n\t\t On a 2 solutions : %.4f et %.4f",a,b);
        }
        else if(fabs(f(a))<= tolerence)
        {
            printf("\n\t\tLa solution est x = %.4f",a);
        }
        else if(fabs(f(b))<= tolerence)
        {
            printf("\n\t\tLa solution est x = %.4f",b);
        }
        else if(f(a)*f(b)>0)
        {
            printf("\n\t\tCette equation admet un nombre paire de solutions");
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
                if(i == 1000){
                    printf("\nLa methode ne connverge pas apres %d iterations", i);
                    exit(10);                }
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


///*********************************************************************
///************************************************                     *********
///Definition de fonctions de systeme d equations non-lineaire                   *********
///************************************************                     *********
///*********************************************************************

///Fonctions particulieres
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

float norme(float x[19],int n)
{
    float ref;
    int i;
    ref=0;
    for(i=0; i<n; i++) if (x[i]>ref) ref=x[i];
    return(ref);
}



void saisirMatrice(float a[19][19],float b[19],int n)
{
    int i,j,retour =0;
    printf("\n\t\tSaisie de la matrice A:\n\n");
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            do
            {
                printf("\t\tA[%d][%d] = ",i+1,j+1);
                retour = scanf("%f",&a[i][j]);
                if(retour == 0)printf("La valeur doit etre reelle.");
            }
            while(retour == 0);

        }
    }
    retour=0;
    printf("\n\t\tSaisir le vecteur B:\n\n");
    for(i=0; i<n; i++)
    {
        do
        {
            printf("\t\tB[%d] = ",i+1);
            retour = scanf("%f",&b[i]);
            if(retour == 0)printf("La valeur doit etre reelle.");
        }
        while(retour == 0);
    }
}


void afficheSysteme(float a[19][19],float b[19],int n)
{
    int i,j;
    printf("\n\n");
    for (i=0; i<n; i++)
    {
        printf("\t\t |");
        for (j=0; j<n; j++)
        {
            printf("%.2f ",a[i][j]);
        }
        printf("| | %.2f |",b[i]);
        printf("\n");
    }
}

// fonction d'affichage matrice

void afficheMatrice(float a[19][19],int n)
{
    int i,j;
    printf("\n\n");
    for (i=0; i<n; i++)
    {
        printf("\t\t |");
        for (j=0; j<n; j++)
        {
            printf("%.2f ",a[i][j]);
        }
        printf("|\n");
    }
}

// Mettre � Zero les elements qui doivent etre des z�ro

void rechercheZero(float a[19][19],float b[19],int n)
{
    int i,j;
    float eps=1e-4;
    for(i=0; i<n; i++)
    {
        for (j=0; j<n; j++) if (fabs(a[i][j])<eps) a[i][j]=0;
        if (fabs(b[i])<eps) b[i]=0;
    }
}

void coMatrices(float a[19][19],float c[19][19],int i,int j,int n)
{
    int l,k;
    for(l=0; l<n; l++) for(k=0; k<n; k++)
        {
            if ((l<i)&&(k<j)) c[l][k]=a[l][k];
            if ((l>i)&&(k<j)) c[l-1][k]=a[l][k];
            if ((l<i)&&(k>j)) c[l][k-1]=a[l][k];
            if ((l>i)&&(k>j)) c[l-1][k-1]=a[l][k];
        }
}
// calcul du determinant

float determinant(float a[19][19],int n)
{
    int k,j;
    float c[19][19],s;

    k=n-1;

    if(n==0) return(1);

    s=0;
    for(j=0; j<n; j++)
    {
        coMatrices(a,c,k,j,n);
        s=s+pow(-1,k+j)*a[k][j]*determinant(c,k);
    }
    return(s);
}


//****************
///Methodes
//****************
//Gauss sans pivot
void gauss(float A[19][19],float B[19],int n)
{
    float a[19][19],b[19];
    int cpt,cpt1;
    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }
    system("cls");
    printf("\n\t\t              --------------------------------------           *\n");
    printf("\t\t              --------------------------------------           *\n");
    printf("\t\t              | RESOLUTION PAR LA METHODE DE GAUSS |           *\n");
    printf("\t\t             --------------------------------------           *\n");
    printf("\t\t              --------------------------------------           *\n");

    float x[19],p,s;
    int i,j,k;

    for(k=0; k<n; k++)
    {
        if (a[k][k]==0)
        {
            printf("\n\n\tUn des pivots est nul alors la methode de Gauss simple est non applicable, utilisez une autre methode...\n\n");
            exit(10);
        }

        //reduction
        for(i=k+1; i<n; i++)
        {
            p=a[i][k]/a[k][k];
            for (j=k; j<n; j++) a[i][j]=a[i][j]-p*a[k][j];
            b[i]=b[i]-p*b[k];
        }
    }

    //Resolution
    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++)s=s+a[i][j]*x[j];
        x[i]=(b[i]-s)/a[i][i];
    }
    rechercheZero(a,b,n);
    printf("\n\t\tMatrice apres triangularisation :");
    afficheSysteme(a,b,n);
    printf("\n\t\tResultat apres la remontee :\n\n");
    for (i=0; i<n; i++) printf("\t\tX[%d] = %.2f ;\n",i+1,x[i]);
}

//Gauss avec pivot partiel
void gaussPivot(float A[19][19],float B[19],int n)
{
    system("cls");
    printf("\n\t\t              ---------------------------------------------------           \n");
    printf("\t\t             ---------------------------------------------------           \n");
    printf("\t\t              | RESOLUTION PAR LA METHODE DE GAUSS PIVOT PARTIEL |           \n");
    printf("\t\t             ---------------------------------------------------           \n");
    printf("\t\t              ---------------------------------------------------           \n");

    float a[19][19],b[19];
    int cpt,cpt1;
    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }

    float x[19],p,s,ref,temp;
    int i,j,k,ligne;

    for(k=0; k<n; k++)
    {
    // pivot maximum
        ref=0;
        for(i=k; i<n; i++) if(fabs(a[i][k])>ref)
            {
                ref=fabs(a[i][k]);
                ligne=i;
            }

    // pivot
        for(j=k; j<n; j++)
        {
            temp=a[k][j];
            a[k][j]=a[ligne][j] ;
            a[ligne][j]=temp;
        }

        temp=b[k];
        b[k]=b[ligne];
        b[ligne]=temp;

        if (a[k][k]==0)
        {
            printf("\n\n\t ??? Un des pivots est nul alors la methode de Gauss pivot partiel est non applicable, utilisez une autre methode...\n\n");
            exit(10);
        }

        //triangularisation
        for(i=k+1; i<n; i++)
        {
            p=a[i][k]/a[k][k];
            for (j=k; j<n; j++) a[i][j]=a[i][j]-p*a[k][j];
            b[i]=b[i]-p*b[k];
        }
    }

    //Remontee
    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++) s=s+a[i][j]*x[j];
        x[i]=(b[i]-s)/a[i][i];
    }
    rechercheZero(a,b,n);
    printf("\n\t\tMatrice apres triangularisation :");
    afficheSysteme(a,b,n);
    printf("\n\t\tResultat apres la remontee :\n\n");
    for (i=0; i<n; i++) printf("\t\tX[%d] = %.2f ;\n",i+1,x[i]);
    printf("\n");

}

//Gauss Jordan
void gaussJordan(float A[19][19],float B[19],int n)
{
    system("cls");
    printf("\n\t\t              ---------------------------------------------          \n");
    printf("\t\t              ---------------------------------------------           \n");
    printf("\t\t              | RESOLUTION PAR LA METHODE DE GAUSS JORDAN |           \n");
    printf("\t\t              ---------------------------------------------           \n");
    printf("\t\t              ---------------------------------------------          \n");

    float a[19][19],b[19];
    int cpt,cpt1;
    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }

    float p;
    int i,j,k;

    for(k=0; k<n; k++)
    {
        if (a[k][k]==0)
        {
            printf("\n\n\t ??? Un des pivots est nul alors la methode de Gauss jordan est non applicable, utilisez une autre methode...\n\n");
            exit(10);
        }

        p=a[k][k];

        //normalisation
        for (j=k; j<n; j++) a[k][j]=a[k][j]/p;
        b[k]=b[k]/p;

        //Remontee
        for(i=0; i<n; i++)
        {
            if (i!=k)
            {
                p=a[i][k];
                for (j=k; j<n; j++) a[i][j]=a[i][j]-p*a[k][j];
                b[i]=b[i]-p*b[k];
            }
        }
    }
    rechercheZero(a,b,n);
    printf("\n\t\tMatrice apres triangularisation :");
    afficheSysteme(a,b,n);
    printf("\n\t\tResultat apres la remontee :\n\n");
    for(i=0; i<n; i++) printf("\t\tX[%d] = %.2f ;\n",i+1,b[i]);
    printf("\n");

}

//Crout
void crout(float A[19][19],float B[19],int n)
{
    system("cls");
    printf("\n\t\t             ----------------------------------------          \n");
    printf("\t\t              ----------------------------------------           \n");
    printf("\t\t              | RESOLUTION PAR LA METHODE DE LU CROUT |           \n");
    printf("\t\t              ----------------------------------------           \n");
    printf("\t\t              ----------------------------------------         \n");

    float a[19][19],b[19];
    int cpt,cpt1;
    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }

    float L[19][19],U[19][19],x[19],y[19],s;
    int i,j,k,m;

    for (i=0; i<n; i++) for (j=0; j<n; j++)
        {
            if(i==j) U[i][j]=1;
            else U[i][j]=0;
            L[i][j]=0;
        }

    for (m=0; m<n; m++)
    {
        for (i=m; i<n; i++)
        {
            s=0;
            for (k=0; k<m; k++) s=s+L[i][k]*U[k][m];
            L[i][m]=a[i][m]-s;
        }

        if (L[k][k]==0)
        {
            printf("\n\n\t\t* Un mineur nul ! => methode de LU n'est pas applicable. Utilisez une autre methode.\n\n");
            exit(10);
        }

        for (j=m+1; j<n; j++)
        {
            s=0;
            for (k=0; k<m; k++) s=s+L[m][k]*U[k][j];
            U[m][j]=(a[m][j]-s)/L[m][m];
        }
    }

    // Resolution
    for(i=0; i<n; i++)
    {
        s=0;
        for(j=0; j<i; j++) s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++)s=s+U[i][j]*x[j];
        x[i]=(y[i]-s)/U[i][i];
    }

    printf("\n\t\t* A = L * U \n");
    printf("\n\t\t* Matrice Lower L :");
    afficheMatrice(L,n);
    printf("\n\t\t* Matrice Upper U :");
    afficheMatrice(U,n);
    printf("\n\t\t* Resultat :\n\n");
    for (i=0; i<n; i++) printf("\t\tX[%d] = %.2f ;\n",i+1,x[i]);
}

//Doolittle
void doolittle(float A[19][19],float B[19],int n)
{
    system("cls");
    printf("\n\t\t              ---------------------------------------------        \n");
    printf("\t\t              ---------------------------------------------           \n");
    printf("\t\t              | RESOLUTION PAR LA METHODE DE LU DOOLITTLE |           \n");
    printf("\t\t              ---------------------------------------------           \n");
    printf("\t\t              ---------------------------------------------          \n");

    float a[19][19],b[19];
    int cpt,cpt1;

    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }

    float L[19][19],U[19][19],x[19],y[19],s;
    int i,j,k,m;

    for (i=0; i<n; i++) for (j=0; j<n; j++)
        {
            if(i==j) L[i][j]=1;
            else L[i][j]=0;
            U[i][j]=0;
        }

    for (m=0; m<n; m++)
    {
        for (j=m; j<n; j++)
        {
            s=0;
            for (k=0; k<m; k++) s=s+L[m][k]*U[k][j];
            U[m][j]=a[m][j]-s;
        }
        if (U[k][k]==0)
        {
            printf("\n\n\t\tUn mineur nul ! => methode de LU non applicable\n\n");
            exit(10);
        }

        for (i=m+1; i<n; i++)
        {
            s=0;
            for (k=0; k<m; k++) s=s+L[i][k]*U[k][m];
            L[i][m]=(a[i][m]-s)/U[m][m];
        }
    }

    // resolution
    for(i=0; i<n; i++)
    {
        s=0;
        for(j=0; j<i; j++) s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++)s=s+U[i][j]*x[j];
        x[i]=(y[i]-s)/U[i][i];
    }

    printf("\n\t\t * A = L * U \n");
    printf("\n\t\t* La matrice Lower L :");
    afficheMatrice(L,n);
    printf("\n\t\t* La matrice Upper U :");
    afficheMatrice(U,n);
    printf("\n\n\t\tResultat :\n\n");
    for (i=0; i<n; i++) printf("\t\tX[%d] = %f ;\n",i+1,x[i]);
}


//Cholesky
void cholesky(float A[19][19],float B[19],int n)
{
    system("cls");
    printf("\n\t\t             -----------------------------------------          \n");
    printf("\t\t              -----------------------------------------           \n");
    printf("\t\t              | RESOLUTION PAR LA METHODE DE CHOLESKY |           \n");
    printf("\t\t              -----------------------------------------           \n");
    printf("\t\t              -----------------------------------------         \n");

    float a[19][19],b[19];
    int cpt,cpt1;
    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }

    float L[19][19],Lt[19][19],x[19],y[19],s,p;
    int i,j,k;

    // verification de le symetrie
    for (i=0; i<n; i++) for (j=0; j<n; j++)
            if (a[i][j]!=a[j][i])
            {
                printf("\n\n\t\t* Matrice non symetrique => methode de Cholesky non applicable; saisir une autre matrice ou changer de methode\n\n");
                exit(10);
            }

    for (i=0; i<n; i++) for (j=0; j<n; j++) L[i][j]=0;

    for (i=0; i<n; i++)
    {
        s=0;
        for (k=0; k<i; k++) s=s+pow(L[i][k],2);
        p=a[i][i]-s;

        if (p<=0)
        {
            printf("\n\n\t\t* Matrice non definie positive => methode de Cholesky non applicable; saisir une autre matrice ou changer de methode\n\n");
            exit(10);
        }

        L[i][i]=sqrt(p);

        for(j=i+1; j<n; j++)
        {
            s=0;
            for (k=0; k<i; k++) s=s+L[i][k]*L[j][k];
            L[j][i]=(a[j][i]-s)/L[i][i];
        }
    }

    for (i=0; i<n; i++) for (j=0; j<n; j++) Lt[i][j]=L[j][i];

    // resolution
    for(i=0; i<n; i++)
    {
        s=0;
        for(j=0; j<i; j++) s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++) s=s+Lt[i][j]*x[j];
        x[i]=(y[i]-s)/Lt[i][i];
    }

    printf("\n\t\t A = L * Lt \n");
    printf("\n\t\t La matrice L :");
    afficheMatrice(L,n);
    printf("\n\t\t La matrice Lt :");
    afficheMatrice(Lt,n);
    printf("\n\t\t* Resultat :\n\n");
    for (i=0; i<n; i++) printf("\t\tX[%d] = %.2f ;\n",i+1,x[i]);
}

//Jacobi
void jacobie(float A[19][19],float B[19],int n)
{
    system("cls");
    printf("\n\t\t*            ----------------------------------------          *\n");
    printf("\t\t*              ----------------------------------------           *\n");
    printf("\t\t*              | RESOLUTION PAR LA METHODE DE JACOBIE |           *\n");
    printf("\t\t*              ----------------------------------------           *\n");
    printf("\t\t*              ----------------------------------------          *\n\n");

    float a[19][19],b[19];
    int cpt,cpt1;
    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }

    float x[19],x1[19],x2[19],s,eps=1e-4;
    int i,j,k,iter=0;

    //saisie x0
    printf("\n\t\t  Veuillez saisir le vecteur solution X0 : \n\n");
    for (i=0; i<n; i++)
    {
        printf("\t\tX(0)[%d]= ",i+1);
        scanf("%f",&x1[i]);
    }

    do
    {
        for(i=0; i<n; i++)
        {
            s=0;
            for (j=0; j<n; j++) if (i!=j) s=s+a[i][j]*x1[j];
            x2[i]=(b[i]-s)/a[i][i];
        }
        for (k=0; k<n; k++)
        {
            x[k]=fabs(x1[k]-x2[k]);
            x1[k]=x2[k];
        }

        iter++;
    }
    while (norme(x,n)>eps) ;

    printf("\n\t\t* Resultat :\n\n");
    for (i=0; i<n; i++) printf("\t\tX[%d] = %f ;\n",i+1,x2[i]);
    printf("\n\t\t* %d iterations,  10^-4 pres \n",iter);
}

//Gauss Seidel
void gaussSeidel(float A[19][19],float B[19],int n)
{
    system("cls");
    printf("\n\t\t*            --------------------------------------------      *\n");
    printf("\t\t*              --------------------------------------------           *\n");
    printf("\t\t*              | RESOLUTION PAR LA METHODE DE GAUSS SEIDEL |           *\n");
    printf("\t\t*              --------------------------------------------           *\n");
    printf("\t\t*              ---------------------------------------------      *\n\n");

    float a[19][19],b[19];
    int cpt,cpt1;
    for(cpt=0; cpt<n; cpt++)
    {
        for(cpt1=0; cpt1<n; cpt1++)
        {
            a[cpt][cpt1] = A[cpt][cpt1];
        }
    }

    for(cpt=0; cpt<n; cpt++)
    {
        b[cpt] = B[cpt];
    }

    float x[19],x1[19],x2[19],s,p,eps=1e-4;
    int i,j,k,iter=0;

    //initialisation du vecteur
    printf("\n\t\t* Veuillez saisir le vecteur solution X0 : \n\n");
    for (i=0; i<n; i++)
    {
        printf("\t\tX(0)[%d]= ",i+1);
        scanf("%f",&x1[i]);
    }

    do
    {
        for(i=0; i<n; i++)
        {
            s=0;
            p=0;
            for (j=i+1; j<n; j++) s=s+a[i][j]*x1[j];
            for (j=0; j<i; j++) p=p+a[i][j]*x2[j];
            x2[i]=(b[i]-s-p)/a[i][i];
        }
        for (k=0; k<n; k++)
        {
            x[k]=fabs(x1[k]-x2[k]);
            x1[k]=x2[k];
        }

        iter++;
    }
    while (norme(x,n)>eps) ;

    printf("\n\t\t* Resultat :\n\n");
    for (i=0; i<n; i++) printf(" X_%d = %.2f ;\n",i+1,x2[i]);
    printf("\n\t\t %d iterationS, 10^-4 pres. \n",iter);
}


///*********************************************************************
///************************************************                     *********
///FIN S-E-L                                                               *********
///************************************************                     *********
///*********************************************************************

///*********************************************************************
///************************************************                     *********
///Definition de fonctions d interpolation lineaire                            *********
///************************************************                     *********
///*********************************************************************

float* AllocatVect(int n)
{
    float* Vect = (float*) malloc(sizeof(float)*n);
    if(Vect == NULL)
    {
        printf("Memoire insuffisante pour continuer!\n");
        exit(EXIT_FAILURE);
    }
    return Vect;
}

float** AllocatMat(int ligne, int colonne)
{
    float** Mat = malloc(sizeof(float*)*ligne);
    int i;
    for(i = 0; i<ligne; i++)
    {
        Mat[i] = (float*)malloc(sizeof(float)*colonne);
        if(Mat[i]==NULL)
        {
            printf("Memoire insuffisante pour continuer !");
            exit(EXIT_FAILURE);
        }
    }
    return Mat;
}

void SaisieDeDonnees (int n,float* x, float* y)
{
    int i;
    for(i=0; i<n; i++)
    {
        printf("Saisir x%d:\t",i);
        x[i] = ctrsaisiefloat();
        printf("Saisir y%d:\t",i);
        y[i] = ctrsaisiefloat();
        printf("\n");
    }

    printf("\n\n\t Coordonn�es saisie :\n\n");
    for(i=0; i<n; i++)
    {
        printf("\n\tx%d = %2.2f\t|\ty%d = %2.2f\n\n-------------------------------------------------------\n", i, x[i], i, y[i]);
    }
}

float* choleskyInter(float** Mat,float* b,int n)
{
    float** L=AllocatMat(n,n);
    float** Lt=AllocatMat(n,n);
    float* x=AllocatVect(n);
    float* y=AllocatVect(n);
    float s,p;
    int i,j,k;

    // v�rification de le sym�trie
    for (i=0; i<n; i++) for (j=0; j<n; j++)
            if (Mat[i][j]!=Mat[j][i])
            {
                printf("\n\n * La matrice saisie n est pas symetrique ,la methode n est donc applicable\n");
                exit(EXIT_FAILURE);
            }

    for (i=0; i<n; i++) for (j=0; j<n; j++) L[i][j]=0;

    for (i=0; i<n; i++)
    {
        s=0;
        for (k=0; k<i; k++) s=s+pow(L[i][k],2);
        p=Mat[i][i]-s;

        if (p<=0)
        {
            printf("\n\n * LA matrice saisie n est pas definie positive par consequent la methode de Cholesky n est pas applicable\n\n");
            exit(EXIT_FAILURE);
        }

        L[i][i]=sqrt(p);

        for(j=i+1; j<n; j++)
        {
            s=0;
            for (k=0; k<i; k++) s=s+L[i][k]*L[j][k];
            L[j][i]=(Mat[j][i]-s)/L[i][i];
        }
    }

    for (i=0; i<n; i++) for (j=0; j<n; j++) Lt[i][j]=L[j][i];

    // resolution
    for(i=0; i<n; i++)
    {
        s=0;

        for(j=0; j<i; j++) s=s+L[i][j]*y[j];
        y[i]=(b[i]-s)/L[i][i];
    }

    for(i=n-1; i>=0; i--)
    {
        s=0;
        for(j=i+1; j<n; j++) s=s+Lt[i][j]*x[j];
        x[i]=(y[i]-s)/Lt[i][i];
    }
    return x;
}
int ctrsaisieint(void)
{

    int ok;
    char f[100];
    do
    {
        int i;
        fflush(stdin);
        ok=scanf("%[^\n]",f);
        int taille=strlen(f);
        fflush(stdin);
        for(i=0; i<taille; i++)
        {

            if(!isdigit(f[i]))
            {
                ok = 0;
                printf("Saisie invalide\n");
                break;
            }
        }
        fflush(stdin);
    }
    while(ok!=1);
    return (int) atoi(f);
}

float ctrsaisiefloat(void)
{
    int ok;
    char f[100];
    do
    {
        int point = 0, tiret = 0;
        int i;

        ok=scanf("%[^\n]",f);
        int taille=strlen(f);
        fflush(stdin);
        for(i=0; i<taille; i++)
        {

            if(!(isdigit(f[i]) || f[i]=='.' || f[i]=='f' || f[i]=='-'))
            {
                ok = 0;
                printf("Saisie invalide, ressaisir\n");
                break;
            }

            if(f[i]=='.')
            {
                point++;
            }
            if(f[i]=='-')
            {
                tiret++;
            }

            if(point>1 || tiret>1 || !(isdigit(f[i]) || f[i]=='.' || f[i]=='f' || f[i]=='-'))
            {
                printf("Saisie invalide, ressaisir\n");
                ok=0;
                break;
            }
        }
        fflush(stdin);
    }
    while(ok!=1);
    return (float) atof(f);
}

float* CalculCoefDePhi(int n,float* x, float* y)
{
    float* tab = AllocatVect(n);
    int i,j;
    for(i=0; i<n; i++)
    {
        tab[i]=1;
        for(j=0; j<n; j++)
        {
            if(i!=j)
            {
                tab[i] *= (x[i] - x[j]);
            }
        }
        tab[i] = y[i]/tab[i];

    }
    return tab;
}

float* CalculCoefDep(int n, float* x, float* y,int p)
{
    float** P = AllocatMat(p+1,p+1);
    float* A = AllocatVect(p+1);
    float* B = AllocatVect(p+1);
    int i, j, k;
    int q = 2*p, m, l = p;
    for(i=0;i<p+1;i++)
    {
        m = q;

        for(j=0;j<p+1;j++)
        {
            P[i][j] = 0;
            for(k=0;k<n;k++)
            {
                P[i][j] += pow(x[k],m);
            }
            m--;
        }
        q--;
    }
    for(i=0;i<p+1;i++)
    {
        B[i] = 0;
        for(k=0;k<n;k++)
        {
            B[i] += pow(x[k],l)*y[k];
        }
        l--;
    }
    A = choleskyInter(P, B, p+1);
    return A;
}

float** CalculDiffDiv(int n, float* x, float* y)
{
    float** tab = AllocatMat(n,n);
    int i,j;
    for(i = 0; i<n; i++)
    {
        tab[i][0] = y[i];
    }
    for(j = 1; j<n; j++)
    {
        for(i = j; i < n; i++)
        {
            tab[i][j] = (tab[i-1][j-1] - tab[i][j-1])/(x[0] - x[j]);
        }

    }

    return tab;
}

float* CalculCoefDeN(int n, float* x, float* y)
{
    float** diffDIV = AllocatMat(n,n);
    float* tab = AllocatVect(n);
    diffDIV = CalculDiffDiv(n, x, y);
    int i,j;
    for(j = 0; j<n; j++)
    {
        for(i = j; i < n; i++)
        {
            tab[i] = diffDIV[i][j];
        }
    }
    return tab;
}

//****************
///Methodes
//****************

///Lagrange
void Lagrange(int n,float* x, float* y)
{

    printf("\n\t\t-------------------------------------------------------\n");
    printf("\t\t-------------------------------------------------------\n");
    printf("\t\t  INTERPOLATION LINEAIRE PAR LA METHODE DE LAGRANGE   \n");
    printf("\t\t-------------------------------------------------------\n");
    printf("\t\t-------------------------------------------------------\n\n");

    float* CoefDePhi = CalculCoefDePhi(n,x,y);
    int i,j;
    printf(" \nPolynome d'interpolation de LAGRANGE :\n\n");
    printf("P%d(x)  = ",n);

    for(i=0; i<n; i++)
    {
        if(CoefDePhi[i]>0)
        {
            printf("(%.2f)",CoefDePhi[i]);
            for(j=0; j<n; j++)
            {
                if(i!=j)
                {
                    printf("(x-(%.2f))",x[j]);
                }

            }
            if(i!=(n-1))
            {
                printf("+");
            }
        }
        else
        {
            printf("(%.2f)",(-CoefDePhi[i]));
            for(j=0; j<n; j++)
            {
                if(i!=j)
                {
                    printf("(-x+(%.2f))",x[j]);
                }

            }
            if(i!=(n-1))
            {
                printf("+");
            }
        }
    }
}

//Newton
void Newton(int n, float* x, float* y)
{

    printf("\n\t\t-------------------------------------------------------\n");
    printf("\t\t-------------------------------------------------------\n");
    printf("\t\t  INTERPOLATION LINEAIRE PAR LA METHODE DE NEWTON   \n");
    printf("\t\t-------------------------------------------------------\n");
    printf("\t\t-------------------------------------------------------\n\n");
    float* tabAlpha = AllocatVect(n);
    tabAlpha = CalculCoefDeN(n, x, y);
    int i,j;
    printf("\nLe Polynome d'interpolation :\n\n");
    printf("P(x)=(%.2f)+",tabAlpha[0]);
    for(i = 1; i < n; i++)
    {
        if(tabAlpha[i]>0)
        {
            printf("(%.2f)", tabAlpha[i]);
            for(j = 0; j < i; j++)
            {
                printf("(x-(%.2f))",x[j]);
            }
            if(i!=(n-1))
            {
                printf("+");
            }
        }
        else
        {
            printf("%.2f", tabAlpha[i]);
            for(j = 0; j < i; j++)
            {
                printf("(-x+(%.2f))",-x[j]);
            }
            if(i!=(n-1))
            {
                printf("+");
            }
        }

    }
}

//Moindres carres
void Moindre(int n, float* x, float* y)
{
    printf("\n\t\t------------------------------------------------------------\n");
    printf("\t\t-------------------------------------------------------------\n");
    printf("\t\t  INTERPOLATION LINEAIRE PAR LA METHODE DE MOINDRES CARREES   \n");
    printf("\t\t------------------------------------------------------------\n");
    printf("\t\t------------------------------------------------------------\n\n");

    int i,p,puis;
    do
    {
        printf("Saisir le degre du polynome\n");
        p = ctrsaisieint();
        if(p<1 || p>n)printf("p ne peut prendre cette valeur!");
    }while(p<1 || p>n);
    puis = p;
    float* A = CalculCoefDep(n, x, y,p);
    printf("\nLe Polynome d'interpolation des moindres carrees: \n");
    printf("\tP(x) =");
    for(i=0;i<p+1;i++)
    {
        printf("%.2fx^%d+", A[i], puis);
        puis--;
    }
}

///*********************************************************************
///************************************************                     *********
///FIN Interpolation-L                                                      *********
///************************************************                     *********
///*********************************************************************


///*********************************************************************
///*************************************************** *           *********
///Definition de fonctions d'equations differentielles     ************************
///**************************************************** *          *********
///*********************************************************************

float SuiteDeleur(float h, float xn, float yn)
{
    return yn + h*FonctionDerivee(xn,yn);
}

float SuiteDeKunta(float h, float xn, float yn)
{
    return yn + (h/6)*((FonctionDerivee(xn,yn)) + 2*(FonctionDerivee(xn+h/2, yn + (h/2)*FonctionDerivee(xn,yn))));
}

float FonctionDerivee(float x, float y)
{
    return pow(x,2) - pow(y,2);
}

//****************
///Methodes
//****************

//Euler
void TraitementEuler(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h)
{

        printf("\n\t\t--------------------------------------------\n");
        printf("\t\t----------------------------------------------\n");
        printf("\t\t  EQUATION LINEAIRE PAR LA METHODE DE EULER   \n");
        printf("\t\t----------------------------------------------\n");
        printf("\t\t----------------------------------------------\n\n");

    int i;
    TableauDesXi[0] = xo;
    TableauDesYi[0] = yo;
    for(i = 0; i < n; i++)
    {
        TableauDesXi[i+1] = h + TableauDesXi[i];
        TableauDesYi[i+1] = SuiteDeleur(h, TableauDesXi[i], TableauDesYi[i]);
    }
    printf("\nXi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesXi[i]);
    }
    printf("\n\t--------------------------------------------------\n");
    printf("Yi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesYi[i]);
    }
}

//Kunta
void TraitementKunta(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h)
{

    printf("\n\t\t------------------------------------------------\n");
    printf("\t\t--------------------------------------------------\n");
    printf("\t\t  EQUATION LINEAIRE PAR LA METHODE DE RUNGE KUNTA   \n");
    printf("\t\t--------------------------------------------------\n");
    printf("\t\t--------------------------------------------------\n\n");

    int i;
    TableauDesXi[0] = xo;
    TableauDesYi[0] = yo;
    for(i = 0; i < n; i++)
    {
        TableauDesXi[i+1] = h + TableauDesXi[i];
        TableauDesYi[i+1] = SuiteDeleur(h, TableauDesXi[i], TableauDesYi[i]);
    }
    printf("\nXi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesXi[i]);

    }
    printf("\n\t--------------------------------------------------\n");
    printf("Yi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesYi[i]);
    }
}

//Euler modifie
void TraitementEulerModifie(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h)
{
    printf("\n\t\t--------------------------------------------------\n");
    printf("\t\t----------------------------------------------------\n");
    printf("\t\t  EQUATION LINEAIRE PAR LA METHODE DE EULER MODIFIEE   \n");
    printf("\t\t----------------------------------------------------\n");
    printf("\t\t----------------------------------------------------\n\n");

    int i;
    TableauDesXi[0] = xo;
    TableauDesYi[0] = yo;
    for(i = 0; i < n; i++)
    {
        TableauDesXi[i+1] = h + TableauDesXi[i];
        TableauDesYi[i+1] = SuiteDeKunta(h, TableauDesXi[i], TableauDesYi[i]);
    }
    printf("\nXi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesXi[i]);

    }
    printf("\n\t--------------------------------------------------\n");
    printf("Yi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesYi[i]);
    }
}

//Kunta ordre 4
void TraitementKuntaOrdre4(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h)
{
    printf("\n\t\t--------------------------------------------------------\n");
    printf("\t\t----------------------------------------------------------\n");
    printf("\t\t  EQUATION LINEAIRE PAR LA METHODE DE RUNGE KUNTA ORDRE 4   \n");
    printf("\t\t----------------------------------------------------------\n");
    printf("\t\t----------------------------------------------------------\n\n");

    int i;
    TableauDesXi[0] = xo;
    TableauDesYi[0] = yo;
    h = h/2;
    for(i = 0; i < n; i++)
    {
        TableauDesXi[i+1] = h + TableauDesXi[i];
        TableauDesYi[i+1] = SuiteDeleur(h, TableauDesXi[i], TableauDesYi[i]);
    }
    printf("\nXi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesXi[i]);

    }
    printf("\n\t--------------------------------------------------\n");
    printf("Yi\t|");
    for(i = 0; i <n; i++)
    {
        printf("\t%.2f |",TableauDesYi[i]);
    }

}

///*********************************************************************
///************************************************                     *********
///FIN equation-differentiel                                                               *********
///************************************************                     *********
///*********************************************************************



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
    int choix_met, retour;
    char rep;
    system("cls");
    printf("\t\t*       LES METHODES DE RESOLUTION DES SYSTEMES EQUATIONS LINEAIRES      *\n\n");
    do
    {
        printf("\n\n\t\tVeuillez saisir le nombre de ligne de la matrice : ");
        retour = scanf("%d",&n);
        fflush(stdin);
        if(retour==0) printf("\t\tLe nombre de ligne doit etre un reel. Veuillez ressaisir : ");
    }
    while(retour == 0);

    saisirMatrice(A,B,n);
    printf("\n\n\t\tAffichage du systeme : ");
    afficheSysteme(A,B,n);
    do
    {
        printf("\n\t\t*       LES METHODES DE RESOLUTION DES SYSTEMES D'EQUATIONS LINEAIRES      *\n\n");
        printf("\n\t\t\t1- Gauss sans pivot");
        printf("\n\t\t\t2- Gauss avec pivot");
        printf("\n\t\t\t3- Gauss Jordan");
        printf("\n\t\t\t4- CROUT");
        printf("\n\t\t\t5- Doolittle");
        printf("\n\t\t\t6- Cholesky");
        printf("\n\t\t\t7- Jacobi");
        printf("\n\t\t\t8- Gauss Seidel");

        do
        {
            printf("\n\n\t\tVeuillez choisir une methode : ");
            scanf("%d", &choix_met);
            fflush(stdin);
        }
        while( choix_met < 1 || choix_met > 8);

        switch(choix_met)
        {
        case 1 :
            gauss(A,B,n);
            break;
        case 2 :
            gaussPivot(A,B,n);
            break;
        case 3 :
            gaussJordan(A,B,n);
            break;
        case 4 :
            crout(A,B,n);
            break;
        case 5 :
            doolittle(A,B,n);
            break;
        case 6 :
            cholesky(A,B,n);
            break;
        case 7 :
            jacobie(A,B,n);
            break;
        case 8 :
            gaussSeidel(A,B,n);
            break;
        }

        do
        {
            fflush(stdin);
            printf("\n\n\t\tVoulez-vous revenir au menu des systemes d equations non-lineaires (O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tSaisissez o/O pour Oui ou n/N pour Non : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }
        while(rep != 'O' && rep != 'N');
        system("cls");
    }
    while(rep=='O');
}

///**********main() - intepolation lineaire
void interpolation()
{
    int n, choix_met, retour;
    char rep;
    system("cls");
    setlocale(LC_CTYPE,"");
    printf("\n\t\t      INTERPOLATION LINEAIRES      \n\n");
    do
    {
        printf("\n\n\t\tVeuillez saisir le degre du systeme : ");
        retour = scanf("%d",&n);
        fflush(stdin);
        if(retour==0) printf("\t\tSaisir une valeur reelle : ");
    }
    while(retour == 0);
    n++;
    float* x = AllocatVect(n);//tableau de xi
    float* y = AllocatVect(n);//tableau de yi
    SaisieDeDonnees(n,x,y);
    fflush(stdin);
    do
    {
        printf("\n\t\t*       LES METHODES D INTERPOLATIONS LINEAIRES      *\n");
        printf("\n\t\t\t1- METHODE DE LAGRANGE");
        printf("\n\t\t\t2- METHODE DE NEWTON");
        printf("\n\t\t\t3- METHODE DES MOINDRES CARREES");

        do
        {
            printf("\n\n\t\tVeuillez choisir une methode : ");
            scanf("%d", &choix_met);
            fflush(stdin);
        }
        while( choix_met < 1 || choix_met > 3);

        switch(choix_met)
        {
        case 1 :
            Lagrange(n, x, y);
            break;
        case 2 :
            Newton(n, x, y);
            break;
        case 3 :
            Moindre(n, x, y);
            break;
        }

        do
        {
            fflush(stdin);
            printf("\n\n\t\tVoulez-vous revenir au menu des interpolationslineaires (O/N) ? : ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tSaisissez o/O pour Oui ou n/N pour Non : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }
        while(rep != 'O' && rep != 'N');
        system("cls");
    }
    while(rep=='O');
}

///**********main() - equation differentielle
void equaDif()
{
    int n;//dimension du systeme
    int choix_met, retour;
    char rep;
    float yo, xo, h;
    float* TableauDesXi;
    float* TableauDesYi;
    system("cls");
    setlocale(LC_CTYPE,"");
    printf("\t\t\t\tEQUATIONS DIFFERENTIELLES\n");
    printf("\t\t\t\t-------------------------\n");

    do
    {
        printf("\n\n\t\tVeuillez saisir le degre du systeme : ");
        retour = scanf("%d",&n);
        fflush(stdin);
        if(retour==0) printf("\t\tSaisir une valeur reelle : ");
    }
    while(retour == 0);
    n++;
    TableauDesXi = AllocatVect(n);
    TableauDesYi = AllocatVect(n);
    printf("\nSaisir yo: \n");
    yo = ctrsaisiefloat();
    printf("\nSaisir xo: \n");
    xo = ctrsaisiefloat();
    printf("\n\nSaisir le pas h: \n");
    h = ctrsaisiefloat();

    do
    {
        printf("\n\t\t*       LES METHODES D EQUATIONS DIFFERENTIELLES      *\n");
        printf("\n\t\t\t1- METHODE DE EULER");
        printf("\n\t\t\t2- METHODE DE RUNGE KUNTA");
        printf("\n\t\t\t3- METHODE DE EULER MODIFIE");
        printf("\n\t\t\t4- METHODE DE RUNGE KUNTA ORDRE 4");

        do
        {
            printf("\n\n\t\tVeuillez choisir une methode : ");
            scanf("%d", &choix_met);
            fflush(stdin);
        }
        while( choix_met < 1 || choix_met > 4);

        switch(choix_met)
        {
        case 1 :
            TraitementEuler(TableauDesXi, TableauDesYi, n, xo, yo, h);
            break;
        case 2 :
            TraitementKunta(TableauDesXi, TableauDesYi, n, xo, yo, h);
            break;
        case 3 :
            TraitementEulerModifie(TableauDesXi, TableauDesYi, n, xo, yo, h);
            break;
        case 4 :
            TraitementKuntaOrdre4(TableauDesXi, TableauDesYi, n, xo, yo, h);
            break;
        }

        do
        {
            fflush(stdin);
            printf("\n\n\t\tVoulez-vous revenir au menu des equations differentielles (O/N) ?  ");
            scanf("%c", &rep);
            fflush(stdin);
            rep = toupper(rep);
            while(rep != 'O' && rep != 'N')
            {
                printf("\t\tSaisissez o/O pour Oui ou n/N pour Non : ");
                scanf("%c", &rep);
                fflush(stdin);
                rep = toupper(rep);
            }
        }
        while(rep != 'O' && rep != 'N');
        system("cls");
    }
    while(rep=='O');
}

///*************************************************************************************************************************
///*************************************************************************************************************************
///Main principal-----------------------------------------------------------------------------------------------------------
///*************************************************************************************************************************
///*************************************************************************************************************************



int main()
{
    setlocale(LC_CTYPE,"");
    //Variables
    char choix_ini;
    char rep;
    do{
    system("cls");
    printf("\n\t\t*                LES METHODES NUMERIQUES                  *\n");
    printf("\n\t\t\tA- Equation non lineaire ");
    printf("\n\t\t\tB- Systeme d'equation lineaire");
    printf("\n\t\t\tC- Interpolations");
    printf("\n\t\t\tD- Equations differentielles");

    do
    {
        printf("\n\n\t\tVotre choix : ");
        scanf("%c", &choix_ini);
        fflush(stdin);
        choix_ini = toupper(choix_ini);

    }
    while(choix_ini != 'A' && choix_ini != 'B' && choix_ini != 'C' && choix_ini != 'D');

    switch(choix_ini)
    {
    case 'A':
        equation_lineaire();
        break;
    case 'B':
        systeme_equation_lineaire();
        break;
    case 'C':
        interpolation();
        break;
    case 'D':
        equaDif();
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
