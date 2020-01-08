#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <locale.h>

/// PROTOTYPE DES FONCTIONS

///Equations non-linéiares
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
//fonctions
void gauss();
void gaussPivot();
void gaussJordan();
void crout();
void doolittle();
void cholesky();
void jacobi();
void gaussSeidel();


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
            tolerence= saisirEntier("\t\tEntrer le critère d'arrêt : 10^-");
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
        printf("La convergence n'est pas atteinte après %d iterations", iteration);
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
        printf("La convergence n'est pas atteinte après %d iterations", iteration);
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
        printf("La convergence n'est pas atteinte après %d iterations", iteration);
        exit(10);
    }
}

///FIN E-N-L
///****************************************************************************************************************
///****************************************************************************************************************

///Definition de fonctions d'equations non-lineaire
//Gauss sans pivot
void gauss()
{
    printf("En cours...");
}

//Gauss avec pivot
void gaussPivot()
{
    printf("En cours...");
}

//Gauss Jordan
void gaussJordan()
{
    printf("En cours...");
}

//Crout
void crout()
{
    printf("En cours...");
}

//Doolittle
void doolittle()
{
    printf("En cours...");
}

//Cholesky
void cholesky()
{
    printf("En cours...");
}

//Jacobi
void jacobi()
{
    printf("En cours...");
}

//Gauss Seidel
void gaussSeidel()
{
    printf("En cours...");
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
    int choix_met;
    char rep;
    do{
             system("cls");

    printf("\t\t*       LES METHODES DE RESOLUTION DES SYSTEMES D'EQUATIONS LINEAIRES      *\n");

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
        printf("\n\n\t\tVotre choix : ");
        scanf("%d", &choix_met);
        fflush(stdin);
    }
    while( choix_met < 1 || choix_met > 8);

    switch(choix_met)
    {
    case 1 :
        gauss();
        break;
    case 2 :
        gaussPivot();
        break;
    case 3 :
        gaussJordan();
        break;
    case 4 :
        crout();
        break;
    case 5 :
        doolittle();
        break;
    case 6 :
        cholesky();
        break;
    case 7 :
        jacobi();
        break;
    case 8 :
        gaussSeidel();
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

