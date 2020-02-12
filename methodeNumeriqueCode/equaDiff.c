#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <locale.h>

///Equations différentielles
void equaDif();
//fonctions principales
void TraitementEuler(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h);
void TraitementKunta(float* TableauDesXi, float* TableauDesYi, int n, float xo, float yo,  float h);
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

int main()
{
        equaDif();
        return 0;
}

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

int Menu(void)
{
    int choix;

    printf("1: METHODE D'EULER\n2: METHODE DE RUNGE-KUNTA\n3: TERMINER!!!\n");
    do
    {
        printf("\n\nSaisir votre choix pour continuer\n");
        choix = ctrsaisieint();
        if(choix<1 || choix>3)printf("ressaisir a nouveau\n");
    }
    while(choix<1 || choix>3);

    return choix;
}
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
float FonctionDerivee(float x, float y)
{
    return pow(x,2) - pow(y,2);
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
                printf("Saisie invalide, ressaisir\n");
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
                printf("Saisie invalide, resaisir\n");
                ok=0;
                break;
            }
        }
        fflush(stdin);
    }
    while(ok!=1);
    return (float) atof(f);
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


///*********************************************************************
///************************************************                     *********
///FIN equation-differentiel                                                               *********
///************************************************                     *********
///*********************************************************************


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

        do
        {
            printf("\n\n\t\tVeuillez choisir une methode : ");
            scanf("%d", &choix_met);
            fflush(stdin);
        }
        while( choix_met < 1 || choix_met > 2);

        switch(choix_met)
        {
        case 1 :
            TraitementEuler(TableauDesXi, TableauDesYi, n, xo, yo, h);
            break;
        case 2 :
            TraitementKunta(TableauDesXi, TableauDesYi, n, xo, yo, h);
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