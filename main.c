/*
Noé Aubin-Cadot.
Codé en novembre 2008.
Adapté pour SDL2 en juin 2017.
Déposé sur GitHub en septembre 2019.
C'est une simulation de l'équation d'onde 2D par une surface qui ondule longitudinalement.
*/
// Compiler avec :
// gcc -o main main.c -framework SDL2 -framework SDL2_image -framework OpenGL -framework GLUT -Wno-deprecated
// ou encore :
// gcc $(sdl2-config --cflags) -Wall -o main  main.c $(sdl2-config --libs) -lSDL2_image -framework OpenGL -framework GLUT -Wno-deprecated

// Pour C :
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

// Pour SDL :
#include <SDL2/SDL.h>
//#include <SDL2_image/SDL_image.h>
#include <SDL2/SDL_image.h>

// Pour OpenGL :
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>

void Dessiner_sin(double *EQondex,double *EQondey);

double angleX = 0;
double angleY = 0;
double angleZ = 0;
double vueX = 0.001;
double vueY = 0;
double vueZ = 26;
int L=25;
double ROT = 0;
double SINGULARITE=0.01;
double facteur=2;

// Type d'onde :
int type_onde = 3; // choix du type d'onde longitudinale :
// 1 : onde venant d'une singularité
// 2 : onde sinus
// 3 : onde spirale

SDL_Window *fenetre;
SDL_GLContext contexte;

int main(int argc, char *argv[])
{	
    // Initialisation de SDL
    SDL_Init(SDL_INIT_VIDEO);
    
    // Création d'une fenêtre
    fenetre = SDL_CreateWindow("Eq. d'onde longitudinale 2D", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 640, 480, SDL_WINDOW_OPENGL);
    contexte = SDL_GL_CreateContext(fenetre);

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective(70,(double)640/480,1,1000);

    glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
    
    Uint32 last_time = SDL_GetTicks();
    Uint32 current_time,ellapsed_time;
    Uint32 start_time;
	
	int i,j,k;
	int M=L;
	int N=L;
	int TOURNE=0;
	int ONDEgo=0;
	double FORCE=1;
	double E11,E12,E22,E21;
	int inerationPARplot=100;
	double eps=0.000005;
	double theta=M_PI/5;
	double *EQondex;EQondex=(double *) malloc((M*N)*sizeof(double));
	double *EQonde1x;EQonde1x=(double *) malloc((M*N)*sizeof(double));
	double *EQonde2x;EQonde2x=(double *) malloc((M*N)*sizeof(double));
	double *EQondey;EQondey=(double *) malloc((M*N)*sizeof(double));
	double *EQonde1y;EQonde1y=(double *) malloc((M*N)*sizeof(double));
	double *EQonde2y;EQonde2y=(double *) malloc((M*N)*sizeof(double));
	
	// conditions initiales de l'onde
	if (type_onde==1)
	{
		// ONDE UNE SINGULARITE
		for (i=0;i<M;i++)
			{
			for (j=0;j<N;j++)
				{
				EQondex[i*M+j]=(i-M/2)/facteur;
				EQonde1x[i*M+j]=(i-M/2)/facteur;
				EQonde2x[i*M+j]=(i-M/2)/facteur;
				EQondey[i*M+j]=(j-N/2)/facteur;
				EQonde1y[i*M+j]=(j-N/2)/facteur;
				EQonde2y[i*M+j]=(j-N/2)/facteur;
				}
			}
		EQondex[N*M/2]=SINGULARITE;
		EQondey[N*M/2]=SINGULARITE;
	}
	if (type_onde==2)
	{
		// ONDE SINUSOIDALE
		for (i=0;i<M;i++)
		{
			for (j=0;j<N;j++)
			{
				EQondex[i*M+j] =sin(2*M_PI*(i-M/2)/(M-1))*sin(2*M_PI*(j-N/2)/(N-1))+(i-M/2)/facteur;
				EQonde1x[i*M+j]=sin(2*M_PI*(i-M/2)/(M-1))*sin(2*M_PI*(j-N/2)/(N-1))+(i-M/2)/facteur;
				EQonde2x[i*M+j]=sin(2*M_PI*(i-M/2)/(M-1))*sin(2*M_PI*(j-N/2)/(N-1))+(i-M/2)/facteur;
				EQondey[i*M+j] =sin(2*M_PI*(i-M/2)/(M-1))*sin(2*M_PI*(j-N/2)/(N-1))+(j-N/2)/facteur;
				EQonde1y[i*M+j]=sin(2*M_PI*(i-M/2)/(M-1))*sin(2*M_PI*(j-N/2)/(N-1))+(j-N/2)/facteur;
				EQonde2y[i*M+j]=sin(2*M_PI*(i-M/2)/(M-1))*sin(2*M_PI*(j-N/2)/(N-1))+(j-N/2)/facteur;	
			}
		}
	}
	if (type_onde==3)
	{
		// ONDE SPIRALE
		for (i=0;i<M;i++)
		{
			for (j=0;j<N;j++)
			{
				double x=(double)2*i/(M-1)-1;
				double y=(double)2*j/(N-1)-1;
				double theta2=theta*cos(x*M_PI/2)*cos(y*M_PI/2);
				EQondex[i*M+j]=((i-M/2)*cos(theta2)+(j-N/2)*sin(theta2))/facteur;
				EQonde1x[i*M+j]=EQondex[i*M+j];
				EQonde2x[i*M+j]=EQondex[i*M+j];
				EQondey[i*M+j]=(-(i-M/2)*sin(theta2)+(j-N/2)*cos(theta2))/facteur;
				EQonde1y[i*M+j]=EQondey[i*M+j];
				EQonde2y[i*M+j]=EQondey[i*M+j];	
			}
		}
	}

    SDL_Event event;
    while (1)
    {
        start_time = SDL_GetTicks();
        while (SDL_PollEvent(&event))
        {
            switch(event.type)
            {
                case SDL_QUIT:
					exit(0);
					break;
				case SDL_KEYDOWN:
                switch(event.key.keysym.sym)
				{
                    case SDLK_ESCAPE: // Arreter le jeu
                        exit(0);
						break;
					case SDLK_SPACE: // Stop ou part la rotation du plan d'onde
                        {int EVENT=0;
						if (TOURNE)
							{TOURNE=0;EVENT=1;}
						if (!TOURNE && !EVENT)
							{TOURNE=1;}
						break;
						}
					case SDLK_RETURN: // Stop ou part le fait que l'onde evolue
                        {int EVENT=0;
						if (ONDEgo)
							{ONDEgo=0;EVENT=1;}
						if (!ONDEgo && !EVENT)
							{ONDEgo=1;}
						break;
						}
					case SDLK_UP:
						vueZ-=vueZ/10;
						break;
					case SDLK_DOWN:
						vueZ+=vueZ/10;
						break;
					case SDLK_RIGHT:
						angleZ-=30;
						break;
					case SDLK_LEFT:
						angleZ+=30;
						break;
					case SDLK_o:
						//FORCE--;
						FORCE-=0.1;
						break;
					case SDLK_p:
						//FORCE++;
						FORCE+=0.1;
						break;
					case SDLK_i:
						ROT+=2;
						break;
					case SDLK_k:
						ROT-=2;
						break;
				}
                break;
            }
        }

        current_time = SDL_GetTicks();
        ellapsed_time = current_time - last_time;
        last_time = current_time;
		
		if (TOURNE)
			{
			angleZ += 0.01 * ellapsed_time;
			}
        
        Dessiner_sin(EQondex,EQondey);
		if (ONDEgo)
			{
			for (k=0;k<inerationPARplot;k++)
			{
				for (i=1;i<M-1;i++)
					{
					for (j=1;j<N-1;j++)
						{
						EQonde2x[i*M+j]=EQonde1x[i*M+j];
						EQonde1x[i*M+j]=EQondex[i*M+j];
						EQonde2y[i*M+j]=EQonde1y[i*M+j];
						EQonde1y[i*M+j]=EQondey[i*M+j];
						}
					}
				for (i=1;i<M-1;i++)
					{
					for (j=1;j<N-1;j++)
						{
						// EN X
						E11=EQonde1x[(i+1)*M+j] - EQonde1x[i*M+j];
						E12=EQonde1x[(i-1)*M+j] - EQonde1x[i*M+j];
						E22=EQonde1x[i*M+j+1] - EQonde1x[i*M+j];
						E21=EQonde1x[i*M+j-1] - EQonde1x[i*M+j];
						EQondex[i*M+j]=2*EQonde1x[i*M+j]-EQonde2x[i*M+j]+eps*(
												 (E11) * pow(E11*E11+1,(FORCE-1)/2)
												+(E12) * pow(E12*E12+1,(FORCE-1)/2)
												+(E22) * pow(E22*E22+1,(FORCE-1)/2)
												+(E21) * pow(E21*E21+1,(FORCE-1)/2)
												);
						// EN Y
						E11=EQonde1y[(i+1)*M+j] - EQonde1y[i*M+j];
						E12=EQonde1y[(i-1)*M+j] - EQonde1y[i*M+j];
						E22=EQonde1y[i*M+j+1] - EQonde1y[i*M+j];
						E21=EQonde1y[i*M+j-1] - EQonde1y[i*M+j];
						EQondey[i*M+j]=2*EQonde1y[i*M+j]-EQonde2y[i*M+j]+eps*(
															(E11) * pow(E11*E11+1,(FORCE-1)/2)
															+(E12) * pow(E12*E12+1,(FORCE-1)/2)
															+(E22) * pow(E22*E22+1,(FORCE-1)/2)
															+(E21) * pow(E21*E21+1,(FORCE-1)/2)
															);
						}
					}
				}
			}

        ellapsed_time = SDL_GetTicks() - start_time;
        if (ellapsed_time < 10)
        {
            SDL_Delay(10 - ellapsed_time);
        }
    }
    return 0;
}

void Dessiner_sin(double *EQondex,double *EQondey)
{
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity( );

	int M=L;
	int N=L;
    gluLookAt(vueX,vueY,vueZ,0,0,ROT,0,0,1);

    glRotated(angleZ,0,0,1);
    glRotated(angleY,0,1,0);
	glRotated(angleX,1,0,0);
	
    glBegin(GL_QUADS);
	int i,j;
	for (i=0;i<M-1;i++)
		{
		for (j=0;j<N-1;j++)
			{
			glColor3ub(255,255*(1+pow(-1,i+j))/2,255*(1+pow(-1,i+j))/2);
			glVertex3d(facteur*EQondex[i*M+j],facteur*EQondey[i*M+j],0.0001);
			glVertex3d(facteur*EQondex[i*M+j+1],facteur*EQondey[i*M+j+1],0.0002);
			glVertex3d(facteur*EQondex[(i+1)*M+j+1],facteur*EQondey[(i+1)*M+j+1],0.0003);
			glVertex3d(facteur*EQondex[(i+1)*M+j],facteur*EQondey[(i+1)*M+j],0.0004);
			}
		}
	glEnd();
    glFlush();
    SDL_GL_SwapWindow(fenetre);
}