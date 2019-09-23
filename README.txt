2019-09-21

		☐ ÉQUATION D'ONDE 2D LONGITUDINALE ☐

######################################################################
######################################################################

		DESCRIPTION

C'est une simulation de l'équation d'onde 2D par une surface qui ondule longitudinalement.

Codé en novembre 2008.
Adapté pour SDL2 en juin 2017.
Déposé sur GitHub en septembre 2019.

######################################################################
######################################################################

		OS

Ça marche sur macOS Mojave 10.14.6, mais ça devrait se porter aisément à d'autres OS.

######################################################################
######################################################################

		COMPILATION

D'abord, installer Brew.
Ensuite, utiliser Brew pour installer SDL :
	brew install SDL2
	brew install SDL2_image

Le code devrait compiler avec :

gcc $(sdl2-config --cflags) -Wall -o main  main.c $(sdl2-config --libs) -lSDL2_image -framework OpenGL -framework GLUT -Wno-deprecated

######################################################################
######################################################################

		CONTRÔLES

Espace : active/désactive la rotation continue ou non
Entrer : active/désactive l'équation d'onde

Gauche/droite : joue avec l'orientation de la surface
Haut/bas : s'approcher/s'éloigner de la surface

######################################################################
######################################################################

Noé Aubin-Cadot, 2019.