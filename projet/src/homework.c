#include "fem.h"


double Hermite(double x,double h0,double h,double dist,double d0){
    double a = h0;
    double b = 0;
    double c = (-3*(h0-h))/(d0*d0);
    double d = (2*(h0-h))/(d0*d0*d0);
    return a+b*dist+c*dist*dist+d*dist*dist*dist;
}

double geoSize(double x, double y){
    
    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
    double d0_star = d0;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;
    double dist_hole = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)) - r1;
    double dist_notch = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)) - r0;
    double d1_star = d1;
    return h;
    /*if(dist_notch < d0_star && dist_hole < d1_star){
        return fmin(Hermite(x0,h0,h,dist_notch,d0_star),Hermite(x1,h1,h,dist_hole,d1_star));
    }
    else if (dist_notch < d0_star && dist_hole > d1_star) {
    // Utilise hermiteInterpolation_Notch
        return Hermite(x0,h0,h,dist_notch,d0_star);
    }else if (dist_hole < d1_star && dist_notch > d0 ) {
    // Utilise hermiteInterpolation_Hole
        return Hermite(x1,h1,h,dist_hole,d1_star);
    }else {
        return h;
    }*/
}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    int ierr;
     theGeometry->h = 0.4; 

    // Ajout des points pour définir la géométrie
    double x = 0, y = 0, z = 0;
    gmshModelOccAddPoint(x, y, z, 1, 1, &ierr); // Augmenter le rayon pour réduire le nombre de points
    gmshModelOccAddPoint(x + 9, y, z, 1, 2, &ierr);
    gmshModelOccAddPoint(x + 9, y + 7, z, 1, 3, &ierr);
    gmshModelOccAddPoint(x+8, y + 7, z, 1, 4, &ierr);
    gmshModelOccAddPoint(x, y + 1, z, 1, 5, &ierr);

    gmshModelOccAddPoint(x + 1, y + 1, z, 1, 6, &ierr);
    gmshModelOccAddPoint(x + 4, y + 3.2, z, 1, 7, &ierr);
    gmshModelOccAddPoint(x + 4, y + 1, z, 1, 8, &ierr);

    gmshModelOccAddPoint(x + 5, y + 1, z, 1, 9, &ierr);
    gmshModelOccAddPoint(x + 5, y + 3.9, z, 1, 10, &ierr);
    gmshModelOccAddPoint(x + 8, y + 6, z, 1, 11, &ierr);
    gmshModelOccAddPoint(x + 8, y + 1, z, 1, 12, &ierr);
    // Ajouter d'autres points selon votre géométrie

    // Ajout des lignes pour former les bords de la surface
    gmshModelOccAddLine(1, 2, 1, &ierr);
    gmshModelOccAddLine(2, 3, 2, &ierr);
    gmshModelOccAddLine(3, 4, 3, &ierr);
    gmshModelOccAddLine(4, 5, 4, &ierr);
    gmshModelOccAddLine(5, 1, 5, &ierr);
    gmshModelOccAddLine(6, 7, 6, &ierr);
    gmshModelOccAddLine(7, 8, 7, &ierr);
    gmshModelOccAddLine(8, 6, 8, &ierr);
    gmshModelOccAddLine(9, 10, 9, &ierr);
    gmshModelOccAddLine(10, 11, 10, &ierr);
    gmshModelOccAddLine(11, 12, 11, &ierr);
    gmshModelOccAddLine(12, 9, 12, &ierr);
    // Ajouter d'autres lignes selon votre géométrie

    // Construction de la boucle de la surface
    int curveTags[] = {1,2,3,4,5}; // Utilisation de deux lignes pour former une boucle
    gmshModelOccAddCurveLoop(curveTags, 5, 1, &ierr);

    // Création de la surface
    int wireTags[] = {1}; // Utilisation de la boucle comme fil
    int outersurface = gmshModelOccAddPlaneSurface(wireTags, 1, 1, &ierr); 

    int innerCurveTags[] = {6, 7, 8}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags, 3, 2, &ierr);

    // Créer une nouvelle surface à l'intérieur
    int innerWireTags[] = {2}; // Utilisation de la nouvelle boucle comme fil
    int innersurface = gmshModelOccAddPlaneSurface(innerWireTags, 1, 2, &ierr);

    int innerCurveTags2[] = {9, 10, 11, 12}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags2, 4, 3, &ierr);

    int innerWireTags2[] = {3}; // Utilisation de la nouvelle boucle comme fil
    int innersurface2 = gmshModelOccAddPlaneSurface(innerWireTags2, 1, 3, &ierr);

    // Couper la surface extérieure par la surface intérieure
    int outer[] = {2,outersurface};
    int inner[] = {2,innersurface};
    gmshModelOccCut(outer,2,inner,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    int inner2[] = {2,innersurface2};
    gmshModelOccCut(outer,2,inner2,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    

    // Définition de la taille de maillage
    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);

    // Génération du maillage
    gmshModelMeshGenerate(2, &ierr);
    
}