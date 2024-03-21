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
     theGeometry->h = 0.3; 

    // Ajout des points pour définir la géométrie
    double x = 0, y = 0, z = 0;
    gmshModelOccAddPoint(x, y, z, 1, 1, &ierr); // Augmenter le rayon pour réduire le nombre de points
    gmshModelOccAddPoint(x, y+1, z, 1, 2, &ierr);
    gmshModelOccAddPoint(x+0.58333, y + 5.5, z, 1, 3, &ierr);
    gmshModelOccAddPoint(x+7.5, y + 11 , z, 1, 4, &ierr);
    gmshModelOccAddPoint(x+14.41667, y + 5.5, z, 1, 5, &ierr);
    gmshModelOccAddPoint(x+15, y + 1, z, 1, 6, &ierr);
    gmshModelOccAddPoint(x+15, y, z, 1, 7, &ierr);

    gmshModelOccAddPoint(x+0.5, y, z, 1, 8, &ierr);
    gmshModelOccAddPoint(x+0.5, y+1, z, 1, 9, &ierr);
    gmshModelOccAddPoint(x+0.6622, y + 2.298, z, 1, 10, &ierr);
    gmshModelOccAddPoint(x+4.378, y + 5, z, 1, 11, &ierr);
    gmshModelOccAddPoint(x+10.622, y + 5, z, 1, 12, &ierr);
    gmshModelOccAddPoint(x+14.378, y + 2.298, z, 1, 13, &ierr);
    gmshModelOccAddPoint(x+14.5, y+1 , z, 1, 14, &ierr);
    gmshModelOccAddPoint(x+14.5, y, z, 1, 15, &ierr);

    gmshModelOccAddPoint(x+0.7499, y + 3, z, 1, 16, &ierr);
    gmshModelOccAddPoint(x+3.5, y + 5, z, 1, 17, &ierr);
    gmshModelOccAddPoint(x+1, y + 5, z, 1, 18, &ierr);

    gmshModelOccAddPoint(x+14.25, y+3 , z, 1, 19, &ierr);
    gmshModelOccAddPoint(x+11.5, y+5, z, 1, 20, &ierr);
    gmshModelOccAddPoint(x+14, y+5, z, 1, 21, &ierr);

    gmshModelOccAddPoint(x + 1.4, y + 5.5, z, 1, 22, &ierr);
    gmshModelOccAddPoint(x + 7, y + 5.5, z, 1, 23, &ierr);
    gmshModelOccAddPoint(x + 4.24, y + 7.7, z, 1, 24, &ierr);

    gmshModelOccAddPoint(x + 13.6, y + 5.5, z, 1, 25, &ierr);
    gmshModelOccAddPoint(x + 8, y + 5.5, z, 1, 26, &ierr);
    gmshModelOccAddPoint(x + 10.76, y + 7.7, z, 1, 27, &ierr);

    gmshModelOccAddPoint(x + 7.5, y + 10.31, z, 1, 28, &ierr);
    gmshModelOccAddPoint(x + 5.85, y + 9, z, 1, 29, &ierr);
    gmshModelOccAddPoint(x + 9.15, y + 9, z, 1, 30, &ierr);

    gmshModelOccAddPoint(x + 5.15, y + 8.45, z, 1, 31, &ierr);
    gmshModelOccAddPoint(x + 4.65, y + 8.05, z, 1, 32, &ierr);
    gmshModelOccAddPoint(x + 7.5, y + 5.7, z, 1, 33, &ierr);
    gmshModelOccAddPoint(x + 10.35, y + 8.05, z, 1, 34, &ierr);
    gmshModelOccAddPoint(x + 9.85, y + 8.45, z, 1, 35, &ierr);
    

    

    

    gmshModelOccAddLine(1, 2, 1, &ierr);
    gmshModelOccAddLine(2, 3, 2, &ierr);
    gmshModelOccAddLine(3, 4, 3, &ierr);
    gmshModelOccAddLine(4, 5, 4, &ierr);
    gmshModelOccAddLine(5, 6, 5, &ierr);
    gmshModelOccAddLine(6, 7, 6, &ierr);
    gmshModelOccAddLine(7, 1, 7, &ierr);

    gmshModelOccAddLine(8, 9, 8, &ierr);
    gmshModelOccAddLine(9, 10, 9, &ierr);
    gmshModelOccAddLine(10, 11, 10, &ierr);
    gmshModelOccAddLine(11, 12, 11, &ierr);
    gmshModelOccAddLine(12, 13, 12, &ierr);
    gmshModelOccAddLine(13, 14, 13, &ierr);
    gmshModelOccAddLine(14, 15, 14, &ierr);
    gmshModelOccAddLine(15, 8, 15, &ierr);

    gmshModelOccAddLine(16, 17, 16, &ierr);
    gmshModelOccAddLine(17, 18, 17, &ierr);
    gmshModelOccAddLine(18, 16, 18, &ierr);

    gmshModelOccAddLine(19, 20, 19, &ierr);
    gmshModelOccAddLine(20, 21, 20, &ierr);
    gmshModelOccAddLine(21, 19, 21, &ierr);

    gmshModelOccAddLine(22, 23, 22, &ierr);
    gmshModelOccAddLine(23, 24, 23, &ierr);
    gmshModelOccAddLine(24, 22, 24, &ierr);

    gmshModelOccAddLine(25, 26, 25, &ierr);
    gmshModelOccAddLine(26, 27, 26, &ierr);
    gmshModelOccAddLine(27, 25, 27, &ierr);

    gmshModelOccAddLine(28, 29, 28, &ierr);
    gmshModelOccAddLine(29, 30, 29, &ierr);
    gmshModelOccAddLine(30, 28, 30, &ierr);

    gmshModelOccAddLine(31, 32, 31, &ierr);
    gmshModelOccAddLine(32, 33, 32, &ierr);
    gmshModelOccAddLine(33, 34, 33, &ierr);
    gmshModelOccAddLine(34, 35, 34, &ierr);
    gmshModelOccAddLine(35, 31, 35, &ierr);





    

    





    
    

    

    // Ajouter d'autres lignes selon votre géométrie

    // Construction de la boucle de la surface
    int curveTags[] = {1,2,3,4,5,6,7}; // Utilisation de deux lignes pour former une boucle
    gmshModelOccAddCurveLoop(curveTags, 7, 1, &ierr);

    // Création de la surface
    int wireTags[] = {1}; // Utilisation de la boucle comme fil
    int outersurface = gmshModelOccAddPlaneSurface(wireTags, 1, 1, &ierr); 

    
    int innerCurveTags[] = {8,9,10,11,12,13,14,15}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags, 8, 2, &ierr);

    int innerCurveTags2[] = {16,17,18}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags2, 3, 3, &ierr);

    int innerCurveTags3[] = {19,20,21}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags3, 3, 4, &ierr);

    int innerCurveTags4[] = {22,23,24}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags4, 3, 5, &ierr);

    int innerCurveTags5[] = {25,26,27}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags5, 3, 6, &ierr);

    int innerCurveTags6[] = {28,29,30}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags6, 3, 7, &ierr);

    int innerCurveTags7[] = {31,32,33,34,35}; // Les tags des points intérieurs
    gmshModelOccAddCurveLoop(innerCurveTags7, 5, 8, &ierr);

    
    // Créer une nouvelle surface à l'intérieur
    int innerWireTags[] = {2}; // Utilisation de la nouvelle boucle comme fil
    int innersurface = gmshModelOccAddPlaneSurface(innerWireTags, 1, 2, &ierr);

    int innerWireTags2[] = {3}; // Utilisation de la nouvelle boucle comme fil
    int innersurface2 = gmshModelOccAddPlaneSurface(innerWireTags2, 1, 3, &ierr);

    int innerWireTags3[] = {4}; // Utilisation de la nouvelle boucle comme fil
    int innersurface3 = gmshModelOccAddPlaneSurface(innerWireTags3, 1, 4, &ierr);

    int innerWireTags4[] = {5}; // Utilisation de la nouvelle boucle comme fil
    int innersurface4 = gmshModelOccAddPlaneSurface(innerWireTags4, 1, 5, &ierr);

    int innerWireTags5[] = {6}; // Utilisation de la nouvelle boucle comme fil
    int innersurface5 = gmshModelOccAddPlaneSurface(innerWireTags5, 1, 6, &ierr);

    int innerWireTags6[] = {7}; // Utilisation de la nouvelle boucle comme fil
    int innersurface6 = gmshModelOccAddPlaneSurface(innerWireTags6, 1, 7, &ierr);

    int innerWireTags7[] = {8}; // Utilisation de la nouvelle boucle comme fil
    int innersurface7 = gmshModelOccAddPlaneSurface(innerWireTags7, 1, 8, &ierr);

    
   

    

    

    // Couper la surface extérieure par la surface intérieure
    int outer[] = {2,outersurface};
    int inner[] = {2,innersurface};
    int inner2[] = {2,innersurface2};
    int inner3[] = {2,innersurface3};
    int inner4[] = {2,innersurface4};
    int inner5[] = {2,innersurface5};
    int inner6[] = {2,innersurface6};
    int inner7[] = {2,innersurface7};


    gmshModelOccCut(outer,2,inner,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(outer,2,inner2,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(outer,2,inner3,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(outer,2,inner4,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);    
    gmshModelOccCut(outer,2,inner5,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(outer,2,inner6,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    gmshModelOccCut(outer,2,inner7,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);

    


    
    

    // Définition de la taille de maillage
    geoSetSizeCallback(geoSize);
    gmshModelOccSynchronize(&ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);

    // Génération du maillage
    gmshModelMeshGenerate(2, &ierr);
    
}