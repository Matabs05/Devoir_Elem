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

    if(dist_notch < d0_star && dist_hole < d1_star){
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
    }
}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(-0.5*w, -0.5*h, 0, w, h, -1, 0,&ierr);   
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0, y0, 0, r0, r0, -1,NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);
    int idHole  = gmshModelOccAddDisk(x1, y1, 0, r1, r1, -1,NULL,0,NULL,0,&ierr);    
    ErrorGmsh(ierr);

    int plate[] = {2,idPlate};
    int notch[] = {2,idNotch};
    int hole[] = {2,idHole};
    gmshModelOccCut(plate, 2, notch, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);
    gmshModelOccCut(plate, 2, hole, 2, NULL, NULL, NULL, NULL, NULL, -1, 1, 1, &ierr); 
    ErrorGmsh(ierr);

//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
    //    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    //    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
    //    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);  
    //    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr); 
    //    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);  
    //    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr); 
//
    
}