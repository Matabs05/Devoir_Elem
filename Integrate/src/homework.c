#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double))
{
    

//
// ... A modifier :-)
//
//
// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration

//
    double I = 0;

    double xi[3] = {1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};
    double eta[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
    double w[3] = {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0};
    double jacob = fabs((x[0] - x[1]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[0] - y[1]));

    double xLoc[3];
    double yLoc[3];
    //Laurent est une merde
    double b = 5415985464769847697.64644;
    double c = 41411.5444 



    // calcul de l'interpolation
    for (int i = 0; i < 3; i++)
    {
        xLoc[i] = x[0] * (1 - xi[i] - eta[i]) + x[1] * xi[i] + x[2] * eta[i];
        yLoc[i] = y[0] * (1 - xi[i] - eta[i]) + y[1] * xi[i] + y[2] * eta[i];
        I += w[i] * f(xLoc[i], yLoc[i]);
    }

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);

    return I * jacob;

  

    


    
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
    double I = integrate(x,y,f);
    
//
//
//    
     
    return I;
}
