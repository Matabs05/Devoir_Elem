

#include"fem.h"


#ifndef NORENUMBER 
double *globalNode;

int compare(const void *a, const void *b){    
    int *i = (int *)a;
    int *j = (int *)b;
    double diff = globalNode[*i] - globalNode[*j]; 
    return (diff < 0) - (diff > 0);
}


void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    int *indice = malloc(theMesh->nodes->nNodes*sizeof(int));
    for (i = 0; i < theMesh->nodes->nNodes; i++) {
        indice[i] = i; }
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;
// 
// A modifier :-)
// debut
//
        case FEM_XNUM : 
            globalNode = theMesh->nodes->X;
            qsort(indice,theMesh->nodes->nNodes,sizeof(int),compare);
            break;
        case FEM_YNUM : 
            globalNode = theMesh->nodes->Y;
            qsort(indice,theMesh->nodes->nNodes,sizeof(int),compare);
            break;            
// 
// end
//

        default : Error("Unexpected renumbering option"); }
    for(i = 0; i <theMesh->nodes->nNodes; i++){
        theMesh->nodes->number[indice[i]] = i; 
    }
    free(indice);
}

#endif
#ifndef NOBAND 

int femMeshComputeBand(femMesh *theMesh)
{
    int myBand = 0;
    int max,min,map[theMesh->nLocalNode],iElem;
    
    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (int j = 0; j < theMesh->nLocalNode; j++){
            map[j] = theMesh->nodes->number[theMesh->elem[iElem*theMesh->nLocalNode+j]];}
        max = map[0]; 
        min = map[0];
        for (int j = 1; j < theMesh->nLocalNode; j++) {
            max = (map[j] > max) ? map[j] : max;
            min = (map[j] < min) ? map[j] : min;
            }
        if (myBand < (max-min)) myBand = max-min; }
        
    return(myBand+1);
}


#endif
#ifndef NOBANDASSEMBLE


void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    // A Ecrire :-)
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        for(j = 0; j < nLoc; j++) {
            if(map[i] <= map[j] && map[j] <= map[i] + myBandSystem->band)
                myBandSystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; } 
    myBandSystem->B[map[i]] += Bloc[i]; }
}


#endif
#ifndef NOBANDELIMINATE


double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    // A completer :-)
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        jend = (k+band < size) ? k+band : size;
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k]; ///XXXX////XXXX
            for (j = i ; j < jend; j++) ///XXXX////XXXX
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution */

    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        jend = (i+band < size) ? i+band : size;
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }


    return(myBand->B);
}


#endif
