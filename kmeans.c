#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h> 
#include <ctype.h>

void validity(int k, int n, char *kPlain, int max_iter, int numOfArgs); 
void readInput(char* fileNam, double** vectors);
void intializeCentroids(double**centroids, double**vectors, int k, int dim);
double Distance(double* x, double* y, int dim);
void updateCentroid(double** centroids, double** clusters, int i, int dim, int numOfVec, double**vectors);
void assignVectorToCluster(double** vectors,double** centroids,double** clusters, int vectorIndex, int dim, int numOfVec, int k);
void writeToOutput(char* fileName,double** centroids);
void freeMatrix(double** matrix, int rowNum);
double** matrixAlloc(int rowNum, int colNum , int stage);
int calcNumOfVec(char *filename);
int calcDim(char *filename);
void writeToFile(char* outFileName, int k, int dim, double** centroids);
void printInvalidInput(int stage);
void printAnErrorHasOccured(int stage);
void freeAll();
void freeVectorsCentroidsAndClusters();
void freeVectorsAndCentroids();
void freeFunc(int stage);
void isNumber(char *s);

double kDouble, max_iter_double;
int g;
int row;
int col;
FILE * outFile;

char *input, *output;
int num_of_unchanged_centorids, i, j, r, l, k, max_iter,t;
double **vectors; 
int numOfArgs;
double **centroids;
double** clusters;
double *old_centroid, *new_centroid;
int numOfVec;
int dim;
char *kPlain;
char *max_iter_plain;






/*For successful running return 0 otherwise 1*/
int main(int argc, char* argv[]) {
    numOfArgs = argc;
    kPlain = argv[1];
    kDouble = atof(argv[1]);
    if (fmod(kDouble,1) != 0.0){
        printInvalidInput(0);
    }
    k = atoi(argv[1]);
    if (numOfArgs == 4){
        max_iter = 200;
        input = argv[2];
        output = argv[3];
        }
    else if (numOfArgs == 5){
        max_iter_plain = argv[2];
        max_iter_double = atof(argv[2]);
        if (fmod(max_iter_double,1) != 0.0){
            printInvalidInput(0);
        }
        max_iter = atoi(argv[2]);
        input = argv[3];
        output = argv[4];
        }
    else {
        printInvalidInput(0);
    }
    
    

    numOfVec = calcNumOfVec(input);
    dim = calcDim(input);

    validity(k, numOfVec, kPlain, max_iter, numOfArgs);

    vectors = matrixAlloc(numOfVec, dim, 0);
    centroids = matrixAlloc(k, dim, 1);
    if ((vectors == NULL) || (centroids == NULL)){
        printAnErrorHasOccured(2);
    }

    readInput(input, vectors);

    
    intializeCentroids(centroids, vectors, k, dim);

    for (i=0; i<max_iter; ++i){
        clusters = matrixAlloc(k, numOfVec,2);
        if (clusters == NULL){
            printAnErrorHasOccured(3);
        }
        num_of_unchanged_centorids = 0;
        for (j=0; j<numOfVec; ++j){
            assignVectorToCluster(vectors, centroids, clusters, j, dim, numOfVec, k);
            }
        for (l=0; l<k; ++l){
            old_centroid = malloc(dim * sizeof(double));
            if (old_centroid == NULL){
                printAnErrorHasOccured(4);
            }
            for (g=0; g<dim; g++)
            {
                old_centroid[g] = centroids[l][g];
            }
            updateCentroid(centroids, clusters, l, dim, numOfVec, vectors);
            if (sqrt(Distance(old_centroid, centroids[l], dim))<0.001){
                num_of_unchanged_centorids += 1;
                }
            }
            free(old_centroid);
        freeMatrix(clusters, k);
        if (num_of_unchanged_centorids == k){
            break;
        }
    }

    /*write results to file*/
    writeToFile(output, k, dim, centroids);

    freeVectorsAndCentroids();
    
    exit(0);      
}

void freeAll(){
    freeMatrix(vectors, numOfVec);
    freeMatrix(centroids, k);
    freeMatrix(clusters, k);
    free(old_centroid);
}

void freeVectorsCentroidsAndClusters(){
    freeMatrix(vectors, numOfVec);
    freeMatrix(centroids, k);
    freeMatrix(clusters, k);
}

void freeVectorsAndCentroids(){
    freeMatrix(vectors, numOfVec);
    freeMatrix(centroids, k);
}

void freeFunc(int stage){
    if (stage == 1){
        freeMatrix(vectors, numOfVec);
    }
    else if (stage == 2){
        freeVectorsAndCentroids();
    }
    else if (stage == 3){
        freeVectorsCentroidsAndClusters();
    }
    else if (stage == 4){
        freeAll();
    }
}

void printInvalidInput(int stage) {
    freeFunc(stage);
    printf("Invalid Input!\n");
    exit(1);
}

void printAnErrorHasOccured(int stage) {
    freeFunc(stage);
    printf("An Error Has Occurred\n");
    exit(1);
}

void isNumber(char *s){
    t = 0;
    while(s[t] != '\0'){
        if (!isdigit(s[t])){
            printInvalidInput(0);
        }
    t++;
    }
}
/* checking if the size of k is valid, if valid return 1 else 0*/
void validity(int k, int n, char *kPlain, int max_iter, int numOfArgs) {
    /* check of k*/ 
    isNumber(kPlain);
    if (atoi(kPlain) <= 1 ) {
        printInvalidInput(0);
    }
    if ((k<=1) || (k>=n)){
        printInvalidInput(0);
    }
    /* check of N (numofvec)*/ 
    if (n <= 0){
        printInvalidInput(0);
    }
    /* check of max_iter*/ 

    if (numOfArgs == 5){
        isNumber(max_iter_plain);
        if (max_iter <= 0) {
            printInvalidInput(0);
        }
    }
}


int calcNumOfVec(char *filename){
    char ch, ifeof;
    FILE* f = fopen(filename, "r");
    int numOfVec = 0;
    do {
        ch = fgetc(f);
        if(ch == '\n')
            numOfVec++;
    } while (ch != EOF);
    ifeof = fclose(f);
    if (ifeof == EOF) {
        printAnErrorHasOccured(0);
    }
    return numOfVec;
}

int calcDim(char *filename) {
    char ch, ifeof;
    FILE* f = fopen(filename, "r");
    int commasCount = 0;
    do {
        ch = fgetc(f);
        if(ch == ',')
            commasCount++;
    } while (ch != '\n');
    ifeof = fclose(f);
    if (ifeof == EOF) {
        printAnErrorHasOccured(0);
    }
    return (commasCount + 1);
}


/* func to read data from file, return a list of lstsc(vectors = [],vectorCooirdinates = [])*/
void readInput(char *fileName, double **vectors)
{
    double coordinate;
    char comma;
    FILE *fptr;
    int row = 0, col = 0;
    vectors[0][0]=0;
    if (!(fptr = fopen(fileName, "r"))) {
        printInvalidInput(2);  
    }
    while (fscanf(fptr, "%lf%c", &coordinate, &comma) == 2) {
        vectors[row][(col)++] = coordinate;
        if (comma == '\n' || comma == '\r') {
            /*fscanf(fptr, "%c", &comma);*/
            row ++;
            col = 0;
        }
        else if (comma != ',') {
            printAnErrorHasOccured(2); 
        }
    }
    fclose(fptr);
}

/* creates lst in length k, initialize with first k vectors */
void intializeCentroids(double**centroids, double**vectors,int k, int dim)
{
    int i;
    for (i=0; i<k; i++){
        int j;
        for (j=0; j<dim; j++)
            centroids[i][j] = vectors[i][j];
        }
}

/* calculates euclideanDistance between 2 vectors, returns dis */
double Distance(double* x, double* y,int dim)
{
    double sumOfSquaredDiffrence = 0;
    int i;
    for (i=0; i<dim; ++i){
        sumOfSquaredDiffrence += pow(x[i]-y[i],2);}
    return sumOfSquaredDiffrence;
}

/* updates centroid */
void updateCentroid(double** centroids, double** clusters, int i, int dim, int numOfVec, double**vectors){
    int j;
    for (j=0 ; j<dim;++j){
        int cluster_size = 0;
        double sum = 0.0;
        int r=0;
        int vector_ind;
        while((clusters[i][r] != -1.0) && (r < numOfVec)){
            vector_ind = clusters[i][r];
            sum += vectors[vector_ind][j];
            cluster_size+=1;
            r++;
        }
        if (cluster_size == 0){
            printAnErrorHasOccured(4);
        }
        centroids[i][j] = sum/cluster_size;
    }
}

/* assigns vector to closest cluster */
void assignVectorToCluster(double** vectors,double** centroids,double** clusters, int vectorIndex, int dim, int numOfVec, int k){
    double min = DBL_MAX;
    int index = 0;
    int i;
    int r;
    double distance;
    for (i=0; i<k;++i){
        distance = Distance(vectors[vectorIndex], centroids[i], dim);
        if (distance < min){
            min = distance;
            index = i;}}
    for (r=0; r<numOfVec; ++r){
        if (clusters[index][r] == -1.0){
            clusters[index][r] = vectorIndex;
            break;
        }
    }

}

/* Memory allocation for 2D doubles matrix */
double** matrixAlloc(int rowNum, int colNum, int stage) {
    double **p = (double **)malloc(rowNum * sizeof(double*));
    if (p == NULL) {
        printAnErrorHasOccured(stage);
    }
    for (i=0; i < rowNum; i++) {
        p[i] = (double*) malloc(colNum * sizeof(double));
        if (p[i] == NULL){
            printAnErrorHasOccured(stage);
        }
        for (j=0; j < colNum; j++) {
            p[i][j] = -1.0;
        }

    }
    return p;
}


/* freeing the memory of a given 2D matrix */
void freeMatrix(double** matrix, int rowNum) {
    for (i=0; i < rowNum; i++) {
        free(matrix[i]);
    }
    free(matrix);
}


void writeToFile(char* outFileName, int k, int dim, double** centroids){
    outFile = fopen(outFileName, "w");
    for (row = 0; row < k; row++) {
        for (col = 0; col < dim; col++) {
            if (col != (dim - 1)) {
                fprintf(outFile, "%.4f,", centroids[row][col]);
            } else {
                fprintf(outFile, "%.4f\n", centroids[row][col]);
            }
        }
    }

    if (fclose(outFile) == EOF) {
        printAnErrorHasOccured(2);
    }
}
