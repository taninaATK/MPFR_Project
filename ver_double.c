/*
    Authors :
    KAOUCH Abdelssamad
    AIT KHELIFA Tanina
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/*
    Returns a random double value in range min ___ max
*/
double random_double_in(int min, int max, int d){
    if(d){
    //Pick a random double between 0 ___ 1
    double pick = ((double) rand()) / RAND_MAX;

    //From range 0 ___ 1, is now between 0 ____ max - min
    pick = (max - min) * pick;
    
    //Pick is now between min ___ max
    pick = (min + pick);

    //Returns the value of pick
    return pick;
    } else {
        double pick = (rand() % (max - min + 1)) + min;
        return pick;
    }
}

/*
    Initializes m, a size n*n matrix, with random double values in range min ___ max
*/
void randMatrix(double* m, int n, double min, double max, int d){
    for(int i = 0; i < n*n; i++){
        m[i] = random_double_in(min, max, d);
    }
}

/*
    Copy the content of a matrix to an other one.
    After execution, to is overwritten with from.
*/
void copyMatrix(double* from, double* to, int l, int c){
    for(int i = 0; i<l*c; i++){
        to[i] = from[i];
    }
}

/*
    Function used to print out matrices of size lines*columns with 4 decimals places. 
*/
void printMatrix(double* matrix, int lines, int columns){
    for(int line = 0; line<lines; line++){
        for(int col = 0; col<columns; col++){
            printf(" %4f ", matrix[columns*line + col]);
        }
        printf("\n");
    }
    printf("\n");
}

/*
    Function used to compute the transpose of a lines*columns matrix.
    After execution, to is overwritten with the transpose matrix of from. 
*/
void transpose(double* from, double* to, int lines, int columns){
    for(int i = 0; i < lines; i++){
        for(int j = 0; j < columns; j ++){
            to[j*lines + i] = from[i*columns + j];
        }   
    }
}

/* 
    Computes the naive general matrix product. a is (m*n), b is (n*p) c is (m*p).
    After execution :
    - c is overwritten with the matrix product ab
*/
void matrix_product(double* a, double* b, double*c, int m, int n, int p){
    double sum;
    for(int i = 0; i < m; i++){
        for (int j = 0; j < p; j++){
            sum = 0.;
            for(int k = 0; k < n; k++){
                sum += a[i*n + k]*b[k*p + j];
            }
            c[i*p + j] = sum;
        }
    }
}

/*
    Initialize the ID matrix l*c
*/
void initId(double* m, int l, int c){
    for(int i = 0; i < l; i++){
        for(int j = 0; j < c; j++){
            if(i == j){
                m[i*c + j] = 1;
            } else {
                m[i*c + j] = 0;
            }
        }
    }
}

/*
    Overwrites g with the Givens' matrix used to annihilate the coefficient (i,j) in the
    l*c matrix m.
*/
void givens_matrix(double* m, double* g, int l, int c, int i, int j){
    initId(g, l, c);
    double down = sqrt((pow(m[j*c + j], 2) + pow(m[i*c + j], 2)));
    double coef_c = m[j*c + j]/down;
    double coef_s = m[i*c + j]/down;
    g[i*c + i] = coef_c;
    g[j*c + j] = coef_c;
    g[j*c + i] = coef_s;
    g[i*c + j] = -coef_s;
}

/*
    Classical Givens' algorithm on square matrices of sizes n
    (Meaning we are only going to nullify the elements on the subdiagonal of m)
    On return :
    - q is overwritten with the product of all the Givens' matrices used
    - r is overwritten with an upper triangular matrix
*/
void givens(double* m, double* q, double* r,int n){
    double* q_tmp = malloc(n*n*sizeof(double)); //Stocker toutes les givens' G_x...G_0
    double* tmp = malloc(n*n*sizeof(double)); //Stocker le résultats des produits de matrices
    initId(q, n, n);
    copyMatrix(m, r, n, n);

    for(int j = 0; j < n; j++){
        for(int i = j+1; i < n; i++){
            if(r[i*n + j] != 0.){
                givens_matrix(r, q_tmp, n, n, i, j);  //q_tmp = Gij
                matrix_product(q_tmp, r, tmp, n, n, n); //tmp = GijR
                copyMatrix(tmp, r, n, n); //r = GijR
                transpose(q_tmp, tmp, n, n); //tmp = Gij*
                matrix_product(q, tmp, q_tmp, n, n, n); //q_tmp = QGij*
                copyMatrix(q_tmp, q, n, n); // q = QGij*
            }
        }
    }

    free(q_tmp);
    free(tmp);
}

/*
    Overwrites r with the roation matrix used to annihilate the coefficient (i,j) in the
    l*c matrix m and so that rmr* also has a 0 in position (i,j).
*/
void rotation_matrix(double* m, double* r, int l, int c, int i, int j){
    initId(r, l, c);
    double down = sqrt((pow(m[i*c + j], 2) + pow(m[(i-1)*c + j], 2))); //A[i][j] et A[i-1][j]
    double coef_c = m[(i-1)*c + j]/down;
    double coef_s = m[i*c + j]/down;
    
    r[i*c + i] = coef_c;
    r[(i-1)*c + i-1] = coef_c;
    r[i*c + i-1] = 0.-coef_s;
    r[(i-1)*c + i] = coef_s;

}

/*
    Algorithm that makes de n*n matrix m upper hessenberg
    After execution :
    - r is overwritten with an upper hessenberg matrix
    - q is overwritten with the product of all the rotation matrices used
    Computing (r)*m should give us the matrix m we had at the start
*/
void algo_hessenberg(double* m, double* q, double* r, int n){
    double* q_tmp = malloc(n*n*sizeof(double)); //Stocker toutes les rotation R_x...R_0
    double* tmp = malloc(n*n*sizeof(double)); //Stocker le résultats des produits de matrices
    copyMatrix(m, r, n, n);
    for(int j = 0; j < n-2; j++){
        for(int i = n-1; i > j + 1; i--){
            //If the coefficient to annihilate isnt 0 already
            if(r[i*n + j] != 0.){
                rotation_matrix(r, q_tmp, n, n, i, j);

                //QR
                matrix_product(q_tmp, r, tmp, n, n, n);
                copyMatrix(tmp, r, n, n);

                //QRQ*
                transpose(q_tmp, tmp, n, n);
                matrix_product(r, tmp, q_tmp, n, n, n);
                copyMatrix(q_tmp, r, n, n);
                //printf("Matrice q :\n");
                //printMatrix(q, n, n);
            }

        }
    }
    free(q_tmp);
    free(tmp);
}

/*
    Returns non zero when NOT all the coefficient of the first sudiagonal of the n*n
    matrix a are strictly smaller than epsilon.
*/
int cant_end(double* a, int n, double epsilon){
    for(int i = 1; i < n-1; i++){
        if(fabs(a[i*n + (i - 1)]) >= epsilon){
            return 1;
        }
    }
    return 0;
}

void eigenvalues(double* a, double* q, double* qt, double* r, int n, int max_it, double epsilon){
    int it = 0;
    algo_hessenberg(a, q, r, n);
    copyMatrix(r, a, n, n);
    printf("Matrice m hessenberg :\n");
    printMatrix(r, n, n);
    time_t start = time(NULL);
    while((cant_end(a, n, epsilon)) && (it < max_it)){
        givens(a, q, r, n);
        transpose(q, qt, n, n); //QT contains the transpose of Q
        matrix_product(qt, a, r, n, n, n); //A = Q*A
        matrix_product(r, q, a, n, n, n); //A contains the product Q*AQ
        it++;
    }
    //printf("Ended after %d iterations.\n", it);
    time_t end = time(NULL);
    printf("Chrono : %ld\n", (end - start));
}

/*
    Fonction qui initialise une matrice 4*4, utilisée pour du debug 
*/
void _initTest(double* m){
    //Init de m (poly page 21)
    m[0] = 1.;
    m[1] = 2.;
    m[2] = 3.;
    m[3] = 4.;
    m[4] = 5.;
    m[5] = 6.;
    m[6] = 7.;
    m[7] = 8.;
    m[8] = 9.;
    m[9] = 1.;
    m[10] = 2.;
    m[11] = 3.;
    m[12] = 4.;
    m[13] = 5.;
    m[14] = 6.;
    m[15] = 7.;
}

void checkTest(){
}

int main(int argc, char** argv){
    //Seeding for inits
    srand(time(NULL));
    double epsilon;
    int min, max, n, d;
    if(argc == 6){
        d = atoi(argv[1]);;
        n = atoi(argv[2]);
        min = atoi(argv[3]);
        max = atoi(argv[4]);
        epsilon = strtod(argv[5], NULL);
    } else {
        min = 1;
        max = 10;
        epsilon = 0.0001;
        n = 20;
        d = 0;
    }

    

    int max_it = 10000;


    double* m = malloc(n * n * sizeof(double));
    double* r = malloc(n * n * sizeof(double));
    double* q = malloc(n * n * sizeof(double));
    double* qt = malloc(n * n * sizeof(double));

    //_initTest(m);
    randMatrix(m, n, min, max, d);
    printf("Matrice m :\n");
    printMatrix(m, n, n);
    eigenvalues(m, q, qt, r, n, max_it, epsilon);
    printf("Matrice m :\n");
    printMatrix(m, n, n);
    free(m);
    free(r);
    free(qt);
    free(q);
}