/*
    Authors :
    KAOUCH Abdelssamad
    AIT KHELIFA Tanina
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpfr.h>

/*
    Overwrites v with a random double value in range min ___ max
*/
void random_mpfr_in(mpfr_t v, int min, int max, int d){
    if(d){
        //Pick a random double between 0 ___ 1
        mpfr_set_d(v, (double) rand()/ RAND_MAX, MPFR_RNDN);

        //From range 0 ___ 1, is now between 0 ____ max - min
        double mm = max-min;
        mpfr_mul_d(v, v, mm, MPFR_RNDN);
        
        //Pick is now between min ___ max
        mpfr_add_d(v, v, min, MPFR_RNDN);
    } else {
        double pick = (rand() % (max - min + 1)) + min;
        mpfr_set_d(v, pick, MPFR_RNDN);
    }
    
}

/*
    Initialize each mpfr_t of a matrix to chosen precision
*/
void init_mpfr_matrix(mpfr_t* m, int n, mpfr_prec_t prec){
    for(int i = 0; i < n*n; i++){
            mpfr_init2(m[i], prec);
    }
    //printf("Out init\n");
}

/*
    Desallocates each mpfr_t of a matrix
*/
void clear_mpfr_matrix(mpfr_t* m, int n){
    for(int i = 0; i < n*n; i++){
        mpfr_clear(m[i]);
    }
    //printf("Out clear\n");
}

/*
    Copy the content of a matrix to an other one.
    After execution, to is overwritten with from.
*/
void copyMatrix(mpfr_t* from, mpfr_t* to, int l, int c){
    for(int i = 0; i<l*c; i++){
        mpfr_set(to[i], from[i], MPFR_RNDN);
    }
}

/*
    Initializes m, a size n*n matrix, with random double values in range min ___ max
*/
void randMatrix(mpfr_t* m, int n, double min, double max, mpfr_prec_t prec){
    for(int i = 0; i < n*n; i++){
        random_mpfr_in(m[i], min, max, prec);
    }
}

/*
    Function used to print out matrices of size lines*columns with 4 decimals places. 
*/
void printMatrix(mpfr_t* m, int lines, int columns){
    for(int line = 0; line<lines; line++){
        for(int col = 0; col<columns; col++){
            printf("|");
            mpfr_out_str(stdout, 10, 4, m[line*columns + col], MPFR_RNDN);
            printf("|");
        }
        printf("\n");
    }
    printf("\n");
}

/*
    Function used to compute the transpose of a lines*columns matrix.
    After execution, to is overwritten with the transpose matrix of from. 
*/
void transpose(mpfr_t* from, mpfr_t* to, int lines, int columns){
    for(int i = 0; i < lines; i++){
        for(int j = 0; j < columns; j ++){
            mpfr_set(to[j*lines + i], from[i*columns + j], MPFR_RNDN);
        }   
    }
}

/* 
    Computes the naive general matrix product. a is (m*n), b is (n*p) c is (m*p).
    After execution :
    - c is overwritten with the matrix product ab
*/
void matrix_product(mpfr_t* a, mpfr_t* b, mpfr_t*c, int m, int n, int p, mpfr_t t1, mpfr_t t2){
    for(int i = 0; i < m; i++){
        for (int j = 0; j < p; j++){
            mpfr_set_d(t1, 0., MPFR_RNDN); // t1 = sum of the products
            for(int k = 0; k < n; k++){
                mpfr_mul(t2, a[i*n + k], b[k*p + j], MPFR_RNDN); // t2 = element of A * element of B
                mpfr_add(t1, t1, t2, MPFR_RNDN); // t1 += t2
            }
            mpfr_set(c[i*p + j], t1, MPFR_RNDN);
        }
    }
}

/*
    Initialize the ID matrix l*c
*/
void initId(mpfr_t* m, int l, int c){
    for(int i = 0; i < l; i++){
        for(int j = 0; j < c; j++){
            if(i == j){
                mpfr_set_d(m[i*c + j], 1., MPFR_RNDN);
            } else {
                mpfr_set_d(m[i*c + j], 0., MPFR_RNDN);
            }
        }
    }
}

/*
    Overwrites g with the Givens' matrix used to annihilate the coefficient (i,j) in the
    l*c matrix m.
*/
void givens_matrix(mpfr_t* m, mpfr_t* g, int l, int c, int i, int j, mpfr_t t1, mpfr_t t2){
    initId(g, l, c);

    mpfr_sqr(t1, m[j*c + j], MPFR_RNDN); 
    mpfr_sqr(t2, m[i*c + j], MPFR_RNDN);
    mpfr_add(t1, t1, t2, MPFR_RNDN);
    mpfr_sqrt(t1, t1, MPFR_RNDN);

    mpfr_div(g[i*c + i], m[j*c + j], t1, MPFR_RNDN);
    mpfr_div(g[j*c + j], m[j*c + j], t1, MPFR_RNDN);
    mpfr_div(g[j*c + i], m[i*c + j], t1, MPFR_RNDN);
    mpfr_div(g[i*c + j], m[i*c + j], t1, MPFR_RNDN);
    mpfr_neg(g[i*c + j], g[i*c + j], MPFR_RNDN);
}

/*
    Computes the QR decompostion of a n*n matrix m using Givens' algorithm
*/
void givens(mpfr_t* m, mpfr_t* r, mpfr_t* q, int n, mpfr_t* q_tmp, mpfr_t* tmp, mpfr_t t1, mpfr_t t2){
    initId(q, n, n);
    copyMatrix(m, r, n, n);

    for(int j = 0; j < n; j++){
        for(int i = j+1; i < n; i++){
            if(mpfr_regular_p(r[i*n + j])){ //If is not 0, inf, or NaN
                givens_matrix(r, q_tmp, n, n, i, j, t1, t2);  //q_tmp = Gij
                matrix_product(q_tmp, r, tmp, n, n, n, t1, t2); //tmp = GijR
                copyMatrix(tmp, r, n, n); //r = GijR

                transpose(q_tmp, tmp, n, n); //tmp = Gij*
                matrix_product(q, tmp, q_tmp, n, n, n, t1, t2); //q_tmp = QGij*
                copyMatrix(q_tmp, q, n, n); // q = QGij*
            }
        }
    }
}

/*
    Overwrites q with the roation matrix used to annihilate the coefficient (i,j) in the
    l*c matrix m and so that rmr* also has a 0 in position (i,j).
*/
void rotation_matrix(mpfr_t* m, mpfr_t* q, int l, int c, int i, int j, mpfr_t t1, mpfr_t t2){
    initId(q, l, c);
    mpfr_sqr(t1, m[i*c + j], MPFR_RNDN); 
    mpfr_sqr(t2, m[(i-1)*c + j], MPFR_RNDN);
    mpfr_add(t1, t1, t2, MPFR_RNDN);
    mpfr_sqrt(t1, t1, MPFR_RNDN);

    mpfr_div(q[i*c + i], m[(i-1)*c + j], t1, MPFR_RNDN);
    mpfr_div(q[(i-1)*c + i-1], m[(i-1)*c + j], t1, MPFR_RNDN);
    mpfr_div(q[(i-1)*c + i], m[i*c + j], t1, MPFR_RNDN);
    mpfr_div(q[i*c + i-1], m[i*c + j], t1, MPFR_RNDN);
    mpfr_neg(q[i*c + i-1], q[i*c + i-1], MPFR_RNDN);
}

/*
    Algorithm that makes de n*n matrix m upper hessenberg
    After execution :
    - r is overwritten with an upper hessenberg matrix
    - q is overwritten with the product of all the rotation matrix used
    Computing qr should give us the matrix m we had at the start
*/
void algo_hessenberg(mpfr_t* m, mpfr_t* r, mpfr_t* q, int n, mpfr_t* q_tmp, mpfr_t* tmp, mpfr_t t1, mpfr_t t2){
    copyMatrix(m, r, n, n);
    for(int j = 0; j < n-2; j++){
        for(int i = n-1; i > j + 1; i--){
            mpfr_set_d(t1, 0., MPFR_RNDN);
            //If the coefficient to annihilate isnt 0 already
            if(!mpfr_equal_p(m[i*n + j], t1)){
                rotation_matrix(r, q_tmp, n, n, i, j, t1, t2);
                
                //QR
                matrix_product(q_tmp, r, tmp, n, n, n, t1, t2);
                copyMatrix(tmp, r, n, n);

                //QRQ*
                transpose(q_tmp, tmp, n, n);
                matrix_product(r, tmp, q_tmp, n, n, n, t1, t2);
                copyMatrix(q_tmp, r, n, n);
            }

        }
    }
}

int cant_end(mpfr_t* a, int n, double epsilon, mpfr_t t){
    for(int i = 1; i < n-1; i++){
        mpfr_abs(t, a[i*n + (i - 1)], MPFR_RNDN);
        //if abs(element) >= epsilon
        if(mpfr_cmp_d(t, epsilon) >= 0){
            return 1;
        }
    }
    return 0;
}

void eigenvalues(mpfr_t* a, mpfr_t* q, mpfr_t* qt, mpfr_t* r, mpfr_t* tmp, int n, int max_it, double epsilon, mpfr_t t1, mpfr_t t2){
    int it = 0;
    algo_hessenberg(a, r, q, n, qt, tmp, t1, t2);
    copyMatrix(r, a, n, n);
    // printf("Matrice a hessenberg :\n");
    // printMatrix(a, n, n);

    time_t start = time(NULL);
    while((cant_end(a, n, epsilon, t1)) && (it < max_it)){
        //printf("%d\n", it);
        givens(a, r, q, n, qt, tmp, t1, t2);
        transpose(q, qt, n, n); //tmp = Gij* //QTMP contains the transpose of Q
        matrix_product(qt, a, r, n, n, n, t1, t2); //R contains Q*A
        matrix_product(r, q, a, n, n, n, t1, t2); //A contains Q*AQ
        it++;
    }
    //printf("Ended after %d iterations.\n", it);
    time_t end = time(NULL);
    printf("Chrono : %lds\n", (end - start));
}

int main(int argc, char** argv){
    //Seeding for tests init
    srand(time(NULL));
    double epsilon;
    int min, max, n, d, prec;
    if(argc == 7){
        d = atoi(argv[1]);;
        n = atoi(argv[2]);
        min = atoi(argv[3]);
        max = atoi(argv[4]);
        epsilon = strtod(argv[5], NULL);
        prec = atoi(argv[6]);
    } else {
        min = 1;
        max = 10;
        epsilon = 0.0001;
        n = 20;
        d = 0;
        prec = 100;
    }


    int max_it = 10000;

    mpfr_t t1;
    mpfr_t t2;
    mpfr_inits2(prec, t1, t2, NULL);
    mpfr_set_d(t2, 1., MPFR_RNDN);

    mpfr_t *m = malloc(n * n * sizeof(*m));
    mpfr_t *r = malloc(n * n * sizeof(*r));
    mpfr_t *q = malloc(n * n * sizeof(*q));
    mpfr_t *q_tmp = malloc(n * n * sizeof(*q_tmp));
    mpfr_t *tmp = malloc(n * n * sizeof(*tmp));

    init_mpfr_matrix(m, n, prec);
    init_mpfr_matrix(r, n, prec);
    init_mpfr_matrix(q, n, prec);
    init_mpfr_matrix(q_tmp, n, prec);
    init_mpfr_matrix(tmp, n, prec);

    randMatrix(m, n, min, max, d);
    printf("Matrice m\n");
    printMatrix(m, n, n);

    eigenvalues(m, q, q_tmp, r, tmp, n, max_it, epsilon, t1, t2);
    printf("Matrice m :\n");
    printMatrix(m, n, n);

    clear_mpfr_matrix(m, n);
    clear_mpfr_matrix(r, n);
    clear_mpfr_matrix(q, n);
    clear_mpfr_matrix(q_tmp, n);
    clear_mpfr_matrix(tmp, n);
    mpfr_clears(t1, t2, NULL);
    mpfr_free_cache();
}