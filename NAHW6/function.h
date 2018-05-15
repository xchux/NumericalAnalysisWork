#define EPSILON 0.000001
#define MAXSTEP 1000
#define PI 3.141592653

void gen_mtx(double **A, int n);

void print_mtx(double **A, int n);

void copy_mtx(double **A, double **cA, int n);

void make_identity_mtx(double **A, int n);

void make_rotate_mtx(double **A, int p, int q, double c, double s, int n);

double **alloc_mtx(int n);

void pre_mtx_mult(double **R, double **A, int n);

void post_mtx_mult(double **R, double **A, int n);

void transpose_mtx(double **A, int n);

/*
	傳回找到A矩陣中具最大絕對值的非對角線元素
	並用 p, q 記錄其index 
*/
double max_off_diag_entry(double **A, int *p, int *q, double *offdiag_sum, int n);

double  inner_product(double *a, double *b, int n);

double *alloc_vec(int n);

void mtx_vec_mult(double *a, double **A, double *b, int n);

double vec_norm2(double *a, int n);

/*
	利用旋轉矩陣對A不斷做similiar transform
	最後得 A' = pT * A * P
	A' is diagonal mtx
	eigenvalue is at the main diagonal
	eigenvector = A' * P
*/
int jacobian_method(double **A, double **P, int flag, int n);

//extract engenvalue & eigenvector
void extract_eigen(double **A, double **P, double *v, int n);

//calculate the 2 norm of A * v - lamda * v and return
double cal_norm(double **A, double *vector, double value, int n);

//print eigen vector inner product
void eigen_vector_inner_product(double **A, int n);

//ouput csv file
void create_marks_csv(char *filename, char *item1 ,char *item2,  double *offDiag, int n);
