#define MAX_NP 1000000
#define NUM_PARAMS 3
#define MAX_ITERATIONS 50

static double jacT[NUM_PARAMS][MAX_NP];
static double jac[MAX_NP][NUM_PARAMS];
static double H[NUM_PARAMS][NUM_PARAMS];
static double d[NUM_PARAMS][NUM_PARAMS];
static double r[MAX_NP];
static double delta[NUM_PARAMS];
static double func[MAX_NP];
static double lambda;
static double ssr;
static int final = 0;

static double lambda_0 = 0.001; //pow(10,-2);		// initial value of lambda

int fitData(double y_dat[], int len, double p[NUM_PARAMS], int print, int plot);
int updateJac(int len, double p[NUM_PARAMS]);
int solveDelta(double y_dat[], int len);
int checkX2(double y_dat[], int len, double p[NUM_PARAMS], int iteration);
int computeErrors(double y_dat[], int len);
double dotProd(double a[], double b[], int len);
