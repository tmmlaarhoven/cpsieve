//#################################################################################################
//
// This is a proof-of-concept implementation of the CPSieve algorithm, as originally described 
// in the paper:
//
//			 "Efficient (ideal) lattice sieving using cross-polytope LSH"
//
// This implementation was written by Thijs Laarhoven [mail at thijs dot com] in 2015. The 
// CPSieve algorithm is an extension of the GaussSieve algorithm, introduced by Micciancio
// and Voulgaris at SODA'10, and uses similar techniques to those presented in the HashSieve
// algorithm (see related code elsewhere). The idea of lattice sieving further dates back to 
// work of Ajtai, Kumar, and Sivakumar at STOC'01. 
//
// This code is released under the MIT license (see https://opensource.org/licenses/MIT).
//
//#################################################################################################

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

//########################################
//  		SCHEME PARAMETERS
//########################################
// N:   Dimension of the lattice and the Euclidean space (full-rank lattices)
// SEED: Seed of the input challenge basis (usually 0)
// K: 	Hash length; number of cross-polytopes used for a partition in a single hash table
// T: 	Number of hash tables
// Note that K = 0 and T = 1 corresponds to an unoptimized implementation of the GaussSieve
//########################################

#define N 50								// Lattice dimension
#define SEED 0								// Seed of the intput challenge basis	
#define K 2									// Hash length, AND-composition
#define T 20								// Hash tables, OR-composition
#define BUCKETS 10000 // = (2*N)^k			// Number of buckets in each HyperSimplexSieve hash table (assuming k = 1)

//########################################
//  		GLOBAL PARAMETERS
//########################################
#define MAX_VECTORS 1000000					// Maximum total number of LatticeVectorors in the system
#define MAX_ITERATIONS 100000000000			// Maximum number of "iterations"
#define MAX_COLLISIONS_1 10000000			// Maximum number of type-1 collisions: sampling the 0-vector
#define MAX_COLLISIONS_2 10000000			// Maximum number of type-2 collisions: after reductions, obtaining the 0-vector 
#define COS45 0.7071067811					// Cos(Pi/4) = Cos(45 degrees) = sqrt(2)/2 = 0.707106...
#define PI 3.14159265359

// For each vector we store the entries and its norm squared
struct LatticeVector {
	long long int Crd[N];					// Vector coordinates
	unsigned long long int Norm2;			// Vector squared norm
};

// For the Gram-Schmidt basis we also store the entries and their squared norms
struct RealVector {
	double Crd[N];							// Vector coordinates
	double Norm2;							// Vector squared norm
};

// Each hash table bucket is stored in memory as a "Bucket"
struct Bucket {
	LatticeVector** Pointers;				// List of pointers to list vectors in this bucket
	unsigned long long int Length;			// Number of vectors in the bucket, occupied entries in Pointers
	unsigned long long int Size;			// Internal length of the Pointers-array, usually a power of 2
};

// HashSieve-specific variables
Bucket HashTables[T][BUCKETS];				// Locality-sensitive hash tables

// CPSieve: random rotation matrices for each T and K
RealVector App[T][K][N];

// GaussSieve variables
LatticeVector B[N];							// Lattice basis
RealVector Bs[N];							// Gram-Schmidt basis
double mu[N][N];							// Gram-Schmidt coefficients
LatticeVector Vectors[MAX_VECTORS];			// All lattice vectors in the system
LatticeVector* Stack[MAX_VECTORS];			// The algorithm's stack of pointers to vectors
unsigned long long int VectorsLength = 0;	// Total number of vectors in the system
unsigned long long int StackLength = 0;		// Number of vector pointers on the stack
unsigned long long int Reductions1 = 0;		// Count number of v-reductions
unsigned long long int Reductions2 = 0;		// Count number of w-reductions
unsigned long long int Collisions1 = 0;		// Count number of collisions occurring from sampling
unsigned long long int Collisions2 = 0;		// Count number of list collisions
unsigned long long int Comparisons = 0;		// Count number of inner products between v and w
unsigned long long int Hashes = 0;			// Count number of "hashes" computed
unsigned long long int MinNorm2; 			// Current minimum of squared vector norms
unsigned long long int Target;				// The target length for the current dimension
unsigned long long int Iteration;			// The current iteration count

// Klein sampler variables
RealVector xKlein;
double AKlein;

double randgauss()
{
	double u1 = 1.0 * ((double) rand() / RAND_MAX);
	double u2 = 1.0 * ((double) rand() / RAND_MAX);
	double res = sqrt(-2 * log(u1)) * cos(2 * PI * u2);
	return res;
}

void printvector(LatticeVector* v)
{
	printf("[");
	for(int i = 0; i < 3; i++){
		printf("%8d ", v->Crd[i]);
	}
	printf("...");
	for(int i = N-3; i < N; i++){
		printf(" %8d", v->Crd[i]);
	}
	printf("]\n");
}

void randsphericalvector(RealVector* v)
{
	for(int i = 0; i < N; i++){
		v->Crd[i] = randgauss();
		v->Norm2 += v->Crd[i] * v->Crd[i];
	}
	double n = sqrt(v->Norm2);
	for(int i = 0; i < N; i++){
		v->Crd[i] = v->Crd[i] / n;
	}
	v->Norm2 = 1;
}

void randgaussvector(RealVector* v)
{
	for(int i = 0; i < N; i++){
		v->Crd[i] = randgauss() / sqrt((double)N);
		v->Norm2 += v->Crd[i] * v->Crd[i];
	}
}

// Compute the inner product between different integer (lattice) vectors v and w
long long int ip(LatticeVector* v, LatticeVector* w)
{
	long long int res = 0;
	for(int i = 0; i < N; i++){
		res += v->Crd[i] * w->Crd[i];
	}
	return res;
}

// Compute the inner product between a lattice vector and a real vector
double ip(LatticeVector* v, RealVector* w)
{
	double res = 0;
	for(int i = 0; i < N; i++){
		res += (double)v->Crd[i] * w->Crd[i];
	}
	return res;
}

// Compute the inner product between two real vectors v and w
double ip(RealVector* v, RealVector* w)
{
	double res = 0;
	for(int i = 0; i < N; i++){
		res += v->Crd[i] * w->Crd[i];
	}
	return res;
}

// Add a vector v to a hash table bucket b
void bucketAdd(Bucket* b, LatticeVector* v)
{
	// If the bucket overflows, make a new bucket of twice the size
	if(b->Length == b->Size){
		b->Size <<= 1;
		LatticeVector** NewPointers;
		NewPointers = new LatticeVector*[b->Size];
		for(short i = 0; i < b->Length; i++){
			NewPointers[i] = b->Pointers[i];
		}
		delete [] b->Pointers;
		b->Pointers = NewPointers;
	}
	
	// Insert v into the bucket
	b->Pointers[b->Length] = v;
	b->Length++;
}

// Remove a vector v from a hash table bucket b
void bucketRemove(Bucket* b, LatticeVector* v)
{
	// Find w's position in the hash bucket
	int vPos = 0;
	while(b->Pointers[vPos] != v && vPos < b->Length){
		vPos++;
	}
	if(vPos >= b->Length){
       	perror("Vector not found in bucket...\n");
       	exit(-1);
	}		
	// Make the bucket shorter		
	b->Length--;
	b->Pointers[vPos] = b->Pointers[b->Length];
}

// Compute the locality-sensitive hash of a vector v for table t
int lshashk(LatticeVector* v, int t, int k)
{
	// Apply rotation App[t][k]
	RealVector vrot;
	for(int i = 0; i < N; i++){
		vrot.Crd[i] = 0.;
		for(int j = 0; j < N; j++){
			vrot.Crd[i] += (double) v->Crd[j] * App[t][k][i].Crd[j];
			vrot.Norm2 += vrot.Crd[i] * vrot.Crd[i];
		}
	}
	
	// Use Terasawa and Tanaka's method using orthoplexes
	int maxpos = 0;			// Position of maximum so far
	int tmppos;
	double maxval = vrot.Crd[0] * vrot.Crd[0];
	double tmpval;
	short maxsgn = (vrot.Crd[0] > 0 ? 1 : -1);
		
	// Find the largest entry (including sign)
	for(int i = 1; i < N; i++){
		tmpval = vrot.Crd[i] * vrot.Crd[i];
		if(tmpval > maxval){ // New global maximum
			maxpos = i;
			maxval = tmpval;
			maxsgn = (vrot.Crd[i] > 0 ? 1 : -1);
		}
	}
	
	// Compute hash value
	int res;
	if(maxsgn = 1)
		res = maxpos; 		// between 0 and n - 1		(for positive largest coordinates)
	else
		res = N + maxpos;	// between n and 2n - 1 	(for negative largest coordinates)
	return res;
}

// Compute the locality-sensitive hash of a vector v for table t
int lshash(LatticeVector* v, int t)
{
	int res = 0;
	for(int kk = 0; kk < K; kk++){
		res *= 2 * N;
		res += lshashk(v, t, kk);
	}
	return res;
}

// Add/subtract the vector w to/from v
void add(LatticeVector* v, LatticeVector* w, long long int vw)
{
	if(vw > 0){
		// Subtract w from v
		for(int i = 0; i < N; i++)
			v->Crd[i] -= w->Crd[i];
		v->Norm2 += w->Norm2 - 2 * vw;
	}
	else{
		// Add w to v
		for(int i = 0; i < N; i++)
			v->Crd[i] += w->Crd[i];
		v->Norm2 += w->Norm2 + 2 * vw;
	}
}

// Generate a new, randomly sampled vector v using a naive method
void sampleSimple(LatticeVector* v)
{
	int i, j;
	long long int coeff;
	for (j = 0; j < N; j++){
		v->Crd[j] = 0;
	}
	for (i = 0; i < N; i++){
		coeff = (rand() % 2);
		//coeff = (1.0 * rand() / RAND_MAX > 0.02 ? 0 : 1);
		for (j = 0; j < N; j++){
			v->Crd[j] += (long long int)coeff * (B[i].Crd[j]);
		}
	}
	v->Norm2 = 0;
	for(j = 0; j < N; j++){
		v->Norm2 += v->Crd[j] * v->Crd[j];
	}
}

// Import the basis from the given text file filestring
void importBasis(char* filestring)
{
	int i, j, dgt, busy, crd, crdsgn;
    FILE* input_file;
    
    //char filestring[500];
    //snprintf(filestring, 500, "C:\\Google Drive\\Cpp\\SVP\\dim%usd%u-LLL.txt", N, SEED);
    input_file = fopen(filestring, "r");
    if (input_file == 0){
        perror("Cannot open input file...\n");
        exit(-1);
    }
    else{
    	i = 0;
    	j = 0;
        busy = 0;			// Currently reading a coordinate?
		crd = 0; 			// Coordinate value
		crdsgn = 1; 		// Sign of coordinate
        while ((dgt = fgetc(input_file)) != EOF)
        {
        	if (dgt == '-'){
        		// Start fresh coordinate
        		busy = 1;
        		crd = 0; 
        		crdsgn = -1;
        	}
            else if (isdigit(dgt)){
            	if (busy > 0){
            		// Append digit to coordinate
            		crd *= 10;
            		crd += dgt - '0';
            	}
            	else {
            		// Start fresh coordinate
            		busy = 1;
            		crd = dgt - '0'; 
					crdsgn = 1;           		
            	}
        	}
        	else {
        		if (busy > 0){
            		// Write coordinate to basis
            		B[i].Crd[j] = crd * crdsgn;
            		j++;
            		if (j == N){
						B[i].Norm2 = 0;
						for(int i1 = 0; i1 < N; i1++){
							B[i].Norm2 += B[i].Crd[i1] * B[i].Crd[i1];
						}
						j = 0;
            			i++;
            		}
					busy = 0;
            	}
        	}
		}
	}
    fclose(input_file);
}

// The randRound algorithm as described by Klein
int randRound(double c, double r)
{
	int p = floor(r);
	int q = p + 1;
	double a = r - (double)p;
	double b = 1 - a;
	double spos = 0;
	double sneg = 0;
	for(int i = 0; i < 10; i++){
		spos += exp(-c * (i + b) * (i + b));
		sneg += exp(-c * (i + a) * (i + a));
	}
	double s;
	s = spos + sneg;
	spos = spos / s;
	sneg = sneg / s;
	
	double rr;
	rr = 1.0 * ((double)rand() / RAND_MAX);
	int i = 0;
	if (rr < spos){
		// Integer lies on the positive side of r
		double spos2;
		spos2 = exp(-c * (i + b) * (i + b)) / s;
		while(rr > spos2){
			i++;
			spos2 += exp(-c * (i + b) * (i + b)) / s;
		}
		return q + i;
	}
	else{
		// Integer lies on the negative side of r
		rr = 1 - rr; // Total weight on this side of the curve is 1 - rr
		double sneg2;
		sneg2 = exp(-c * (i + a) * (i + a)) / s;
		while(rr > sneg2){
			i++;
			sneg2 += exp(-c * (i + a) * (i + a)) / s;
		}
		return p - i;
	}
}

// Compute the Gram-Schmidt basis and store it in Bs
void gramSchmidt()
{
	int i,j,k;
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			mu[i][j] = 0;
			Bs[i].Crd[j] = (double) B[i].Crd[j];
		}
		for(k = 0; k < i - 1; k++){
			mu[i][k] = ip(&B[i], &Bs[k]) / Bs[k].Norm2;
			for(j = 0; j < N; j++){
				Bs[i].Crd[j] -= mu[i][k] * Bs[k].Crd[j];
			}
		}
		Bs[i].Norm2 = ip(&Bs[i], &Bs[i]);
	}
}


// Klein's near algorithm -- Call this with d = n - 1 and x = (0, ..., 0)
void nearA(LatticeVector* res, double A, RealVector* x, int d)
{
	if(d == -1){
		for(int i = 0; i < N; i++){
			res->Crd[i] = 0;
		}
		res->Norm2 = ip(res, res);
	}
	else{
		double rd = 0;
		for(int i = 0; i < N; i++){
			rd += (double)x->Crd[i] * Bs[d].Crd[i];
		}
		rd = rd / Bs[d].Norm2;
		//printf("%u\n", Bs[d][N]);
		double cd = A * Bs[d].Norm2;
		double ld = randRound(cd, rd);
		RealVector xp;
		for(int i = 0; i < N; i++){
			xp.Crd[i] = x->Crd[i] + (double)((ld - rd) * Bs[d].Crd[i] - ld * B[d].Crd[i]);
		}
		xp.Norm2 = ip(&xp, &xp);
		nearA(res, A, &xp, d - 1);
		for(int i = 0; i < N; i++){
			res->Crd[i] += (long long int)(ld * B[d].Crd[i]);
		}
		res->Norm2 = ip(res, res);
	}
}

// Initialize Klein sampler; initialize zero-vector x and value A
void initKlein()
{
	for(int i = 0; i < N; i++){
		xKlein.Crd[i] = 0;
	}
	xKlein.Norm2 = 0;
	
	AKlein = log(N) * log(N);
	double minval = Bs[0].Norm2;
	for(int i = 1; i < N; i++){
		if(Bs[i].Norm2 < minval){
			minval = Bs[i].Norm2;
		}
	}
	AKlein /= minval;
	AKlein /= 70.;
}


// Initialize the hash vectors and the hash tables
void initHashes()
{
	// Initialize hash tables as empty
	for(int t = 0; t < T; t++){
		
		// Initialize random rotations
		for(int j = 0; j < N; j++){
			for(int k = 0; k < K; k++){
				randsphericalvector(&App[t][k][j]);
			}
		}

		// Initialize empty hash Buckets
		for(int b = 0; b < BUCKETS; b++){
			HashTables[t][b].Length = 0;
			HashTables[t][b].Size = 4;
			HashTables[t][b].Pointers = new LatticeVector*[4];
		}			
	}
}

// Initialize the stack (and the vectors list) by adding basis vectors to it
void initStack()
{
	// Push all basis vectors to the stack
	for(int i = 0; i < N; i++){
    	for(int j = 0; j < N; j++){
    		Vectors[i].Crd[j] = B[i].Crd[j];
    	}
    	Vectors[i].Norm2 = B[i].Norm2;
    	VectorsLength++;
    	Stack[i] = &Vectors[i];
    	StackLength++;
    }
}

void initParams()
{
	// Hardcoded shortest vector lengths / SVP records according to the SVP challenge database or own experiments
	// If nothing is indicated, these records are for seed 0
	unsigned long long int Targets[200];
	Targets[30] = 2091662;
    Targets[31] = 2117044;
    Targets[32] = 2147531;
    Targets[33] = 2301849;
    Targets[34] = 2302448;
    Targets[35] = 2637604;
    Targets[36] = 2535727;
    Targets[37] = 2470204;
    Targets[38] = 2662328;
    Targets[39] = 3022037;
	Targets[40] = 2898390;
    Targets[41] = 2834609;
    Targets[42] = 2414615;
    Targets[43] = 3037224;
    Targets[44] = 2825373;
    Targets[45] = 3098331;
    Targets[46] = 2989372;
    Targets[47] = 3187572;
    Targets[48] = 2963891; //2964946; //3068302; //3148900; //3222705;
    Targets[49] = 3454355;
    Targets[50] = 3584095;
    Targets[51] = 3551524;
    Targets[52] = 3633605;
    Targets[53] = 3496843;
    Targets[54] = 3694084;
	Targets[55] = 3773021;
	Targets[56] = 3900625;
	Targets[57] = 3815991;
	Targets[58] = 4072324;
	Targets[59] = 3781187;
	Targets[60] = 3779136;
	Targets[61] = 4464769;
	Targets[62] = 4380649;
	Targets[63] = 4228565;
	Targets[64] = 4906284;//4426816;
	Targets[65] = 4396757;
	Targets[66] = 4405628;
	Targets[67] = 4787344;
	Targets[68] = 4588164;
	Targets[69] = 4778537;
	Targets[70] = 4596736;
	Targets[71] = 4963938;
	Targets[72] = 4752400;
	Targets[73] = 4800481;
	Targets[74] = 5085025;
	Targets[75] = 5202961;
	Targets[76] = 5026564;
	Targets[77] = 5500000;
	Targets[78] = 5171076;
	Targets[79] = 5508409;
	Targets[80] = 5166529;
	
	Target = Targets[N];
	StackLength = 0;
	VectorsLength = 0;
	Reductions1 = 0;
	Reductions2 = 0;
	Comparisons = 0;
	Collisions1 = 0;
	Collisions2 = 0;
	Hashes = 0;
	Iteration = 0;
	MinNorm2 = 100000000000;
}

// Sample a new vector using Klein's sampler
void sampleKlein(LatticeVector* res)
{
	nearA(res, AKlein, &xKlein, N - 1);
}


//###############################################################################
//###############################################################################
//###############################################################################

// The main execution
int main(void)
{	
	int bla;
	bla = 0;
	
Start:
	srand(time(0));
	//srand(0);

    char filestring[500];
    snprintf(filestring, 500, "dim%usd%u-LLL.txt", N, SEED);
	
	importBasis(filestring);
	gramSchmidt();
	initKlein();
	initParams();
	initHashes();
	initStack();
    
	printf("=====  HashSieve  ======\n" );
	printf("Dimension (N):  %8d\n", N);
	printf("Hash length (K): %7d\n", K);
	printf("Hash tables (T): %7d\n", T);
	printf("Random seed:    %8d\n", SEED);
	//printf("Probe level:    %8d\n", PROBE);
	//printf("Hypersimplex:   %8d\n", SIMPLEX);
	printf("Target norm^2:  %8d\n", Target);
	printf("------------------------\n");
	
	// Some dummy variables used in the algorithm
	int i,j,k,m,n,t;
	LatticeVector** Candidates;
	LatticeVector* v;
	LatticeVector* w;
	int vHash[T];			// "Temporary" variable for hashes of a target vector v and its rotations
	LatticeVector vrot[N];
	int vHashp;				// Shifted hashes of v used for probing only
	int wHash[T];			// "Temporary" variable for hashes of a candidate vector w
	long long int vw;
	long long int vwAbs;
	long long int vwAbsDouble;
	long long int NCandidates;
	int vReduced;
	time_t start = time(NULL);
	time_t now;
	time_t end;
	int rot;
	
	LatticeVector *shortest;
	
	// The main algorithm loop
	while(Iteration < MAX_ITERATIONS && Collisions2 < MAX_COLLISIONS_2){
		Iteration++;
		
		// Get vector from stack, or sample a new one if the stack is empty
		if(StackLength == 0){
			if(VectorsLength == MAX_VECTORS){
				perror("Vector list overflow...\n");
        		goto End;
			}
			sampleKlein(&Vectors[VectorsLength++]);
			//sampleSimple(&Vectors[VectorsLength++]);
			while(Vectors[VectorsLength-1].Norm2 == 0 && Collisions1 < MAX_COLLISIONS_1){
				Collisions1++;
				sampleKlein(&Vectors[VectorsLength-1]);
				//sampleSimple(&Vectors[VectorsLength-1]);
			}
			v = &Vectors[VectorsLength-1];
		}
		else{
			v = Stack[--StackLength];
		}
		
		vReduced = 0;
		// Check each table for candidate near list vectors
		for(t = 0; t < T; t++){
				
			// Compute v's hash value
			vHash[t] = lshash(v, t);
			Hashes += K;
			Candidates = HashTables[t][vHash[t]].Pointers;
			NCandidates = HashTables[t][vHash[t]].Length;
			
			// Go through the list to find reducing vectors
			for(j = NCandidates - 1; j >= 0; j--){
				w = Candidates[j];
				vw = ip(v, w);				
				Comparisons++;
				vwAbs = (vw > 0 ? vw : -vw);
				vwAbsDouble = (vwAbs << 1);
			
				// Reduce v with w if possible
				if(vwAbsDouble > w->Norm2){
					add(v, w, vw);
					Reductions1++;
					vReduced = 1;
					goto vEnd;
				}
			
				// Reduce w with v if possible
				if(vwAbsDouble > v->Norm2){
					
					// Remove w from the hash tables
					for(int tt = 0; tt < T; tt++){
						wHash[tt] = lshash(w, tt);
						Hashes += K;
						bucketRemove(&HashTables[tt][wHash[tt]], w);
					}
					
					// Reduce w with v
					add(w, v, vw);
					
					Reductions2++;
					if(w->Norm2 > 0)
						Stack[StackLength++] = w;
					else
						Collisions2++;
				}
			}	
		}
		
vEnd:	// We have reached a decision for vector v
		// Push v to stack, list, or delete altogether
		if(vReduced == 0){
			if(v->Norm2 > 0){
				
				// Push v to the hash tables
				for(t = 0; t < T; t++)
					bucketAdd(&HashTables[t][vHash[t]], v);

				// Check for new minimum
				if(v->Norm2 < MinNorm2){
					now = time(NULL);
					printf("New minimum: %11llu (%5d sec)\n", v->Norm2, (now - start)); 
					MinNorm2 = v->Norm2;
				}
				if(v->Norm2 <= Target){
					now = time(NULL);
					printf("Target found: %10llu (%5d sec)\n", v->Norm2, (now - start));
					break;
				}
			}
			else{
				Collisions2++;
			}
		}
		else{
			// Append v to the stack
			Stack[StackLength++] = v;
		}
	}
	
End:
	end = time(NULL);
	printf("------------------------\n");
	
	// Formatting the time taken
	int Time = end - start;
	int Timesec = Time % 60;
	int Timemin = (Time / 60) % 60;
	int Timehr = (Time / 3600) % 24;
	int Timeday = (Time / 86400);
	if(Timeday == 0)					printf("Time:           %02u:%02u:%02u (hh:mm:ss)\n", Timehr, Timemin, Timesec);
	else								printf("Time:    (%3ud) %02u:%02u:%02u (hh:mm:ss)\n", Timeday, Timehr, Timemin, Timesec);
	
	// Formatting the main space costs
	double Space = (T * VectorsLength * sizeof(void*)) + (VectorsLength * N * sizeof(Vectors[0].Crd[0]));
	if(Space < 1000)					printf("Est. space: %12f (bytes)\n", Space);
	else if (Space < 1000000)			printf("Est. space: %12f (kB)\n", Space / 1000);
	else if (Space < 1000000000)		printf("Est. space: %12f (MB)\n", Space / 1000000);
	else if (Space < 1000000000000)		printf("Est. space: %12f (GB)\n", Space / 1000000000);
	else								printf("Est. space: %12f (TB)\n", Space / 1000000000000);
	
	printf("------------------------\n");
	printf("Iterations:   %10llu\n", Iteration);
	printf("Inner products\n");
	printf("- Comparing:%12llu\n", Comparisons);
	printf("- Hashing: %13llu\n", Hashes);
	printf("- Total:   %13llu\n", (Comparisons + Hashes));
	printf("Reductions v: %10llu\n", Reductions1);
	printf("Reductions w: %10llu\n", Reductions2);
	printf("Vectors:      %10llu\n", VectorsLength);
	printf("List length:  %10llu\n", VectorsLength - Collisions2 - StackLength);
	printf("Stack length: %10llu\n", StackLength);
	printf("Collisions\n");
	printf("- Sampling 0: %10llu\n", Collisions1);
	printf("- Reducing:   %10llu\n", Collisions2);
	printf("- Total:      %10llu\n", Collisions1 + Collisions2);
	printf("Shortest:     %10llu\n\n", MinNorm2);	
	
	bla++;
	if(bla < 1){
		goto Start;
	}

}
