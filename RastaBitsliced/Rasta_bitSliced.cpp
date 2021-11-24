/* Code implementing bit-sliced version of Rasta for timing purposes.  Block size might be varied, and constant and key additions are not implemented.  */

#include <stdlib.h>
#include <iostream>
#include <helib/helib.h>
#include <helib/matmul.h>
#include <NTL/BasicThreadPool.h>

using namespace std;
using namespace NTL;
using namespace helib;

#include "matrixFromFile.h"

double timeChi=0, timeLinLayer=0;
int BLOCK_SIZE=351;
int ***Mlist;

void fillMlist(int n_rounds){
  /* Reads in n_rounds matrices from files named "matrix_a.txt" (where a is 1,2,...,). */
  int i, j, m;
  FILE *fp;
  char filename[20];

  Mlist=(int ***)malloc(n_rounds*sizeof(int **));
  for(m=0; m<n_rounds; ++m){
    sprintf(filename,"matrix_%d.txt",m+1);
    fp=fopen(filename,"r");
    Mlist[m]=(int **)malloc(BLOCK_SIZE*sizeof(int *));
    for(i=0; i<BLOCK_SIZE; ++i){
      Mlist[m][i]=(int *)malloc(BLOCK_SIZE*sizeof(int));
      for(j=0; j<BLOCK_SIZE; ++j)
	fscanf(fp,"%d",Mlist[m][i]+j);
    }
    fclose(fp);
  }
}

vector<Ctxt> Chi(vector<Ctxt> input){
  /* Returns Chi(input). */
  vector<Ctxt> output;
  int i;

  for(i=0; i<BLOCK_SIZE; ++i){
    Ctxt tmp(input[(i+1)%BLOCK_SIZE]);
    tmp*=input[(i+2)%BLOCK_SIZE];
    tmp+=input[(i+2)%BLOCK_SIZE];
    tmp+=input[i];
    output.push_back(tmp);
  }
  
  return output;
}

vector<Ctxt> matrixMultiplicationPlain(vector<Ctxt>input, int **M, Ctxt ZeroCtxt){
  /* Computes and returns input*M in the straigh-forward way. */
  int i, j;
  vector<Ctxt> output;

  for(i=0; i<BLOCK_SIZE; ++i){
    Ctxt tmp(ZeroCtxt);
    output.push_back(tmp);
  }
  //output initialized

  for(j=0; j<BLOCK_SIZE; ++j){
    for(i=0; i<BLOCK_SIZE; ++i){
      if(M[i][j]!=0)
	output[j]+=input[i];
    }
  }

  return output;
}

vector<Ctxt> matrixMultiplicationM4R(vector<Ctxt> input, int **M, Ctxt ZeroCtxt){
  /* Computes and returns input*M using the method of the four russians.  For HElib, plain matrix-multiplication seems to be optimal. */
  int i, j, t=2, nTables, tableSize, tableIndex, p, rest;
  vector<Ctxt> output;
 
  //Create lookup-tables
  nTables=BLOCK_SIZE/t;//will be 70 when BLOCKSIZE=351 and t=5
  tableSize=1<<t;
  vector<vector<Ctxt>> ciphertextTable;
  // printf("Ready to build tables\n");
  for(i=0; i<nTables; ++i){
     vector<Ctxt> tab;
    for(j=0; j<tableSize; ++j){
      Ctxt tmp(ZeroCtxt);
      tab.push_back(tmp);
    }
    //ciphertextTable[i] initialized
    for(tableIndex=1; tableIndex<tableSize; ++tableIndex){//nothing to add when tableIndex is 0
      for(p=0; p<t; ++p){
	if((1<<p)&tableIndex)//must add ciphertext p of input block i onto current table index
	  tab[tableIndex]+=input[t*i+p];
      }
    }
    ciphertextTable.push_back(tab);
    //ciphertextTable[i] initialized
  }//all look-up tables initialized
  //printf("tables built\n");
  
  for(i=0; i<BLOCK_SIZE; ++i){
    Ctxt tmp(ZeroCtxt);
    output.push_back(tmp);
  }
  //output initialized

  //Do the actual multiplication
  rest=BLOCK_SIZE%t;
  for(j=0; j<BLOCK_SIZE; ++j){//for every column...
    for(i=0; i<nTables; ++i){//...look up every block of t rows
      tableIndex=0;
       for(p=0; p<t; ++p){
	if(M[(t*i)+p][j]!=0)
	  tableIndex^=(1<<p);
      }//tableIndex indicates the element to add to output
      if(tableIndex!=0)
	output[j]+=ciphertextTable[i][tableIndex];
    }
    //adding the last BLOCK_SIZE mod t elements of input
    for(p=BLOCK_SIZE-rest; p<BLOCK_SIZE; ++p){
      if(M[p][j]!=0)
	output[j]+=input[p];
    }//column j fully processed
  }
    
  return output;
}


vector<Ctxt> Round(vector<Ctxt> input, int r, Ctxt ZeroCtxt){
  /* Returns the state after one round of Rasta. */
  vector<Ctxt> output;
  int i;
 
  clock_t chi_start=clock();
  output=Chi(input);
  clock_t chi_end=clock();
  timeChi+=(double)(chi_end - chi_start)/CLOCKS_PER_SEC;
  
  clock_t lin_start=clock();
  //output=matrixMultiplicationM4R(output,Mlist[r],ZeroCtxt);//select which matrix multiplication to use
  output=matrixMultiplicationPlain(output,Mlist[r],ZeroCtxt);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;
  
  return output;
}


int main(int argc, char* argv[])
{
  /*  Example of BGV scheme  */
  int i, j;
  double secs_taken;
  
  // Plaintext prime modulus
  unsigned long p = 2;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m=30133;
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 200;
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 2;

  clock_t setup_begin = clock();

  cout << "Initialising context object..." << endl;
  // Initialize context
  // This object will hold information about the algebra created from the
  // previously set parameters
  Context context = ContextBuilder<BGV>()
    .m(m)
    .p(p)
    .r(r)
    .bits(bits)
    .c(c)
    .build();
    
  // Print the context
  //context.printout();
  //cout << endl;
  
  // Print the security level
  cout << "Security: " << context.securityLevel() << endl;
  
  // Secret key management
  cout << "Creating secret key..." << endl;
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();
  cout << "Generating key-switching matrices..." << endl;
  // Compute key-switching matrices that we need
  addSome1DMatrices(secret_key);

  // Public key management
  // Set the secret key (upcast: SecKey is a subclass of PubKey)
  const PubKey& public_key = secret_key;

  // Get the EncryptedArray of the context
  const EncryptedArray& ea = context.getEA();

  // Get the number of slots (phi(m))
  //long nslots = ea.size(); Bit-sliced, don't care about slots
  //cout << "Number of slots: " << nslots << endl;

  clock_t setup_end = clock();
  secs_taken = (double)(setup_end - setup_begin) / CLOCKS_PER_SEC;
  cout<<"Set-up took "<<secs_taken<<" seconds"<<endl;

  //Reading in random matrices
  clock_t readMat_start=clock();
  fillMlist(7);
  clock_t readMat_end=clock();
  secs_taken = (double)(readMat_end - readMat_start) / CLOCKS_PER_SEC;
  cout<<"Reading in the 7 matrices took "<<secs_taken<<" seconds"<<endl;
  
  // Create ciphertext vector, fill it with encryptions of random bits
  vector<Ctxt> RastaBlock;
  for(i = 0; i < BLOCK_SIZE; ++i){
     Ptxt<BGV> ptxt(context);
     ptxt[0] = (long)random()%2;
     Ctxt ctxt(public_key);
     public_key.Encrypt(ctxt,ptxt);
     RastaBlock.push_back(ctxt);
  }
  printf("Bit capacities in fresh ciphertexts: ");
  for(i=0; i<10; ++i)
    printf("%d  --  ",(int)RastaBlock[i].bitCapacity());
  printf("\n");
    
  Ptxt<BGV> ptxt(context);
  ptxt[0]=0;
  Ctxt ZeroCtxt(public_key);//used for initialization of tables in M4R
  public_key.Encrypt(ZeroCtxt,ptxt);
  
  /* Homomorphic evaluation of Rasta */
  //Initial linear transformation
  clock_t lin_start=clock();
  //RastaBlock=matrixMultiplicationM4R(RastaBlock,Mlist[0],ZeroCtxt);
  RastaBlock=matrixMultiplicationPlain(RastaBlock,Mlist[0],ZeroCtxt);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;
  cout<<"Time used so far: LinLayer "<<timeLinLayer<<" seconds -- Chi "<<timeChi<<" seconds"<<endl;

  for(i=1; i<=6; ++i){
    cout<<"Entering round "<<i<<endl;
    RastaBlock=Round(RastaBlock,i,ZeroCtxt);
    cout<<"Time used so far: LinLayer "<<timeLinLayer<<" seconds -- Chi "<<timeChi<<" seconds"<<endl;
  }
  printf("Bit capacities in ciphertexts after running Rasta: ");
  for(i=0; i<10; ++i)
    printf("%d  --  ",(int)RastaBlock[i].bitCapacity());
  printf("\n");
 
  cout<<"Time used in Chi-transformations: "<<timeChi<<" seconds"<<endl;
  cout<<"Time used in LinLayer: "<<timeLinLayer<<" seconds"<<endl;
}
