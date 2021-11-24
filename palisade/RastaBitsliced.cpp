/* Bitsliced implementation of Rasta in Palisade/binfhe. 
   Used only for timing purposes, additions of constants and keys not implemented. */
#include <stdlib.h>
#include <time.h>

#include "binfhecontext.h"

using namespace lbcrypto;
using namespace std;

int BLOCK_SIZE=351;
int ***Mlist;
double timeChi=0, timeLinLayer=0;

void fillMlist(int n_rounds){
  /* Reads in n_rounds matrices from files named "matrix_a.txt" (where a is 1,2,...,). */
  int i, j, m;
  FILE *fp;
  char filename[20];

  Mlist=(int ***)malloc(n_rounds*sizeof(int **));
  for(m=0; m<n_rounds; ++m){
    printf("Reading matrix %d\n",m);
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

vector<LWECiphertext> Chi(vector<LWECiphertext> input, BinFHEContext *ctx){
  /* Computes and returns Chi(input). */
  int i;
  LWECiphertext tmp;
  vector<LWECiphertext> output;
  clock_t tic, toc;
  double timeSoFar;

  tic=clock();
  for(i=0; i<BLOCK_SIZE; ++i){
    if(i%10==0){
      toc=clock();
      timeSoFar=((double)(toc-tic))/CLOCKS_PER_SEC;
      cout<<"Produced "<<i<<" bits of Chi in "<<timeSoFar<<" seconds"<<endl;
    }

    tmp=ctx->EvalNOT(input[(i+1)%BLOCK_SIZE]);//tmp = x_i+1
    tmp=ctx->EvalBinGate(AND,tmp,input[(i+2)%BLOCK_SIZE]);//tmp = x_i+2(x_i+1 + 1)
    tmp=ctx->EvalBinGate(XOR,tmp,input[i]);//tmp=x_i+2(x_i+1 + 1) + x_i
    output.push_back(tmp);
  }

  return output;
}

vector<LWECiphertext> matrixMultiplicationPlain(vector<LWECiphertext> input, int **M, BinFHEContext *ctx, LWEPrivateKey secKey){
  /* Computes and returns input*M in the straight-forward way. */
  int i, j;
  vector<LWECiphertext> output;
  LWECiphertext tmp;
  clock_t tic, toc;
  double timeSoFar;

  tic=clock();
  for(j=0; j<BLOCK_SIZE; ++j){
    if(j%10==0){
      toc=clock();
      timeSoFar=((double)(toc-tic))/CLOCKS_PER_SEC;
      cout<<" passing column "<<j<<", time since start: "<<timeSoFar<<" seconds"<<endl;
    }
    tmp=ctx->Encrypt(secKey,0);
    for(i=0; i<BLOCK_SIZE; ++i){
      if(M[i][j]!=0)
	tmp=ctx->EvalBinGate(XOR,tmp,input[i]);
    }
    output.push_back(tmp);
  }

  return output;
}

vector<LWECiphertext> matrixMultiplicationM4R(vector<LWECiphertext> input, int **M, BinFHEContext *ctx, LWEPrivateKey secKey){
  /* Computes and returns input*M using the method of the four russians.  For BLOCK_SIZE=351, t=5 is optimal for binfhe. */
  int i, j, t=5, nTables, tableSize, tableIndex, p, rest;
  LWECiphertext tmp;
  vector<LWECiphertext> output;
  vector<vector<LWECiphertext>> ciphertextTable;
  double time_table, time_10col, time_Tab;

  //Create lookup-tables
  nTables=BLOCK_SIZE/t;//will be 70 when BLOCKSIZE=351 and t=5
  tableSize=1<<t;
  clock_t tabStart=clock();
  for(i=0; i<1/*nTables*/; ++i){
    vector<LWECiphertext> tmpTab;
    for(j=0; j<tableSize; ++j){
      tmp=ctx->Encrypt(secKey,0);
      tmpTab.push_back(tmp);
    }
    //tmpTab initialized
    for(tableIndex=1; tableIndex<tableSize; ++tableIndex){//nothing to add when tableIndex is 0
      for(p=0; p<t; ++p){
	if((1<<p)&tableIndex)//must add ciphertext p of input block i onto current table index
	  tmpTab[tableIndex]=ctx->EvalBinGate(XOR_FAST,tmpTab[tableIndex],input[t*i+p]);
      }
    }//tmpTab initialized
    ciphertextTable.push_back(tmpTab);
    clock_t singleTab_end=clock();
    time_Tab=((double)(singleTab_end-tabStart))/CLOCKS_PER_SEC;
    cout<<"Table "<<i+1<<" made in "<<time_Tab<<" seconds"<<endl;
  }//all look-up tables initialized
  clock_t tabEnd=clock();
  time_table=((double)(tabEnd-tabStart))/CLOCKS_PER_SEC;
  cout<<"   M4R: constructed "<<nTables<<" ciphertext tables of size "<<tableSize<<" each in "<<time_table<<" seconds"<<endl;
  
  for(i=0; i<BLOCK_SIZE; ++i){
    tmp=ctx->Encrypt(secKey,0);
    output.push_back(tmp);
  }
  //output initialized

  //Do the actual multiplication
  rest=BLOCK_SIZE%t;
  clock_t tic=clock();
  for(j=0; j<BLOCK_SIZE; ++j){//for every column...
    if(j%2==0){
      clock_t toc=clock();
      time_10col=((double)(toc-tic))/CLOCKS_PER_SEC;
      cout<<"       passing column "<<j<<" spent "<<time_10col<<" seconds in main loop"<<endl;
    }
    for(i=0; i<nTables; ++i){//...look up every block of t rows
      tableIndex=0;
      for(p=0; p<t; ++p){
	if(M[(t*i)+p][j]!=0)
	  tableIndex^=(1<<p);
      }//tableIndex indicates element from ciphertexTable to add to output
      if(tableIndex!=0)
	output[j]=ctx->EvalBinGate(XOR_FAST,output[j],ciphertextTable[0/*i*/][tableIndex]);
    }
    //adding the last BLOCK_SIZE mod t elements of input
    for(p=BLOCK_SIZE-rest; p<BLOCK_SIZE; ++p){
      if(M[p][j]!=0)
	output[j]=ctx->EvalBinGate(XOR_FAST,output[j],input[p]);
    }//column j fully processed
  }
    
  return output;
}

vector<LWECiphertext> RastaRound(vector<LWECiphertext> input, int r, BinFHEContext *ctx, LWEPrivateKey secKey){
  /* One Rasta round, with Chi and matrix multiplication (no key or constant addition). */
  vector<LWECiphertext> output, afterChi;

  clock_t chi_start=clock();
  afterChi=Chi(input,ctx);
  clock_t chi_end=clock();
  timeChi+=(double)(chi_end - chi_start)/CLOCKS_PER_SEC;

  clock_t lin_start=clock();
  output=matrixMultiplicationM4R(afterChi,Mlist[r],ctx,secKey);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;

  return output;
}


int main() {
  auto cryptoContext = BinFHEContext();
  vector<LWECiphertext> cipher_block;
  LWECiphertext tmp;
  int i, c;

  cryptoContext.GenerateBinFHEContext(STD128);

  // Generate the secret key
  auto sk = cryptoContext.KeyGen();
  // Generate the bootstrapping keys (refresh and switching keys)
  cryptoContext.BTKeyGen(sk);

  fillMlist(7);
  printf("Matrices read in\n");

  printf("Plaintext: [");
  for(i=0; i<BLOCK_SIZE; ++i){
    c=rand()%2;
    printf("%d ",c);
    tmp=cryptoContext.Encrypt(sk,c);
    cipher_block.push_back(tmp);
  }
  printf("]\n");

  clock_t lin_start=clock();
  cipher_block=matrixMultiplicationM4R(cipher_block,Mlist[0],&cryptoContext,sk);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;
  printf("Full LinLayer took %1.3f seconds\n",timeLinLayer);
  exit(0);

  clock_t chi_start=clock();
  cipher_block=Chi(cipher_block,&cryptoContext);
  clock_t chi_end=clock();
  timeChi+=(double)(chi_end - chi_start)/CLOCKS_PER_SEC;
  printf("Full Chi took %1.3f seconds\n",timeChi);
  exit(0);
  
  for(i=1; i<=6; ++i){
    printf("Entering round %d...\n",i);
    cipher_block=RastaRound(cipher_block,i,&cryptoContext,sk);
    printf("Out of round %d.  Time taken so far: Chi %1.3f seconds  -  Matrix multiplication %1.3f seconds\n",i,timeChi,timeLinLayer);
  }
}
