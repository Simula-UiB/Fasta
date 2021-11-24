/* Implementation of Rasta using TFHE.  
   Used only for timing purposes, key and constant additions are not implemented. */

#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <stdlib.h>
#include <time.h>

int BLOCK_SIZE=351;
TFheGateBootstrappingParameterSet *context;
TFheGateBootstrappingSecretKeySet *secret_key;
int ***Mlist;
double timeChi=0, timeLinLayer=0;

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

LweSample *Chi(LweSample *input, const TFheGateBootstrappingCloudKeySet cloud_key){
  /* Computes and returns Chi(input). */
  int i;
  LweSample *tmp, *output;

  output=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);
  tmp=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);

  for(i=0; i<BLOCK_SIZE; ++i){
    bootsANDYN(tmp+i,input+((i+2)%BLOCK_SIZE),input+((i+1)%BLOCK_SIZE),&cloud_key);
    bootsXOR(output+i,tmp+i,input+i,&cloud_key);
  }

  return output;
}

LweSample *matrixMultiplicationPlain(LweSample *input, int **M, const TFheGateBootstrappingCloudKeySet cloud_key){
  /* Computes and returns input*M in the straigh-forward way. */
  int i, j;
  LweSample *output;

  output=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);
  for(i=0; i<BLOCK_SIZE; ++i)
    bootsSymEncrypt(output+i,0,secret_key);
  //output initialized
  for(j=0; j<BLOCK_SIZE; ++j){
    for(i=0; i<BLOCK_SIZE; ++i){
      if(M[i][j]!=0)
	bootsXOR(output+j,output+j,input+i,&cloud_key);
    }
  }

  return output;
}

LweSample *matrixMultiplicationM4R(LweSample *input, int **M, const TFheGateBootstrappingCloudKeySet cloud_key){
  /* Computes and returns input*M using the method of the four russians.  For BLOCK_SIZE=351, t=5 is optimal. */
  int i, j, t=5/*optimal*/, nTables, tableSize, tableIndex, p, rest;
  LweSample *output, **ciphertextTable;

  //Create lookup-tables
  nTables=BLOCK_SIZE/t;//will be 70 when BLOCKSIZE=351 and t=5
  tableSize=1<<t;
  ciphertextTable=(LweSample **)malloc(nTables*sizeof(LweSample *));
  for(i=0; i<nTables; ++i){
    ciphertextTable[i]=new_gate_bootstrapping_ciphertext_array(tableSize,context);
    for(j=0; j<tableSize; ++j)
      bootsSymEncrypt(ciphertextTable[i]+j,0,secret_key);
    //ciphertextTable[i] initialized
    for(tableIndex=1; tableIndex<tableSize; ++tableIndex){//nothing to add when tableIndex is 0
      for(p=0; p<t; ++p){
	if((1<<p)&tableIndex)//must add ciphertext p of input block i onto current table index
	  bootsXOR(ciphertextTable[i]+tableIndex,ciphertextTable[i]+tableIndex,input+(t*i+p),&cloud_key);
      }
    }//ciphertextTable[i] initialized
  }//all look-up tables initialized
  
  output=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);
  for(i=0; i<BLOCK_SIZE; ++i)
    bootsSymEncrypt(output+i,0,secret_key);
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
	bootsXOR(output+j,output+j,ciphertextTable[i]+tableIndex,&cloud_key);
    }
    //adding the last BLOCK_SIZE mod t elements of input
    for(p=BLOCK_SIZE-rest; p<BLOCK_SIZE; ++p){
      if(M[p][j]!=0)
 	bootsXOR(output+j,output+j,input+p,&cloud_key);
    }//column j fully processed
  }
    
  return output;
}

LweSample *RastaRound(LweSample *input, int r, const TFheGateBootstrappingCloudKeySet cloud_key){
  /* One Rasta round, with Chi and matrix multiplication (no key or constant addition). */
  LweSample *output, *afterChi;

  output=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);
  afterChi=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);

  clock_t chi_start=clock();
  afterChi=Chi(input,cloud_key);
  clock_t chi_end=clock();
  timeChi+=(double)(chi_end - chi_start)/CLOCKS_PER_SEC;

  clock_t lin_start=clock();
  output=matrixMultiplicationM4R(afterChi,Mlist[r],cloud_key);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;

  return output;
}


int main(int argc, char *argv[]){
  LweSample *roundInput, *roundOutput;
  int i, c;

  context=new_default_gate_bootstrapping_parameters(128);//want 128-bit security;
  secret_key=new_random_gate_bootstrapping_secret_keyset(context);
  const TFheGateBootstrappingCloudKeySet cloud_key=secret_key->cloud;

  fillMlist(7);
  
  roundInput=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);
  roundOutput=new_gate_bootstrapping_ciphertext_array(BLOCK_SIZE,context);

  printf("Plaintext: [");
  for(i=0; i<BLOCK_SIZE; ++i){
    c=rand()%2;
    printf("%d ",c);
    bootsSymEncrypt(roundInput+i,c,secret_key);
  }
  printf("]\n");
  
  /* Homomorphic evaluation of Rasta */
  clock_t lin_start=clock();
  roundOutput=matrixMultiplicationM4R(roundInput,Mlist[0],cloud_key);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;
  printf("Initial linear transformation took %1.3f seconds\n",timeLinLayer);
  
  roundInput=roundOutput;
  for(i=1; i<=6; ++i){
    printf("Entering round %d...\n",i);
    roundOutput=RastaRound(roundInput,i,cloud_key);
    printf("Out of round %d.  Time taken so far: Chi %1.3f seconds  -  Matrix multiplication %1.3f seconds\n",i,timeChi,timeLinLayer);
    roundInput=roundOutput;
  }
  
  return 0;
}
