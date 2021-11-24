/* Code implementing Rasta with 329-bit block in HElib, using one ciphertext for whole block.
 Used only for timing purposes, addition of constants and keys not implemented. */

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
vector<MatMulFullExec> MAT;


Ctxt Chi(Ctxt input){
  /* Returns Chi(input). */
  Ctxt i0(input), i1(input), i2(input);
  const EncryptedArray& ea = input.getContext().getEA();

  ea.rotate(i1,1);
  ea.rotate(i2,2);

  i1*=i2;
  i1+=i2;
  i1+=i0;
  
  return i1;
}

Ctxt LinLayer(Ctxt input, int round){
  /* Linear layer of Rasta, multiply with random matrix i. */
  Ctxt out(input);

  MAT[round].mul(out);

  return out;
}

Ctxt Round(Ctxt input, int r){
  /* Returns the state after one round of Rasta. */
  Ctxt out(input);
  int i;
 
  clock_t chi_start=clock();
  out=Chi(input);
  clock_t chi_end=clock();
  timeChi+=(double)(chi_end - chi_start)/CLOCKS_PER_SEC;
  
  clock_t lin_start=clock();
  out=LinLayer(out,r);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;
  
  return out;
}


int main(int argc, char* argv[])
{
  /*  Example of BGV scheme  */
  int i, j;
  double secs_taken;
  
  // Plaintext prime modulus
  unsigned long p = 2;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m=30269;
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 500;
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
  long nslots = ea.size();
  cout << "Number of slots: " << nslots << endl;

  clock_t setup_end = clock();
  secs_taken = (double)(setup_end - setup_begin) / CLOCKS_PER_SEC;
  cout<<"Set-up took "<<secs_taken<<" seconds"<<endl;

  //Reading in random matrices
  char filename[80];

  clock_t readMat_start=clock();
  for(i=0; i<7; ++i){
    sprintf(filename,"matrix_%d.txt",i+1);
    matrixFromFile nyM(ea,filename);
    MatMulFullExec MMFE(nyM);
    MAT.push_back(MMFE);
  }
  clock_t readMat_end=clock();
  secs_taken = (double)(readMat_end - readMat_start) / CLOCKS_PER_SEC;
  cout<<"Reading in the 7 matrices took "<<secs_taken<<" seconds"<<endl;
  
  
  // Create plaintext vector of long with nslots elements each
  Ptxt<BGV> ptxt(context);
  // Fill them with random bits
  for (i = 0; i < ptxt.size(); ++i) 
    ptxt[i] = (long)random()%2;
  //cout<<"ptxt = "<<ptxt<<endl<<endl;

  // Create a ciphertext object
  Ctxt ctxt(public_key);
  // Encrypt the plaintext using the public_key
  public_key.Encrypt(ctxt,ptxt);
  cout<<"Initial bit capacity of word: "<<ctxt.bitCapacity()<<endl;

  /* Homomorphic evaluation of Rasta */
  //Initial linear transformation
  clock_t lin_start=clock();
  ctxt=LinLayer(ctxt,0);
  clock_t lin_end=clock();
  timeLinLayer+=(double)(lin_end - lin_start)/CLOCKS_PER_SEC;

  for(i=1; i<=6; ++i){
    cout<<"Entering round "<<i<<endl;
    ctxt=Round(ctxt,i);
  }
  cout<<"Bit capacity of word after Rasta: "<<ctxt.bitCapacity()<<endl;

  cout<<"Time used in Chi-transformations: "<<timeChi<<" seconds"<<endl;
  cout<<"Time used in LinLayer: "<<timeLinLayer<<" seconds"<<endl;
}
