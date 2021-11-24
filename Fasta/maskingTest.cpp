/* Test for checking how expensive plaintext-ciphertext multiplication is in HElib. */

#include <stdlib.h>
#include <iostream>

using namespace std;

#include <helib/helib.h>

using namespace helib;

int main(int argc, char* argv[])
{
  double secs_taken;
  int i;
  
  // Plaintext prime modulus
  unsigned long p = 2;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m=30269;
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 300;
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

  Ptxt<BGV> p1(context), mask(context);
  // Fill p1 with random bits
  for (i = 0; i < p1.size(); ++i) 
    p1[i] = (long)random()%2;
  //Fill mask with 1's in first 100 slots, zero in rest
  for(i=0; i<100; ++i)
    mask[i]=1;
  for(i=100; i<mask.size(); ++i)
    mask[i]=0;

  // Create a ciphertext object
  Ctxt c1(public_key);
  // Encrypt p1 using the public_key
  public_key.Encrypt(c1,p1);
  cout<<"Initial bit capacity of c1: "<<c1.bitCapacity()<<endl;

  for(i=0; i<10; ++i){
    c1*=mask;
    cout<<"bit capacity after masking "<<i+1<<" times: "<<c1.bitCapacity()<<endl;
    //c1*=c1;
    //cout<<"bit capacity after squaring c1 "<<i+1<<" times: "<<c1.bitCapacity()<<endl;
  }
}
