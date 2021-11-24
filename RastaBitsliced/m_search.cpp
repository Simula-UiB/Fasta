#include <stdlib.h>
#include <iostream>
#include <stdio.h>

using namespace std;

#include <helib/helib.h>

using namespace helib;


int main(int argc, char* argv[])
{
  /*  Example of BGV scheme  */
  int i, j;
  double secs_taken;
  
  // Plaintext prime modulus
  unsigned long p = 2;
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m, start;
  // Hensel lifting (default = 1)
  unsigned long r = 1;
  // Number of bits of the modulus chain
  unsigned long bits = 500;
  // Number of columns of Key-Switching matrix (default = 2 or 3)
  unsigned long c = 2;

  FILE *fp;

  start=atoi(argv[1]);

  //for(m=start; m>start+20; m+=2){
  m=start;
  while(m<start+10000){
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

    double secLev=context.securityLevel();

    if(secLev>=126.0)
      printf("m=%d gives security level %1.3f\n",(int)m,secLev);
    m+=2;
    /*
    if((nslots&1) && (nslots>=5)){
      fp=fopen("potential_m_values.txt","a");
      fprintf(fp,"m=%d gives %d slots and security level %1.3f\n",(int)m,(int)nslots,context.securityLevel());
      fclose(fp);
      }*/
  }
}
