/* Example program that encrypts a plaintext file using Fasta and writes the ciphertext to file.
Then it reads in the ciphertext file and decrypts it to verify that we get back the same plaintext.
Key is generated at random on every run. */

#include "Fasta.c"

struct plaintext *readPlaintextFile(const char *filename){
  /* Reads in plaintext from the given file and returns it in a plaintext object. */
  FILE *fp;
  struct plaintext *PT;
  int i, j;

  fp=fopen(filename,"r");
  if(fp==NULL){printf("Could not open the file %s\n",filename); exit(0);}

  PT=(struct plaintext *)malloc(sizeof(struct plaintext *));
  fscanf(fp,"%d\n",&PT->numBlocks);
  PT->block=(struct bitVector **)malloc(PT->numBlocks*sizeof(struct bitVector *));
  for(i=0; i<PT->numBlocks; ++i){
    PT->block[i]=newBitVector(1645);
    for(j=51; j>=0; --j)
      fscanf(fp,"%8x ",PT->block[i]->v+j);
    fscanf(fp,"\n");
  }
  fclose(fp);

  return PT;
}

struct ciphertext *readCiphertextFile(const char *filename){
  /* Reads in ciphertext from the given file and returns it in a ciphertext object. */
  FILE *fp;
  struct ciphertext *CT;
  int i, j;

  fp=fopen(filename,"r");
  if(fp==NULL){printf("Could not open the file %s\n",filename); exit(0);}

  CT=(struct ciphertext *)malloc(sizeof(struct ciphertext *));
  fscanf(fp,"%d\n",&CT->numBlocks);
  CT->block=(struct bitVector **)malloc(CT->numBlocks*sizeof(struct bitVector *));
  for(i=0; i<CT->numBlocks; ++i){
    CT->block[i]=newBitVector(1645);
    for(j=51; j>=0; --j)
      fscanf(fp,"%8x ",CT->block[i]->v+j);
    fscanf(fp,"\n");
  }
  fclose(fp);

  return CT;
}

void printPlaintextFile(struct plaintext *PT, const char *filename){
  /* Prints the plaintext PT to a file with the given name. */
  FILE *fp;
  int i, j;

   fp=fopen(filename,"w");
   fprintf(fp,"%d\n",PT->numBlocks);
   for(i=0; i<PT->numBlocks; ++i){
     for(j=51; j>=0; --j)
       fprintf(fp,"%08x ",PT->block[i]->v[j]);
     fprintf(fp,"\n");
   }
   fclose(fp);
}

void printCiphertextFile(struct ciphertext *CT, const char *filename){
  /* Prints the ciphertext CT to a file with the given name. */
  FILE *fp;
  int i, j;

   fp=fopen(filename,"w");
   fprintf(fp,"%d\n",CT->numBlocks);
   for(i=0; i<CT->numBlocks; ++i){
     for(j=51; j>=0; --j)
       fprintf(fp,"%08x ",CT->block[i]->v[j]);
     fprintf(fp,"\n");
   }
   fclose(fp);
}


int main(int argc, char *argv[]){
  struct plaintext *ptxt, *decrypted;
  struct ciphertext *ctxt;
  struct bitVector *K;
  
  K=randomBitVector(329);
  ptxt=readPlaintextFile("originalPlaintext.txt");
  printf("Plaintext read in\n");

  ctxt=encrypt(ptxt,K);  
  printCiphertextFile(ctxt,"ciphertext.txt");
  printf("Encryption of plaintext written to ciphertext.txt\n");

  decrypted=decrypt(ctxt,K);
  printPlaintextFile(decrypted,"decryptedPlaintext.txt");
  printf("Ciphertext decrypted and written to decryptedPlaintext.txt\n"); 
}

