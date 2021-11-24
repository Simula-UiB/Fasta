/* For generating random 329 x 329 matrices.  NOTE: not checking for invertibility just yet! */

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[]){
  int i, j, mat_nr, bit;
  FILE *fp;
  char filename[80];

  for(mat_nr=7; mat_nr<=7; ++mat_nr){
    sprintf(filename,"matrix_%d.txt",mat_nr);
    fp=fopen(filename,"w");
    for(i=0; i<329; ++i){
      for(j=0; j<329; ++j){
	bit=random()%2;
	fprintf(fp,"%d ",bit);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  } 
}
