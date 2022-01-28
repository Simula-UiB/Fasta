/* File for creating, deleting and doing operations on bitVectors */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

typedef unsigned int u32;
typedef unsigned char u8;

u8 weight[256]={0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,
		2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,
		4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,
		3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
		4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

struct bitVector{
  u32 *v;
  int length, wl;
  u8 con;
};

struct bitVector *newBitVector(int len){
  /* Returns a pointer to a bitVector of length len, initialized to 0-vector. */
  struct bitVector *bv;

  bv=(struct bitVector *)malloc(sizeof(struct bitVector));
  bv->length=len;
  bv->wl=(bv->length+31)>>5;
  bv->v=(u32 *)calloc(bv->wl,sizeof(u32));
  bv->con=0;
  return bv;
}

void deleteBitVector(struct bitVector *bv){
  /* Frees all memory allocated to bv. */
  free(bv->v);
  free(bv);
}

void setBit(struct bitVector *bv, int index){
  /* Sets bit index in bv to 1. */
  if(index>bv->length-1 || index<0){printf("(bitVector->setBit)Can't set bit %d in vector of length %d\n",index,bv->length);exit(0);}
  bv->v[index>>5]|=1<<(index&0x1f);
}

struct bitVector *randomBitVector(int len){
  /* Returns a random bitVector of length len. */
  struct bitVector *bv;
  int i;
  u32 mask, rr;

  bv=newBitVector(len);
  for(i=0; i<bv->wl; ++i){
    bv->v[i]=(u32)random();
    if(32*i+31<len){
      rr=(u32)random();
      if(rr&1)
	setBit(bv,32*i+31);//because highest set bit from random() is always 0
    }
  }
  mask=(1<<(len&0x1f))-1;
  if(mask)
    bv->v[bv->wl-1]&=mask;//setting bits higher than len to 0

  return bv;
}

u8 isSet(struct bitVector *bv, int index){
  /* Returns 1 if bit index in bv is set, 0 otherwise. */
  if(bv->length<index){
    printf("(bitVector->isSet)Trying to check bit %d in vector of length %d\n",index,bv->length);
    return 0;
  }
  if(bv->v[index>>5]&(1<<(index&0x1f)))
    return 1;
  else
    return 0;
}

void unSetBit(struct bitVector *bv, int index){
  /* Sets bit index in bv to 0. */
  if(index>bv->length-1 || index<0){printf("(bitVector->unSetBit)Can't un-set bit %d in vector of length %d\n",index,bv->length);exit(0);}
  if(isSet(bv,index))//ensures the bit really is set
     bv->v[index>>5]^=1<<(index&0x1f);
}

void flipBit(struct bitVector *bv, int index){
  /* Flips bit index in bv. */
  if(index>bv->length-1){printf("(bitVector->flipBit)Can't flip bit %d in vector of length %d\n",index,bv->length);exit(0);}
  if(isSet(bv,index))
    unSetBit(bv,index);
  else
    setBit(bv,index);
}


void wipe(struct bitVector *bv){
  /* Sets all bits in bv to 0. */
  int i;

  for(i=0; i<bv->wl; ++i)
    bv->v[i]=0;
  bv->con=0;
}

void embed(struct bitVector *bv, struct bitVector *ebv){
  /* Copies bv of shorter or equal length than ebv, into ebv. */
  int i;

  if(bv->length>ebv->length){printf("(embed)Can't embed vector of length %d into vector of length %d.\n",bv->length,ebv->length); exit(0);}
  wipe(ebv);
  for(i=0; i<bv->wl; ++i)
    ebv->v[i]=bv->v[i];
}

int highestSetBit(struct bitVector *bv){
  /* Returns the index of the highest set bit in bv.  Returns -1 if bv=0. */
  int i, j;
  u32 m1, m2;

  i=bv->wl-1;
  while(i>=0 && bv->v[i]==0)
    i--;
  if(i<0)
    return -1;

  j=0;
  m1=bv->v[i]&0xffff0000;
  if(m1){
    j+=16;
    m1>>=16;
  }
  else
    m1=bv->v[i]&0xffff;

  m2=m1&0xff00;
  if(m2){
    j+=8;
    m2>>=8;
  }
  else
    m2=m1&0xff;

  m1=m2&0xf0;
  if(m1){
    j+=4;
    m1>>=4;
  }
  else
    m1=m2&0xf;

  m2=m1&0xc;
  if(m2){
    j+=2;
    m2>>=2;
  }
  else
    m2=m1&0x3;

  if(m2>1)
    j++;

  return 32*i+j;
}

int lowestSetBit(struct bitVector *bv){
  /* Returns index of lowest set bit in bv.  Returns -1 if bv=0. */ 
  int i=0, j;
  u32 m1, m2;

  while(i<bv->wl && bv->v[i]==0)
    ++i;
  if(i==bv->wl)
    return -1;

  m1=bv->v[i]&0xffff;
  j=0;
  if(m1==0){
    j+=16;
    m1=(bv->v[i]&0xffff0000)>>16;
  }

  m2=m1&0xff;
  if(m2==0){
    j+=8;
    m2=(m1&0xff00)>>8;
  }

  m1=m2&0xf;
  if(m1==0){
    j+=4;
    m1=(m2&0xf0)>>4;
  }

  m2=m1&0x3;
  if(m2==0){
    j+=2;
    m2=(m1&0xc)>>2;
  }

  if(m2==2)
    j++;

  return 32*i+j;
}

void copyV1toV2(struct bitVector *v1, struct bitVector *v2){
  /* Makes v2 a copy of v1. */
  int i;

  if(v1->length!=v2->length){printf("(copyV1toV2)Trying to copy vector of length %d into vector of legth %d\n",v1->length,v2->length);exit(0);}
  for(i=0; i<v1->wl; ++i)
    v2->v[i]=v1->v[i];
  v2->con=v1->con;
}

u8 equal(struct bitVector *v1, struct bitVector *v2){
  /* Returns 1 if v1 and v2 are equal, and 0 otherwise. */
  int i;

  if(v1->length!=v2->length)
    return 0;
  if(v1->con!=v2->con)
    return 0;
  for(i=0; i<v1->wl; ++i){
    if(v1->v[i]^v2->v[i])
      return 0;
  }
  return 1;
}

int hammingWeight(struct bitVector *bv){
  /* Returns the Hamming weight of bv. */
  int i, w=0;

  for(i=0; i<bv->wl; ++i){
    w+=weight[bv->v[i]&0xff];
    w+=weight[(bv->v[i]>>8)&0xff];
    w+=weight[(bv->v[i]>>16)&0xff];
    w+=weight[bv->v[i]>>24];
  }
  return w;
}


struct bitVector *v1ANDv2(struct bitVector *v1, struct bitVector *v2){
  /* Returns the bitwise AND of v1 and v2. */  
  struct bitVector *ret;
  int i;

  if(v1->length!=v2->length){printf("(bitVector->v1ANDv2)Vectors have different lengths (%d and %d)\n",v1->length,v2->length);exit(0);}
  ret=newBitVector(v1->length);
  for(i=0; i<ret->wl; ++i)
    ret->v[i]=v1->v[i]&v2->v[i];
  ret->con=v1->con&v2->con;

  return ret;
}

void ORv1toV2(struct bitVector *v1, struct bitVector *v2){
  /* ORs v1 onto v2. */
  int i;

  if(v1->length!=v2->length){printf("(bitVector->ORv1toV2)Vectors not of same lengths (%d and %d)\n",v1->length,v2->length);exit(0);}
  for(i=0; i<v1->wl; ++i)
    v2->v[i]|=v1->v[i];
  v2->con|=v1->con;
}

void addV1toV2(struct bitVector *v1, struct bitVector *v2){
  /* Adds v1 to v2 using XOR. */
  int i;

  if(v1->length!=v2->length){printf("(bitVector->addV1toV2)Vectors not of same lengths (%d and %d)\n",v1->length,v2->length);exit(0);}
  for(i=0; i<v1->wl; ++i)
    v2->v[i]^=v1->v[i];
  v2->con^=v1->con;
}

void rotLeft(struct bitVector *bv, int rr){
  /* Rotates bv rr positions to the left, cyclically */
  int i, j;
  struct bitVector *topBits;

  if(rr==0)
    return;
  
  topBits=newBitVector(rr);
  j=0;
  for(i=bv->length-rr; i<bv->length; ++i){
    if(isSet(bv,i))
      setBit(topBits,j);
    j++;
  }
  for(i=bv->length-1; i>=rr; --i){
    if(isSet(bv,i-rr))
      setBit(bv,i);
    else
      unSetBit(bv,i);
  }
  for(i=0; i<rr; ++i){
    if(isSet(topBits,i))
      setBit(bv,i);
    else
      unSetBit(bv,i);
  }
  deleteBitVector(topBits);
}


int hammingDistance(struct bitVector *bv1, struct bitVector *bv2){
  /* Returns the Hamming distance between bv1 and bv2. */
  struct bitVector *s;
  int retd;

  if(bv1->length!=bv2->length){
    printf("(bitVector->hammingDistance)Trying to compute distance between vectors of length %d and %d\n",bv1->length,bv2->length);
    exit(0);
  }

  s=newBitVector(bv1->length);
  copyV1toV2(bv1,s);
  addV1toV2(bv2,s);
  retd=hammingWeight(s);
  deleteBitVector(s);

  return retd;
}

void printVectorHEX(struct bitVector *bv){
  int i;

  for(i=bv->wl-1; i>=0; --i)
    printf("%08x ",bv->v[i]);
}

void printVectorBits(struct bitVector *bv){
  int i;

  printf("[");
  for(i=bv->length-1; i>=0; --i){
    if(isSet(bv,i))
      printf("1");
    else
      printf(".");
  }
  printf("]");
}
