#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#define NUMALPHA 80
#define NUMTEMPL 80
#define NETWIDTH 20
#define NETHEIGHT 120
#define BOXSIZE 30
#define TEMPSIZE 5

#define OVERLAP 0
#define ABSORB 1

double PENALTY=1;
double MUTATION=1e-2;
#define DELETION 1e-2
#define DUPLICATION 1e-2

#define GetIdx(x,y,z) ((x)+(y)*NETWIDTH+(z)*NETWIDTH*NETWIDTH)

int *Box;
char *Light;
int *SolarGrid;
int *AbsGrid;

typedef struct
{
   char *Str;
   int len;
} String;

typedef struct
{
   int type; // 0 = inactive, otherwise its a type of neuron
   int conpt; // binding sites. A +1 binds to a -1, a +2 to a -2, etc
   int idx;
} Node;

typedef struct
{
   Node Array[TEMPSIZE][TEMPSIZE][TEMPSIZE]; // Array of nodes that belong to the template   
} Template;

typedef struct
{
   Node *Array;
} NetworkPreimage;

typedef struct
{
   int templateid[NUMALPHA];
} Mapping;

typedef struct
{
   double act, inp;
   int type;
   int x,y;
} Neuron;

typedef struct
{
   int n1,n2;
} Link;

typedef struct
{
   int NNodes;
   Neuron *N;
   int NLinks;
   Link *L;
} Network;

typedef struct
{
   int idx;
   Network Net;
   NetworkPreimage P;
   String Genome;
   Mapping M;
   double fitness;
   int active;
   double x,y;
   int solar;
   int overlap;
   int cells;
} Organism;

int NumOrg;
Organism *Org;
int oidx=0;

Template BaseTemplates[NUMTEMPL];

String AllocString(int len)
{
   String S;
   
   S.len=len; S.Str=(char*)malloc(len);
   
   return S;
}

void FreeString(String S)
{
   if (S.Str!=NULL)
     free(S.Str);
}

int StringCmp(String S1, String S2)
{
   int i;
   
   if (S1.len!=S2.len) return 0;
   
   for (i=0;i<S1.len;i++)
     if (S1.Str[i]!=S2.Str[i]) return 0;
   
   return 1;
}

String ConcatStrings(String S1, String S2)
{
   String S3=AllocString(S1.len+S2.len);
   int i;
   
   for (i=0;i<S1.len;i++)
     S3.Str[i]=S1.Str[i];
   
   for (i=S1.len;i<S3.len;i++)
     S3.Str[i]=S2.Str[i-S1.len];
   
   return S3;
}

void InsertString(String *Base, String Fragment, int ofs)
{
   String S1,S2,S3,S4;
   int i;
   
   S1=AllocString(ofs);
   for (i=0;i<ofs;i++)
     S1.Str[i]=Base->Str[i];
   
   S2=AllocString(Base->len-ofs);
   for (i=ofs;i<Base->len;i++)
     S2.Str[i-ofs]=Base->Str[i];
   
   S3=ConcatStrings(S1,Fragment);
   S4=ConcatStrings(S3,S2);
   
   FreeString(S1); FreeString(S2); FreeString(S3);
   FreeString(*Base);
   *Base=S4;
}

void ReplacePart(String *S, int ofs, int len, String NewS)
{
   int i,j;
   String Sub1,Sub2;
   String Result;
   
   if (ofs>0)
   {
      Sub1=AllocString(ofs);
   
      for (i=0;i<ofs;i++)
	Sub1.Str[i]=S->Str[i];
      
      Sub2=AllocString(S->len-ofs-len);
      for (i=ofs+len;i<S->len;i++)
	Sub2.Str[i-ofs-len]=S->Str[i];
      
      FreeString(*S);
      if (NewS.len>0)
      {
	 Result=ConcatStrings(Sub1,NewS);
	 *S=ConcatStrings(Result,Sub2);
      
	 FreeString(Result); FreeString(Sub1); FreeString(Sub2);
      }
      else
      {
	 *S=ConcatStrings(Sub1,Sub2);
      
	 FreeString(Sub1); FreeString(Sub2);
      }            
   }
   else
   {
      Sub2=AllocString(S->len-ofs-len);
      for (i=ofs+len;i<S->len;i++)
	Sub2.Str[i-ofs-len]=S->Str[i];
      
      FreeString(*S);
      if (NewS.len>0)
      {
	 *S=ConcatStrings(NewS,Sub2);      
	 FreeString(Sub2);
      }
      else
      {
	 *S=Sub2;      
      }      
   }   
}

int HasSubStr(String S, String SubStr, int n)
{
   int i,j,k;
   int match;
   
   k=0;
   for (i=0;i<S.len;i++)
   {
      match=1;
      for (j=0;(match)&&(j<SubStr.len)&&(i+j<S.len);j++)
      {
	 if (S.Str[i+j]!=SubStr.Str[j])
	   match=0;	
      }      
      if ((j==SubStr.len)&&(match))
      {
	 if (k==n) return i;
	 else k++;
      }      
   }   
   
   return -1;
}

int NumSubStrs(String S, String SubStr)
{
   int i,j,k;
   int match;
   
   k=0;
   for (i=0;i<S.len;i++)
   {
      match=1;
      for (j=0;(match)&&(j<SubStr.len)&&(i+j<S.len);j++)
      {
	 if (S.Str[i+j]!=SubStr.Str[j])
	   match=0;	
      }      
      if ((j==SubStr.len)&&(match))
      {
	 k++;
      }      
   }   
   
   return k;
}

Template GetRotatedTemplate(Template T, int dir)
{
   Template TP;
   int x,y,z,xm,ym,zm;
   
   if (dir==0) return T;
   
   for (z=0;z<TEMPSIZE;z++)
   for (y=0;y<TEMPSIZE;y++)
     for (x=0;x<TEMPSIZE;x++)
     {
	switch (dir)
	{
	 case 1: 
	   xm=y;
	   ym=TEMPSIZE-x-1;
	   zm=z;
	   break;
	 case 2:
	   ym=TEMPSIZE-y-1;
	   xm=TEMPSIZE-x-1;
	   zm=z;
	   break;
	 case 3:
	   xm=TEMPSIZE-y-1;
	   ym=x;
	   zm=z;
	   break;
	 case 4:
	   zm=y;
	   ym=TEMPSIZE-z-1;
	   xm=x;	   
	   break;
	 case 5:
	   zm=x;
	   xm=TEMPSIZE-z-1;
	   ym=y;
	   break;
	}
	
	TP.Array[xm][ym][zm]=T.Array[x][y][z];
     }   
   
   return TP;
}

void AddTemplate(Template T, NetworkPreimage *N)
{
   int x,y,xm,ym,xm2,ym2,xm3,ym3,z,zm,zm2,zm3,bestx,besty,bestz;
   NetworkPreimage N2;
   Template Tprime;
   int dir;
   int overlap, bestdir, bestoverlap;
   
   N2.Array=(Node*)malloc(NETWIDTH*NETWIDTH*NETHEIGHT*sizeof(Node));
   memcpy(N2.Array,N->Array,NETWIDTH*NETWIDTH*NETHEIGHT*sizeof(Node));

   for (z=0;z<NETHEIGHT;z++)
   for (y=0;y<NETWIDTH;y++)
   for (x=0;x<NETWIDTH;x++)
   {
      if (N->Array[GetIdx(x,y,z)].type>0)
	if (N->Array[GetIdx(x,y,z)].conpt!=0)
	{
           bestoverlap=TEMPSIZE*TEMPSIZE;
           bestx=besty=bestz=-1;
	   for (zm=0;zm<TEMPSIZE;zm++)
  	   for (ym=0;ym<TEMPSIZE;ym++)
	   for (xm=0;xm<TEMPSIZE;xm++)
	   {
	      if (T.Array[xm][ym][zm].type>0)
		if (T.Array[xm][ym][zm].conpt==
		    -N2.Array[GetIdx(x,y,z)].conpt)
		{
		   overlap=0;
		   for (zm2=0;zm2<TEMPSIZE;zm2++)
		   for (ym2=0;ym2<TEMPSIZE;ym2++)
		   for (xm2=0;xm2<TEMPSIZE;xm2++)
		   {
		      xm3=(xm2-xm)+x; ym3=(ym2-ym)+y;
		      zm3=(zm2-zm)+z;
		      
		      if ((zm3>=0)&&(xm3>=0)&&(ym3>=0)&&
			  (xm3<NETWIDTH)&&(ym3<NETWIDTH)&&
			  (zm3<NETHEIGHT))
		      {			       
			 if ((N2.Array[GetIdx(xm3,ym3,zm3)].type)&&
			     (T.Array[xm2][ym2][zm2].type))
			   overlap++;
		      }
		   }
		   if (overlap<bestoverlap)
		   {
		      bestoverlap=overlap;
		      bestx=xm; besty=ym; bestz=zm;
		   }
		}		    
	   }	   
		  
	   if (bestoverlap<=3)
	   {			  			  
	      // Add a copy here...
	      for (zm2=0;zm2<TEMPSIZE;zm2++)
		for (ym2=0;ym2<TEMPSIZE;ym2++)
		  for (xm2=0;xm2<TEMPSIZE;xm2++)
		  {
		     xm3=(xm2-bestx)+x; ym3=(ym2-besty)+y;
		     zm3=(zm2-bestz)+z;
		     
		     if ((zm3>=0)&&(xm3>=0)&&(ym3>=0)&&
			 (xm3<NETWIDTH)&&(ym3<NETWIDTH)&&
			 (zm3<NETHEIGHT))
		     {			       
			if (!N2.Array[GetIdx(xm3,ym3,zm3)].type)
			{
			   N2.Array[GetIdx(xm3,ym3,zm3)]=
			     T.Array[xm2][ym2][zm2];
			}
		     }			       
		  }	
	   }		       
	   N2.Array[GetIdx(x,y,z)].conpt=0;
	}		  		  
   }	     
   
   memcpy(N->Array,N2.Array,NETWIDTH*NETWIDTH*NETHEIGHT*sizeof(Node));
   free(N2.Array);
}

void FreeNetwork(Network *N)
{
   free(N->N); N->N=NULL;
   free(N->L); N->L=NULL;
}

void FreeOrganism(Organism *O)
{
//   FreeNetwork(&O->Net);
   FreeString(O->Genome);
   free(O->P.Array);
}

NetworkPreimage GenomeToPreimage(String G, Mapping M)
{
   NetworkPreimage N;
   int idx;
   int templateidx;
   
   memset(&N,0,sizeof(NetworkPreimage));
   
   N.Array=(Node*)malloc(sizeof(Node)*NETWIDTH*NETWIDTH*NETHEIGHT);
   memset(N.Array,0,sizeof(Node)*NETWIDTH*NETWIDTH*NETHEIGHT);
   
   N.Array[GetIdx((NETWIDTH)/2,((NETWIDTH)/2),(NETHEIGHT-1))].conpt=1;
   N.Array[GetIdx((NETWIDTH)/2,((NETWIDTH)/2),(NETHEIGHT-1))].type=1;
   
   for (idx=0;idx<G.len;idx++)
   {
      templateidx=G.Str[idx]-'0';
         
      AddTemplate(BaseTemplates[templateidx],&N);
   }   
   
   return N;
}

	   
void GenerateBaseTemplates(int seed)
{
   int i,x,y,z,xm,ym,zm,n;
   Template T;
   int INVP;
   int count, ncount;
   int ccount;
   
   srand(seed);
   
   for (i=0;i<NUMTEMPL;i++)
   {  
      memset(&T,0,sizeof(Template));
   
      count=1;
      T.Array[1+rand()%(TEMPSIZE-2)][1+rand()%(TEMPSIZE-2)][1+rand()%(TEMPSIZE-2)].type=1;
      ncount=0;
      ccount=3+rand()%4;
      while (ncount<ccount)
      {
	 x=rand()%TEMPSIZE; y=rand()%TEMPSIZE; z=rand()%TEMPSIZE;
	 if (!T.Array[x][y][z].type)
	 {	    
	    n=0;
	    for (zm=z-1;zm<=z+1;zm++)
	    for (ym=y-1;ym<=y+1;ym++)
	      for (xm=x-1;xm<=x+1;xm++)
		if ( (zm>=0)&&(zm<TEMPSIZE)&&(xm>=0)&&(ym>=0)&&(xm<TEMPSIZE)&&(ym<TEMPSIZE))
		  if ((abs(xm-x)+abs(ym-y)+abs(zm-z))==1/*(xm!=x)||(ym!=y)||(zm!=z)*/)
		    n+=(T.Array[xm][ym][zm].type!=0);
	 
	    if (n>0)
	    {
	       T.Array[x][y][z].type=1;
               //if (rand()%5==0) T.Array[x][y][z].type++;

	       if ((z==0)||(x==0)||(y==0)||(z==TEMPSIZE-1)||(x==TEMPSIZE-1)||(y==TEMPSIZE-1))
	       {
		  T.Array[x][y][z].conpt=(rand()%2+1)*(2*(rand()%2)-1);
		  ncount++;
	       }	    
	       count++;
	    }
	 }	 
      }      
      BaseTemplates[i]=T;
   }   
}	   

Organism GenRandomOrganism()
{
   Organism O;
   int i,j,len;
   NetworkPreimage P;
   
   len=16+rand()%17;
   
   O.x=rand()%(BOXSIZE); O.y=rand()%(BOXSIZE);
   
   O.Genome=AllocString(len);
   for (i=0;i<len;i++)
     O.Genome.Str[i]='0'+rand()%NUMALPHA;
   
   for (i=0;i<NUMALPHA;i++)
     O.M.templateid[i]=rand()%NUMTEMPL;
   
   O.idx=oidx; oidx++;
   P=GenomeToPreimage(O.Genome,O.M);
//   O.Net=NetFromPreimage(P);
   O.P=P;
   
   return O;
}
	
int Reproduce(Organism Base)
{
   Organism O;
   int i,j,k,l,len;
   NetworkPreimage P;
   int ofs;
   String SubStr;
   int x,y,z;
   int xm,ym;
   int count=0;
   int pick;
   int occ;

   O.Genome.len=Base.Genome.len;
   len=O.Genome.len;
   
/*   xm=-1; ym=-1;
   for (z=0;z<NETHEIGHT;z++)
     for (y=0;y<NETWIDTH;y++)
       for (x=0;x<NETWIDTH;x++)
       {
	  if (Base.P.Array[GetIdx(x,y,z)].type)
	    count++;
       }
   
   if (count<=0) return 0;

   pick=rand()%count;
   count=0;
   for (z=0;(z<NETHEIGHT)&&(count<=pick);z++)
     for (y=0;(y<NETWIDTH)&&(count<=pick);y++)
       for (x=0;(x<NETWIDTH)&&(count<=pick);x++)
       {
	  if (Base.P.Array[GetIdx(x,y,z)].type)
	    count++;
	  if (count>pick) 
	  {
	     xm=x-NETWIDTH/2+Base.x; ym=y-NETWIDTH/2+Base.y;
	  }
       }
   
   if ((xm<0)||(ym<0)||(xm>=BOXSIZE)||(ym>=BOXSIZE))
     return 0;
      
   for (i=0;i<NumOrg;i++)
   {
      x=xm-Org[i].x+NETWIDTH/2; y=ym-Org[i].y+NETWIDTH/2;
      if ((x>=0)&&(y>=0)&&(x<NETWIDTH)&&(y<NETWIDTH))
        if (Org[i].P.Array[GetIdx(x,y,NETWIDTH-1)].type>0)
           return 0;
   }
*/
   
   xm=1-NETWIDTH+rand()%(BOXSIZE+NETWIDTH);
   ym=1-NETWIDTH+rand()%(BOXSIZE+NETWIDTH);
   
   O.x=xm; O.y=ym;
   
   O.Genome=AllocString(len);
   
   for (i=0;i<len;i++)
   {
      O.Genome.Str[i]=Base.Genome.Str[i];
      if (rand()%1000001<1000000*MUTATION)
      {
	 O.Genome.Str[i]=rand()%NUMALPHA+'0';
      }
   }
   
   if ((rand()%1000001<1000000*DELETION)&&(len>8))
   {
      i=rand()%len; k=rand()%(len-8); if (k+i>=len) k=len-i-1;
      
      for (j=i;j<len-k;j++)
      {
	 O.Genome.Str[j]=O.Genome.Str[j+k];
      }
      O.Genome.len-=k;
      len-=k;
   }      
   
   if (rand()%1000001<1000000*DUPLICATION)
   {
      j=rand()%len;
      k=j+rand()%len;

      SubStr=AllocString(k-j);
      
      for (i=j;i<k;i++)
      {
         SubStr.Str[i-j]=O.Genome.Str[i%len];
      }

      ofs=rand()%len;

      InsertString(&O.Genome,SubStr,ofs);
      FreeString(SubStr);
      len=O.Genome.len;
   }
   
   P=GenomeToPreimage(O.Genome,O.M);
   O.P=P;
   
   O.active=0; O.fitness=0;
   O.idx=oidx; oidx++; 
   NumOrg++;
   Org=(Organism*)realloc(Org,sizeof(Organism)*NumOrg);
   Org[NumOrg-1]=O;
   
   return 1;
}

void FillSpace()
{
   int i,j,x,y,z,xm,ym,zm;
   
   memset(Box,0,sizeof(int)*BOXSIZE*BOXSIZE*NETHEIGHT);
   
   for (i=0;i<NumOrg;i++)
   {
      Org[i].cells=0; Org[i].solar=0; Org[i].fitness=0;
      Org[i].overlap=0;
      
      for (z=0;z<NETHEIGHT;z++)
	for (y=0;y<NETWIDTH;y++)
	  for (x=0;x<NETWIDTH;x++)
	  {
	     xm=x+Org[i].x;
	     ym=y+Org[i].y;
	     zm=z;
	     
	     if (Org[i].P.Array[GetIdx(x,y,z)].type!=0)
	     {
		Org[i].cells++;
	     
		if ((xm<BOXSIZE)&&(ym<BOXSIZE)&&(xm>=0)&&(ym>=0))
		{
		   if (Box[xm+ym*BOXSIZE+zm*BOXSIZE*BOXSIZE]==0)
		     Box[xm+ym*BOXSIZE+zm*BOXSIZE*BOXSIZE]=i+1;
		   else
		   {
		      Org[i].overlap++;
		      j=Box[xm+ym*BOXSIZE+zm*BOXSIZE*BOXSIZE];
		      if (j>0)
			Org[j-1].overlap++;
		      Box[xm+ym*BOXSIZE+zm*BOXSIZE*BOXSIZE]=-1;
		   }		
		}	     
	     }      
	  }      
   }
   
   for (x=0;x<BOXSIZE;x++)
     for (y=0;y<BOXSIZE;y++)
       Light[x+y*BOXSIZE]=1;
   
   for (z=0;z<NETHEIGHT;z++)
   {
      for (y=0;y<BOXSIZE;y++)
	for (x=0;x<BOXSIZE;x++)
	{
	   if (Light[x+y*BOXSIZE])
	   {	      
	      i=Box[x+y*BOXSIZE+z*BOXSIZE*BOXSIZE];
	      	      
	      if (i>0)
	      {
		 Org[i-1].solar++;
	      }
	      if (i!=0)
		Light[x+y*BOXSIZE]=0;	            
	   }	   
	}      
   }   
   
   for (i=0;i<NumOrg;i++)
     Org[i].fitness=ABSORB*Org[i].solar/(double)Org[i].cells-PENALTY;
}

// For this we pick another plant at random and superimpose them
// Then we see how much of each color sunlight reaches this plant's receptors
double CalcFitness(Organism *O, int idx)
{
   int x,y,z,xm,ym,zm;
   int cellcount=0;
   int r,minr;
   double fitness=0;
   int idx2;
   int dlight;
   int overlap=0;
   int abslight;
   int iter;
   int NUMITER=5;
   int ofs;
   int i,j,k;
   double dist;
   int light;
   int solar=0;
   
   memset(SolarGrid,0,NETWIDTH*NETWIDTH*NETHEIGHT*sizeof(int));
   memset(AbsGrid,0,NETWIDTH*NETWIDTH*NETHEIGHT*sizeof(int));
   
   O->cells=0;
   for (i=0;i<NumOrg;i++)
   {
      if (i!=idx)
      {
	 dist=fabs(O->x-Org[i].x);
	 if (fabs(O->y-Org[i].y)>dist) dist=fabs(O->y-Org[i].y);
	 
	 if (dist<NETWIDTH)
	 {
	    for (z=0;z<NETHEIGHT;z++)
	      for (y=0;y<NETWIDTH;y++)
		for (x=0;x<NETWIDTH;x++)
		{
		   xm=x+Org[i].x-O->x;
		   ym=y+Org[i].y-O->y;
		   zm=z;
		   
		   if ((xm>=0)&&(ym>=0)&&(zm>=0)&&
		       (xm<NETWIDTH)&&(ym<NETWIDTH)&&(zm<NETHEIGHT))
		   {
		      j=Org[i].P.Array[GetIdx(x,y,z)].type;
		      
		      if (j>0)
			AbsGrid[GetIdx(xm,ym,zm)]=1;
		   }		   
		}	    	 
	 }	 
      }      
   }     

   for (y=0;y<NETWIDTH;y++)
     for (x=0;x<NETWIDTH;x++)
     {
	SolarGrid[GetIdx(x,y,0)]=1;
     }
   
   abslight=0;
   for (z=0;z<NETHEIGHT;z++) 
   {
     for (y=0;y<NETWIDTH;y++)
      for (x=0;x<NETWIDTH;x++)
       {
	  light=SolarGrid[GetIdx(x,y,z)];	 
          
          dlight=0;
	  if (AbsGrid[GetIdx(x,y,z)])
	    dlight=1;
	  
	  {	     
	     
	     j=O->P.Array[GetIdx(x,y,z)].type;
	     if (j>0)
	     {	     
		cellcount++;
		//if (j==5) cellcount+=2;
		
		if (AbsGrid[GetIdx(x,y,z)])
		  overlap++;
		
		if (!dlight)
                {
                   abslight++;
		   dlight++;
                }
	     }
	  }
	  	
	  if (z<NETHEIGHT-1)
	  {
             light-=dlight; 
	     if (light<0) light=0;
	     SolarGrid[GetIdx(x,y,(z+1))]=light;
	  }	
       }
   }
   O->solar=abslight;
   O->overlap=overlap;
   O->cells=cellcount;
   
   return (abslight)*ABSORB-OVERLAP*overlap-cellcount*PENALTY;
}

double GetHeight(int idx)
{
   int i,j,k;

   for (k=0;k<NETHEIGHT;k++)
   {
      for (i=0;i<NETWIDTH;i++)
        for (j=0;j<NETWIDTH;j++) 
        {
           if (Org[idx].P.Array[GetIdx(i,j,k)].type>0) return (NETHEIGHT-k);
        }
   }
}
	   
void main(int argc, char **argv)
{
   int i,j,k,l;
   char Str[512];
   int iter;
   double fbar;
   int idx;
   int best,bestf;
   double lbar;
   FILE *f;
   double avgover,avgsolar,hbar,cellbar;
   int count;
   char *BASEDIR=argv[1]; 
   int timeidx=time(NULL);
	
   MUTATION=atof(argv[2]);
   PENALTY=atof(argv[3]);
 
   GenerateBaseTemplates(40);  
   
   NumOrg=50;
   Org=(Organism*)malloc(NumOrg*sizeof(Organism));
   SolarGrid=(int*)malloc(NETWIDTH*NETWIDTH*NETHEIGHT*sizeof(int));
   AbsGrid=(int*)malloc(NETWIDTH*NETWIDTH*NETHEIGHT*sizeof(int));

   Box=(int*)malloc(sizeof(int)*BOXSIZE*BOXSIZE*NETHEIGHT);
   Light=(char*)malloc(BOXSIZE*BOXSIZE);
   
   while (1)
   {
      timeidx=time(NULL);
      srand(timeidx);
      for (i=0;i<NumOrg;i++)
      {
         Org[i]=GenRandomOrganism();
         Org[i].active=1;
      }   
   
      iter=0;
      while (iter<1500)
      {
	 iter++;
	 
	 fbar=0; lbar=0;
	 hbar=0;
	 cellbar=0;
	 FillSpace();
	 for (i=0;i<NumOrg;i++)
	 {
	    fbar+=Org[i].fitness;
	    lbar+=Org[i].Genome.len;
	    cellbar+=Org[i].cells;
	 }      
	 fbar/=(double)NumOrg; lbar/=(double)NumOrg; cellbar/=(double)NumOrg;
	 
	 avgover=avgsolar=0; count=0;
	 for (i=0;i<NumOrg;i++)
	 {
	    count++;
	    avgover+=Org[i].overlap; avgsolar+=Org[i].solar;
	    hbar+=GetHeight(i);
	 }
	 if (count>0)
	 { 
	    avgover/=(double)count; avgsolar/=(double)count; hbar/=(double)count;
	 }
	 
	 sprintf(Str,"%s/%d.txt",BASEDIR,timeidx);
      f=fopen(Str,"a");
      fprintf(f,"%d %f %f %f %f %d %f %f\n",iter,fbar,lbar,avgover,avgsolar,NumOrg,hbar,cellbar);
      fclose(f);      

      sprintf(Str,"%s/tree%d.txt",BASEDIR,timeidx);
      f=fopen(Str,"wb");
      for (i=0;i<BOXSIZE*BOXSIZE*NETHEIGHT;i++)
      {
         if (Box[i]!=0)
	   fprintf(f,"%d %d %d\n",i%BOXSIZE,(i/BOXSIZE)%BOXSIZE,(i/(BOXSIZE*BOXSIZE))%NETHEIGHT);
      }
      fclose(f);
      
      for (j=0;j<NumOrg;j++)
      {
	 if ((Org[j].fitness<=0)||(rand()%10==0)) Org[j].active=-1;
/*	 
	 do
	 {
	    i=rand()%NumOrg;
	    j=rand()%NumOrg;
	 } while (i==j);
	 
         if (Org[i].fitness<Org[j].fitness)
	 {
	    Org[i].active=-1;
	 }
	 else Org[j].active=-1;*/
      }
      
      for (i=0;i<NumOrg;i++)
      {
	 if (Org[i].active==-1)
	 {
            k=0;
	    do
	    {	       
	       l=0;
	       do
	       {
		  l++;
		  j=rand()%NumOrg;	       
	       } while ((Org[j].active!=1)&&(l<2*NumOrg));
               k++;
	    } while (!Reproduce(Org[j])&&(k<1000));
	    
	    FreeOrganism(&Org[i]);
	    for (j=i;j<NumOrg-1;j++)
	      Org[j]=Org[j+1];
	    NumOrg--;

            if (k==1000)
            {
               NumOrg++;
               Org=(Organism*)realloc(Org,sizeof(Organism)*NumOrg);
               Org[NumOrg-1]=GenRandomOrganism();
            }
	 }	 
      }

	 for (i=0;i<NumOrg;i++)
	   Org[i].active=1;      
      }   
   }
}

