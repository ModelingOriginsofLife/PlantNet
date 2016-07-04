#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct
{
   int x,y;
   
   double z[8],z2[8];
   int count;
} Entry;

int NEnt=0;
Entry *Ent=NULL;

int FindEnt(int x)
{
   int i;
   
   for (i=0;i<NEnt;i++)
     if ((Ent[i].x==x))
       return i;
   
   return -1;   
}

void BlendNeighbors()
{
   int i,j,k;
   int minj;
   int mindist;
   
   for (k=0;k<NEnt;k++)
   {
      i=rand()%NEnt;
      minj=0;
      mindist=1e8;
      for (j=0;j<NEnt;j++)
      {
	 if ((j!=i)&&(abs(Ent[j].x-Ent[j].y)<mindist))
	 {	    
	    mindist=abs(Ent[j].x-Ent[j].y);
	    minj=j;
         }	 
      }    
      for (j=0;j<8;j++)
      {	   
	 Ent[i].z[j]=(Ent[i].z[j]+Ent[minj].z[j]);
	 Ent[i].z2[j]=Ent[i].z2[j]+Ent[minj].z2[j];
      }      
      Ent[i].count+=Ent[minj].count;
      
      for (j=minj;j<NEnt-1;j++)
	Ent[j]=Ent[j+1];
      NEnt--;
   }  
}

void main(int argc, char **argv)
{
   FILE *f;
   int i,j;
   int a,b;
   double c[8];
   f=fopen(argv[1],"rb");
   
   while (fscanf(f,"%d",&a)!=EOF)
   {
      for (j=0;j<7;j++)
	fscanf(f,"%lf",&c[j]);
	
//      if (a%10==0)
      {	   
	 i=FindEnt(a);
	 
	 if (i==-1)
	 {
	    NEnt++;
	    Ent=(Entry*)realloc(Ent,sizeof(Entry)*NEnt);
	    i=NEnt-1;
	    Ent[i].count=0;
	    Ent[i].x=a; 
	    for (j=0;j<8;j++)
	    {
	       Ent[i].z[j]=0; Ent[i].z2[j]=0;
	    }	 
	 }      
	 
	 Ent[i].count++;
	 for (j=0;j<8;j++)
	 {
	    Ent[i].z[j]+=c[j]; Ent[i].z2[j]+=(double)c[j]*c[j];
	 }      
      }      
   }   
   
   fclose(f);
   
/*   for (i=0;i<100;i++)
   {
      printf("Blending %d\n",i);
      BlendNeighbors();
   }*/   
   
   f=fopen(argv[2],"wb");
   
   for (i=0;i<NEnt;i++)
   {
      fprintf(f,"%d ", Ent[i].x);
      for (j=0;j<8;j++)
      {
	 fprintf(f,"%f %f ", Ent[i].z[j]/(double)Ent[i].count,
		 
	     	sqrt(Ent[i].z2[j]/(double)Ent[i].count-
		  pow(Ent[i].z[j]/(double)Ent[i].count,2))/sqrt(Ent[i].count)
	    );
      }
      
      fprintf(f,"\n");
   }
      
   fclose(f);
}
