/* 
   Function to process genotype data from snpMatrix package, stored as raw R data
*/

#include <stdio.h>
#include <stdlib.h>

void writeGenoDat(unsigned char *geno, int *numprobes, int *numsub, unsigned char *lev, int *numlev, int *sumallele, int *caco, int *levcaco, int *numlevcaco, double *outInv, double *outAleleSum, double *outNoMissCount)
{
   int i,j,l,k,dat,nomiss[*numlevcaco],sum[*numlevcaco];
   double allele1sum,allele2sum,invallele;
   
  	
	
   for(j=0; j< *numprobes; j++)
   {
   
      /*initialize sums for each probe*/
      for(k=0; k<*numlevcaco; k++)
      {
         sum[k]=0;
		 nomiss[k]=0;
      }
	
	
	  allele1sum=0.0;
	  allele2sum=0.0;
	  
	  /*identify minor and mayor alleles*/
	  for (i=0; i<*numsub; i++)
	  {  
		 if (geno[i+j*(*numsub)]==lev[1]) 
		 { 
		    allele1sum=allele1sum+1; 
		 }
		 if (geno[i+j*(*numsub)]==lev[3])
		 {
		     allele2sum=allele2sum+1;
		  }	
	  } 
	  
	  invallele=allele1sum-allele2sum;	
	  /*run over subjects*/
	  for (i=0; i<*numsub; i++)
	  {  
	
	    
    	 for (l=0; l<*numlev; l++)  /* run over numlevels to check which one the data point belongs to and then assign the corresponding integer*/
	     { 
           if (geno[i+j*(*numsub)]==lev[l])
           {
		     dat=l;
			}
         } 
		 
		 outInv[j]=0;
		 if(invallele<0)
		 { 
		   outInv[j]=1;
		   
		   switch(dat){
		   case 1: 
		     dat=3;
			 break;
		   case 3:
		     dat=1;
			 break;
		   }
	     }
		 
		 
		 /*compute Allele Sums by subject level*/
		 if(*sumallele==1)
         {
		   for(k=0; k<*numlevcaco; k++)/* run over cacolevels to check which subject the data point belongs and add it's allele (dat-1) to the subject group*/
		   {
	         if(dat>0) /*discard missings encoded by 0*/
		     { 
						
		       if (caco[i]==levcaco[k])
			   {
                             sum[k]=sum[k]+(dat-1);	/*sum alleles by level*/
			     nomiss[k]=nomiss[k]+1;     /*count not missing data*/
			    }
		      }
		    }
		  } 
		
		
	   }/*end subject run*/
	   
		
	   /*write AlleleSumDat and notMissDat*/
	   if(*sumallele==1)
	   {
	     for(k=0; k<*numlevcaco; k++)
	     {
		    outAleleSum[k+j*(*numlevcaco)]=sum[k];
			outNoMissCount[k+j*(*numlevcaco)]= nomiss[k];
		  }	

	   }
	   
	  
    }

  
   return;

}

  
 
