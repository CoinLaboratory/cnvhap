package lc1.dp.illumina;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringWriter;
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.logging.Logger;

import lc1.dp.data.collection.Info;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.swing.ColorAdapter;
import lc1.stats.Mixture2;
import lc1.stats.OrthogonalProbabilityDistribution;
import lc1.stats.ProbabilityDistribution2;
import lc1.stats.SimpleExtendedDistribution1;
import lc1.stats.TrainableNormal;
import lc1.stats.UniformDistribution;
import lc1.util.Constants;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;

public class IlluminaProbR implements Serializable {
 static final double round = 1000.0;
 
 
 public SimpleExtendedDistribution1[] mixeR;//, mixeB; 
 public SimpleExtendedDistribution1 mixeglobalR;//, mixeglobalB;

 int paramIndex = 1;

 public static void transfer(RegressParams target, 
			DoubleMatrix1D from) {
		for(int i=0; i<target.x0.size(); i++){
			target.x0.set(i, 0, from.get(i));
		}
		
		
	}
 
  Prior  initial = null;//new Prior();  // after switching
 public Prior prior = null;//new Prior(); //  this
  String[] globalRst = new String[4];
 public void calculatePriorsFromLocs( double[] pseudo , double pseudoFill, int it){
	MultipleDataRegression mdr = new MultipleDataRegression(new IlluminaProbR[] {this}, new int[] {-1},pseudo,
			pseudoFill, this.initial);
	mdr.makeRegress(true,true);
	 mdr.rst(globalRst);
	prior.priorR = copy(mdr.paramR);
	prior.priorRVar = copy(mdr.paramRVar);
 }
 
 public static DoubleMatrix1D[] copy(RegressParams[] paramBafVar){
	 if(paramBafVar==null) return null;
	DoubleMatrix1D[] res =  new DoubleMatrix1D[paramBafVar.length];
		for(int i=0; i<paramBafVar.length; i++){
			res[i] =paramBafVar[i]==null? null : paramBafVar[i].regressR.viewColumn(0);
		}
		return res;
 }
     public ProbabilityDistribution2[][] r; //this by snp position and then by alias;  
     public final ProbabilityDistribution2[]rGlobal;
   //    protected int[][] r_alias;  //first is bottom noCopies, second is top noCopies
       public  ColorAdapter ca_r;
     public  boolean[] r_train;
     
      public void setMinMax(double minR, double maxR, double minB, double maxB){
    	  for(int i=0; i<this.r.length; i++){
    		  if(r[i]!=null){
    		  for(int j=0; j<r[i].length; j++){
    			  if(r[i][j]!=null){
    				  ProbabilityDistribution2 di = r[i][j];
    				  di.setMinMax(minR,maxR, minB, maxB);
    				 
    			  }
    		  }
    		  }
    	  }
      }
      public void setMinMax(double minR, double maxR, double minB, double maxB, int i){
    		  if(r!=null && i<r.length && r[i]!=null){
    		  for(int j=0; j<r[i].length; j++){
    			  if(r[i][j]!=null){
    				  ProbabilityDistribution2 di = r[i][j];
    				  di.setMinMax(minR,maxR, minB, maxB);
    				 
    			  }
    		  }
    		  }
      }
  
     
      
     
     

	//public  Boolean[] probeOnly;
        
	
	final int[][] aliasNB;  //first number second b

	 BasisFunction basis_mean;
	BasisFunction basis_var; 
	
	final int index;
//public static boolean startWithZeroMix = false;

	public void addCollection(IlluminaProbR probRB) {
		this.r = (ProbabilityDistribution2[][]) Constants.join(this.r, probRB.r).toArray(new ProbabilityDistribution2[0][]);
		this.mixeR = (SimpleExtendedDistribution1[]) Constants.join(this.mixeR, probRB.mixeR).toArray(new SimpleExtendedDistribution1[0]);
	}
	
	
	public void reverse(){
		Constants.reverse(this.r);
		Constants.reverse(this.mixeR);
		
	}
	
public IlluminaProbR(EmissionStateSpace emstsp,
               boolean train, int ik, int noSNPS,  File clusterFile){
	 index = ik;
		double[] mixt = //startWithZeroMix ? new double[] {1,0} :
    	   new double[] {1-Constants.mixCoeff(ik),Constants.mixCoeff(ik)};
		
    	//   zeroMix = startWithZeroMix;
		this.all_alias = new int[emstsp.genoListSize()];
		if(Constants.trainIndiv()){
			 this.all_indices = new int[emstsp.genoListSize()][];
			 for(int k=0; k<all_indices.length; k++){
				 all_indices[k] = new int[]{k};
				 all_alias[k] = k;
			 }
		}else{ 
			int[][] cnToSplit = Constants.cnToSplit(ik);
		  this.all_indices = new int[cnToSplit.length][];
		  for(int k=0; k<all_indices.length; k++){
            all_indices[k] = emstsp.getGenoForCopyNo(cnToSplit[k]);
            for(int j=0; j<all_indices[k].length; j++){
            	all_alias[all_indices[k][j]] = k;
            }
		  }
         //   all_indices[1] =emstsp.getGenoExcludingCopyNo(0);

            
		}
        
		  this.emstsp = emstsp;
          makeBases();
            this.makePriors(ik, all_indices,emstsp, this.basis_mean.length());
           
    	 this.aliasNB = emstsp.aliasNB();
    	
    	   mixeR = new SimpleExtendedDistribution1[noSNPS];
    	//   this.probeOnly = probeOnly;
    		   this.mixeglobalR = new SimpleExtendedDistribution1(mixt, Double.POSITIVE_INFINITY);
          
       
         /*  this.b_alias = new int[emstsp.size()];
           this.bfrac_alias = new int[emstsp.genoListSize()];
           this.covar_alias = new int[emstsp.genoListSize()];
           this.reg_alias = new int[emstsp.genoListSize()];
        
           SortedSet<Double> bafs = new TreeSet<Double>();
           for(int i=0; i<b_alias.length; i++){
        	    b_alias[i] = emstsp.getGenoForHaplopair(i);
           }
        */  
         
           
       //    this.zeroR = Constants.r_mean(index)[0];
        this.r = new ProbabilityDistribution2[noSNPS][];
        ProbabilityDistribution2[] r_ = makeDists();
           for(int i=0 ; i<r_.length; i++){
        	   r_[i] =getInnerDistribution(i);
        	//   double x_ = this.calcLRR((double)emstsp.getCN(i),(double)Constants.backgroundCount(ik));
        	 //  double fracB = this.calcFracB((double)emstsp.getCN(i), (double)emstsp.getBCount(i));
        	 //  double nocop = emstsp.getCN(i);
        	  // double noB = emstsp.getBCount(i);
        	   //bafs.add(fracB);
        	  // this.reg_alias[i] =  multipleRows[0] ?   emstsp.getCN(i) : 0;
        	  // if(Double.isNaN(fracB)){
        //		   covar_alias[i] = 0;
        //	   }
        	//   else{
        	//	   int k_b = indexOf(fracB, b_mean);
             //      double mean_i_b = transform(b_mean[k_b]);
              //     double meanb = transform(b_mean[k_b]);
               //    this.covar_alias[i]= 1+(   mean_i_b==0 || mean_i_b==1 ? 1:0); 
        	  // }
              
           }
           {
          //baf_list = new double[bafs.size()];
           int k=0;
          /* for(Iterator<Double> it = bafs.iterator(); it.hasNext();k++){
        	   baf_list[k] = it.next();
        	   if(k<b_mean.length && Math.abs(b_mean[k]-baf_list[k])>0.01) throw new RuntimeException("!!");
           }*
           for(int i=0 ; i<r_.length; i++){
        	   double fracB = this.calcFracB((double)emstsp.getCN(i), (double)emstsp.getBCount(i)); 
        	   this.bfrac_alias[i] =  
        		   !IlluminaProbR.multipleRows[1] ?
        				   (Double.isNaN(fracB) ? 1 : 0) : 
        				   (
        		   Constants.singleBClust() ? i :
        		   indexOf(  fracB, this.baf_list ));
           }*/
           }
           this.rGlobal = r_;
           if(clusterFile!=null && clusterFile.exists() && Constants.readClusters()){
                this.readClusterFile(clusterFile, mixt);
        	    }	
           this.mixeglobalR.setProb(mixt);
           for(int i=0; i<noSNPS; i++){
            	  mixeR[i] =  new SimpleExtendedDistribution1( mixt, Double.POSITIVE_INFINITY);
           }
          double mixc = Constants.mixCoeff(ik);
                 inner : for(int i=0; i<r.length; i++){
                	 r[i] = new ProbabilityDistribution2[r_.length];
                			
                	for(int j=0; j<r[i].length; j++){
                		r[i][j] = r_[j]==null ? null :
                			r_[j] instanceof OrthogonalProbabilityDistribution ? 
                					((OrthogonalProbabilityDistribution)r_[j]).clone(Double.POSITIVE_INFINITY, mixeR[i], 
                						 mixeR[i]):
                			   r_[j].clone(Double.POSITIVE_INFINITY, mixeR[i]);
                	}
                 }
         
                 r_train = new boolean[r_.length];
                 Arrays.fill(r_train, train);
            //    this.bfrac_indices =  transform( bfrac_alias);
            //     this.reg_indices = transform(this.reg_alias);
              
            /*  for(int i=0; i<this.covar_alias.length; i++){
            	  skew_count[covar_alias[i]]++;
              }
             for(int i=0; i<covar_indices.length; i++){
            	 covar_indices[i] = new int[skew_count[i]];
             }
             Arrays.fill(skew_count,0);
             for(int i=0; i<this.covar_alias.length; i++){
            	 covar_indices[covar_alias[i]][skew_count[covar_alias[i]]] = i;
           	  skew_count[covar_alias[i]]++;
             }*/
                
     /*   if(initial.priorR==null && (Constants.sum(Constants.numIt)>0)){
        	double[] zerod = new double[2];
        }*/
        if(initial!=null){
	        System.err.println("initial is ");
	       initial.print();
	    
        }
    }



 void makePriors(int ik, int[][] all_indices2, EmissionStateSpace emstsp, int len) {
	 this.prior = new Prior(ik,all_indices, emstsp, len);
 	this.initial =new Prior(prior);
	
}

protected void makeBases() {
	this.basis_mean = new BasisFunction(emstsp, Constants.backgroundCount(index), Constants.basisNme(index));
    this.basis_var = new BasisFunction(emstsp, Constants.backgroundCount(index), Constants.basisNme(index));
	
}

protected ProbabilityDistribution2[] makeDists() {
	  return  new ProbabilityDistribution2[emstsp.copyNumber.size()];
//  return	new ProbabilityDistribution2[emstsp.genoListSize()];
}


protected ProbabilityDistribution2 getInnerDistribution(int index
	  ) {
	int ind1 = this.all_alias[index];
	DoubleMatrix1D v_mean = (new DenseDoubleMatrix1D(this.basis_mean.getVals(index)));
	DoubleMatrix1D v_var = (new DenseDoubleMatrix1D(this.basis_var.getVals(index)));
	Comparable compa = this.emstsp.get(index);
	double mean_i_r = (v_mean.zDotProduct(this.prior.priorR[ind1]));
	double var_i_r = Math.sqrt(v_var.zDotProduct(this.prior.priorRVar[ind1]));
	return  
      		new OrthogonalProbabilityDistribution(
      				new TrainableNormal(compa+":", mean_i_r, var_i_r,var_i_r, round,1.0) ,
      			//	new SkewNormal(compa+":", mean_i_r, var_i_r,10, round,1.0) ,
      				new UniformDistribution(Constants.minB(index), Constants.maxB(index)));
}



public void readClusterFile(File clusterFile, double[] mixt){
	try{
    	BufferedReader br = new BufferedReader(new FileReader(clusterFile));
    	String st = "";
    	this.initial = new Prior(this.index, this.all_indices,emstsp, this.basis_mean.length());
    	st = br.readLine();
    	while(st!=null){
    		if(st.startsWith("mix e global")){
    			double mix = Double.parseDouble(st.split(":")[1].trim());
    			
    				mixt[0] = 1-mix;
    				mixt[1] = mix;
    		 st = br.readLine();	
    		}
    		else if(st.startsWith("prior")){
    			this.initial.parse(st);
    			st = br.readLine();
    		}
    		else if(st.equals("rGlobal")){
    			int kk=0;
    			inner: while((st = br.readLine())!=null ){
    				boolean x = st.startsWith("x");
    				boolean y = st.startsWith("y");
    				if(!(x|| y)) break inner;
    				//String[] string = st.split(":");//[1].split("_");
    				String[] str1 = st.split(":")[2].split("_")[1].split(",");
    				for(int j=0; j<str1.length && j<2; j++){
    					this.rGlobal[kk].setParam(j, x ? 0 : 1, 
    							j==1 ? Math.pow(Double.parseDouble(str1[j]), 2):
    							Double.parseDouble(str1[j]));
    				}
    				if(y) kk++;
    			}
    		}
    		else{
    			st = br.readLine();
    		}
    	}
    	this.prior = new Prior(initial);
    	//this.initial1 = new Prior(initial);
    	}catch(Exception exc){
    		exc.printStackTrace();
    	}
}


	
	


	final public int[][] all_indices; //indices are stratified
	final public int[] all_alias;
              
          
               
       
      
	

      
  public void drop(List<Integer> toDrop) {
	  this.r = (ProbabilityDistribution2[][])Constants.drop(this.r, toDrop);
	this.mixeR = (SimpleExtendedDistribution1[])Constants.drop(this.mixeR,toDrop);
		}
  
  public IlluminaProbR(IlluminaProbR probR, boolean sameMix) {
		this.ca_r = probR.ca_r;
	//	this.zeroR = probR.zeroR;
		this.all_alias = probR.all_alias;
		this.aliasNB = probR.aliasNB;
		this.all_indices = probR.all_indices;
		this.rGlobal = probR.rGlobal;
		this.prior = probR.prior;
		this.initial = probR.initial;
		this.mixeglobalR = probR.mixeglobalR;
		this.basis_mean = probR.basis_mean;
		this.basis_var = probR.basis_var;
	this.emstsp = probR.emstsp;
	this.index = probR.index;
	/*this.covar_indices = probR.covar_indices;
	this.covar_alias = probR.covar_alias;
	this.baf_list = probR.baf_list;
	this.bfrac_alias = probR.bfrac_alias;
	this.reg_alias = probR.reg_alias;
	this.b_alias = probR.b_alias;
	this.bfrac_indices = probR.bfrac_indices;
	this.reg_indices = probR.reg_indices;*/
	
		this.r_train = probR.r_train;
		//this.probeOnly = probR.probeOnly;
		if(sameMix){
			this.mixeR = probR.mixeR;
		//	this.mixeB = probR.mixeB;
		}
		else{
			this.mixeR= new SimpleExtendedDistribution1[probR.mixeR.length];
		//	this.mixeB= new SimpleExtendedDistribution1[probR.mixeB.length];
				for(int i=0; i<mixeR.length; i++){
			    	mixeR[i] = (SimpleExtendedDistribution1)probR.mixeR[i].clone();
			   // 	 if(useBmix) mixeB[i] = (SimpleExtendedDistribution1)probR.mixeB[i].clone();
				}
		}
		this.r = new ProbabilityDistribution2[probR.r.length][];
		inner: for(int i=0; i<r.length;i++){
			r[i] = probR.r[i]==null ? null : new ProbabilityDistribution2[probR.r[i].length];
			if(r[i]==null) continue inner;
			for(int j=0; j<r[i].length; j++){
				if(probR.r[i][j]!=null)
				r[i][j] = //r[i][j] instanceof Mixture2 ? 
					probR.r[i][j] instanceof OrthogonalProbabilityDistribution ? 
							((OrthogonalProbabilityDistribution)probR.r[i][j]).clone(Double.POSITIVE_INFINITY, mixeR[i], null):
						probR.r[i][j].clone(Double.POSITIVE_INFINITY,mixeR[i]);
						//probR.r[i][j].clone();
			}
		}
	}

	private Color mix(float[] col, float[] col2) {
    	   float two = (float) 2.0;
		//  return Color.getHSBColor(
			//	  (col[0]*col[2]+col2[0]*col2[2])/(col[2]+col2[2]), 
				//  (col[1]+col2[1])/two, (col[2]+col2[2])/two);
    	   return new Color(0,0,0);
	}
       
       private Color mix(int[] col, int[] col2) {
    	//   float two = (float) 2.0;
		//  return Color.getHSBColor(
			//	  (col[0]*col[2]+col2[0]*col2[2])/(col[2]+col2[2]), 
				//  (col[1]+col2[1])/two, (col[2]+col2[2])/two);
    	   return new Color(col[0]+col2[0],col[1]+col2[1],col[2]+col2[2]);
	}

	private Color getCol( double red, double green) {
    	 return new Color((int)Math.floor(255*red), (int)Math.floor(255*green),0);
    	   //green.g
		//green.getRGBColorComponents(compArray);
    	//   compArray[0] =(int) Math.floor(green.getRed()*d); 
    	 //  compArray[1] =(int) Math.floor(green.getGreen()*d); 
		//compArray[2]=(int) Math.floor(green.getBlue()*d);
		//return compArray;
		//return Color.getHSBColor(compArray[0], compArray[1],(float)( compArray[2]*d));
	}
final EmissionStateSpace emstsp;
       public void rMaximistation(double pseudo, int pos){
           boolean chnged = false;
           
         if(true ){
        	 regress(pos, pseudo);
        	 chnged = true;
        	 throw new RuntimeException("!!");
         }
         else{   
           double[] pseudoR = Constants.r_train(1, 2);
           double[] pseudoB = Constants.b_train(1, 2);
           double pseudoRho = 1e10;
          // double pseudoRho = Constants.rho_train(1)[2];
           for(int i=0; i<r[pos].length; i++){
               ProbabilityDistribution2 dist_r = r[pos][i];
               if(r_train[i] && dist_r !=null){
                   ((Mixture2) ((Mixture2)dist_r).dist[0])
               .maximise(pseudo*pseudoR[0], pseudo*pseudoR[1],pseudo* pseudoR[2],
            		   pseudo* pseudoB[0],pseudo* pseudoB[1], pseudo*pseudoB[2],pseudo* pseudoRho);
//            		   pseudoB1*meanvarskew[0], pseudoB1*meanvarskew[1], pseudoB1*meanvarskew[2] ,pseudoRho1);
                   chnged = true;
               }
                     // r[i].print(pw);
           }
         }
        	  this.mixeR[pos].transfer(0.0);
        	//  if(useBmix)   this.mixeB[pos].transfer(0.0);
           if(chnged) this.paramIndex++;
       }
       
     /*  public void rMaximistation(int kk, double pseudo, double[] r_train2,
   			double[] b_train, double rho_train, int kk2) {
   		
   		
   	}*/
       
  static  Algebra lg = new Algebra();   
    
   // final Double[] meanPriorR;
    //final Double[] meanPriorB;
    
    //final Set<Integer> toInclude;
   //final private double[] baf_list;
   
   
   
   
 //  public final DoubleMatrix2D x0R_probeOnly;
//final DoubleMatrix2D x0R, x0B;
   /*pseudoR1 is for pseudo counts tying individual clusters to regression based averages */
   public void regress(int pos, double pseudo){
	  if(this.r[pos]==null) return;
	  throw new RuntimeException("!!");
	   /*Regression reg =null;// new Regression(this, pos, pseudo);
	   int st = reg.fillMatrices();
	   reg.regressR(st);
	   if(!probeOnly[pos]) reg.regressBAF(st);
	   reg.setMidPoints();
	   reg.setVariance();*/
       }
   
   
   
   //public final double[] rhoPrior;
   //public final double[] stddevPriorR;
 //  public final double[]  stddevPriorB;
       
       public void getTasks(List tasks, final double pseudo
    		  ) {
         
              
                // if(Constants.plot() && plotBefore)  ip.plotR(ss1, mvf, jpR, null);
             for(int i1=0; i1<r.length; i1++){
            	 final int pos = i1;
            tasks.add(new Callable(){

				@Override
				public Object call() throws Exception {
					rMaximistation(pseudo,
							pos);
					return null;
				}
            	
            });
             }
           
       }
       
      /* public double[] getProbOverRDists(StateDistribution emC) {
           double[] prob_r = new double[this.e];
           Arrays.fill(prob_r, 0.0);
           for(int j0=0; j0<stateIndices.size(); j0++){
               List<Integer> ind =stateIndices.get(j0);
               for(int j1 =0; j1<ind.size(); j1++){
                   prob_r[j0]+=emC.dist[ind.get(j1)];
               }
           }
           return prob_r;
       }*/

       public void initialiseRCounts(){
    		  this.mixeglobalR.initialise();
    		  for(int j=0; j<rGlobal.length; j++){
    			  rGlobal[j].initialise();
    		  }
    	  
           for(int i=0; i<r.length; i++){
        	   if(mixeR[i]!=null)this.mixeR[i].initialise();
        	   if(r[i]!=null){
        	   for(int j=0; j<r[i].length; j++){
        		   ProbabilityDistribution2 dist_r = r[i][j];
        		   dist_r.initialise();
        	   }
        	   }
           }
        }
       public void addRCount( int obj_j, double weight, double r_i, double b_i, int i){
    	   if(true) throw new RuntimeException("!!");
       
       }
      
      
       public void print(int i) {
    	StringWriter sw = new StringWriter();
    	PrintWriter pw = new PrintWriter(sw);
    	pw.print(String.format("%5.3g ", mixeR[i].probs[1]));
 //   	 if(useBmix) pw.print(String.format("%5.3g ", mixeB[i].probs[1]));
    	if(this.r[i]!=null){
    		for(int k=0; k<r[i].length; k++){
    			r[i][k].print(pw);
    		
    			pw.println();
    		}
    	}
    	pw.close();
    	System.err.println(sw.toString());
   	}	
       
       public void print(PrintWriter pw) {
    	   pw.print("mix e global:\t");
    	   if(mixeglobalR!=null){
    		   pw.println(mixeglobalR.probs[1]);
    	   pw.print("mix e global probe Only:\t");
    	   pw.println("glob regression ");
    	   this.prior.print(pw);
    	   pw.println("rGlobal");
    	   for(int i=0; i<this.rGlobal.length; i++){
    		   rGlobal[i].print(pw);
    		   pw.println();
    	   }
    	   }
    	   pw.print("ProbR:\t");
    	   for(int i=0; i<this.mixeR.length; i++){
    		   if(mixeR[i]==null) pw.print("null");
    		   else{
    			   pw.print(String.format("%5.3g ", mixeR[i].probs[1]));
    		   }
    	   }
    	   pw.println();
    	   pw.println();
        }
       
       public double calcR(int obj_j, double r_i, double b_i, int i){
          double res =  r[i][obj_j].probability(r_i, b_i);
          return res;
       }
      
       
       public ProbabilityDistribution2 getDistribution(int noCop, int noB,
   			int pos) {
    	   int al = this.aliasNB[noCop][noB];// b_alias[j];
         if(r[pos]==null  || al >=r[pos].length) {
        	 return null;
         }
         ProbabilityDistribution2 dist =    r[pos][al];
          return dist;
   	}
       
       public ProbabilityDistribution2 getDistributionGlobal(int noCop, int noB
      			) {
    	   if(noB>0) return null;
       	   int al = this.aliasNB[noCop][noB];// b_alias[j];
       	ProbabilityDistribution2 dist =  al<rGlobal.length ? 
           	  rGlobal[al] : null;
             return dist;

      	}
     








//public final double zeroR;

    public int getParamIndex() {
       return paramIndex;
    }
    
   
   /* private Color modify(Color red) {
       return new Color(Math.min(255, red.getRed()+50), 
               Math.min(255, red.getGreen()+50),
               Math.min(255, red.getBlue()+50));
        
    }*/
	public void setToExclude(int pos) {
		for(int k=0; k<this.r[pos].length; k++){
			r[pos][k].setToExclude();
		}
		
	}
	public void append(IlluminaProbR probRB) {
		this.mixeR =(SimpleExtendedDistribution1[]) Constants.append(mixeR, probRB.mixeR, SimpleExtendedDistribution1.class);
		this.r_train = 
			(boolean[]) Constants.append(this.r_train, probRB.r_train, boolean.class);
		this.r = (Mixture2[][]) Constants.append(this.r, probRB.r, Mixture2[].class);
		//this.probeOnly = (Boolean[]) Constants.append(this.probeOnly, probRB.probeOnly, Boolean.class);
	}	
	
	public static IlluminaProbR getMerged(IlluminaProbR[]dists, Info[][] info) throws Exception{
		dists[0].mixeR =(SimpleExtendedDistribution1[]) append(dists,"mixeR", SimpleExtendedDistribution1.class, info);
		dists[0].r = (Mixture2[][]) append(dists,"r", Mixture2[].class, info);
		//dists[0].probeOnly = (Boolean[]) append(dists, "probeOnly", Boolean.class, info);
		return dists[0];
	}
	
	public static Object append(IlluminaProbR[] dists, String type, Class componentType, Info[][] map) throws Exception{
		Object[] left = new Object[dists.length];
		for(int i=0; i<left.length; i++){
			left[i] = IlluminaProbR.class.getField(type).get(dists[i]);
		}
		//	Class componentType = Array.get(left, 0).getClass();
		int[] len = new int[left.length];
		int length = 0;
		for(int i=0; i<len.length; i++){
			len[i] = Array.getLength(left[i]);
			length+=len[i];
		}
			Object res = Array.newInstance(componentType, length);
			outer: for(int i=0; i<map.length; i++){
				for(int j=0; j<map[i].length; j++){
					if(map[i][j]!=null){
						Info inf =map[i][j];
						Object value = Array.get(left[inf.data_index], inf.relative_position);
						Array.set(res, i, value);
						continue outer;
					}
				}
				throw new RuntimeException("!!");
			}
		
			return res;
		}
	public void transferMix(double mp) {
		for(int kk=0; kk< this.mixeR.length; kk++){
			if(mixeR[kk]!=null) {
				mixeR[kk].transfer(mp);
			}
		}
		
	}
	
	
	
	
	
	//public static int change = Constants.changePriorCount();
	
	
	

	//final double[] rVar, rVarPrior;
	//final DoubleMatrix1D[] bafVar, bafVarPrior;  //indexed by bfrac_indices - i.e for diff baf values
	public void addCount(int noCop, int noB, int i, Double r2, Double b,
			double val) {
		if(Constants.allowLocalDist()){
			ProbabilityDistribution2 dist =  this.getDistribution(noCop,noB, i);
		    dist.addCount(r2, b, val);
		}
		if(Constants.measureGlobal()){
				ProbabilityDistribution2 dist =  this.getDistributionGlobal(noCop, noB);
				dist.addCount(r2, b, val);
	}
	}
	public void maximiseGlobal(double[] pseudo,  int it) {
		if(Constants.measureGlobal() && Constants.min(pseudo)<1e10 ){
		this.mixeglobalR.transfer(pseudo[2]);
			/*if(!Constants.maximiseGlobalRegress()){	
				double[] pseudo1 = null;//Constants.r_train2G[0];
				double[] pseudo2 = null;//Constants.b_train2G[0];
			 for(int k=0; k<this.rGlobal.length; k++){
				
				  rGlobal[k].maximise(pseudo1[0], pseudo1[1], pseudo1[2],
		    		        pseudo2[0], pseudo2[1], pseudo2[1],pseudo2[2]);
			  }
			  ///WARNING - THIS IS NOT UPDATING THE PRIOR. JUST CLUSTER POSITIONS
			//  this.calculatePriorsFromLocs(true,new double[] {0.0, 0.0},1.0, -1);
			}
			else{*/
			  this.calculatePriorsFromLocs( pseudo, 0.0, it );   
			//}
			for(int i=0; i<mixeR.length; i++){
				if(mixeR[i]!=null){
					this.mixeR[i].setPriors(this.mixeglobalR);
				}
			}
		}
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		pw.println("initial");
		this.initial.print(pw);
		pw.println("prior");
		this.prior.print(pw);
		Logger.global.info("update global\n"+sw.getBuffer().toString());
	}

	public boolean canMergeWith(IlluminaProbR probRB) {
		return this.basis_mean.equals(probRB.basis_mean) && this.prior.equals(probRB.prior);
	}

	
	
	
}



/*public  ProbabilityDistribution2 makeMixture(ProbabilityDistribution2 prob, SimpleExtendedDistribution1 mixeR, 
SimpleExtendedDistribution1 mixeB, boolean orthogMix, int index){
if(orthogMix && prob instanceof OrthogonalProbabilityDistribution){
OrthogonalProbabilityDistribution opd = (OrthogonalProbabilityDistribution) prob;
opd.distx = makeMixture(opd.distx, mixeR,  Constants.minR(index), Constants.maxR(index));
opd.disty = makeMixture(opd.disty, useBmix ? mixeB : mixeR, min(),max());
return opd;
}
else{
return    new Mixture2( prob, Constants.minR(index), Constants.maxR(index), min(),max(), mixeR);
}


}

public static ProbabilityDistribution makeMixture(ProbabilityDistribution prob, SimpleExtendedDistribution1 mixeR, double min, double max){
if(prob instanceof UniformDistribution) return prob;
else if(prob instanceof Mixture){
throw new RuntimeException("!!");
}
return new Mixture(	prob,min, max,mixeR);
}*/





/*  private DoubleMatrix1D getBafVar(double[][] b_var, List<Integer> nocops,  double[]  bfrac) {
	  int len = nocops.size() * bfrac.length;
	  int[] rows = new int[len];
	  DoubleMatrix2D mat = new DenseDoubleMatrix2D(len,5);
	 
    	DoubleMatrix2D Y = new DenseDoubleMatrix2D(len,1);
    	int k=0;
    	for(int i=0; i<nocops.size(); i++){
    		 int noCop = nocops.get(i);
    		 double ncop = trans(noCop);
    	for(int i1=0; i1<bfrac.length; i1++){
    		
    		 
    		
    		double bcopt = bfrac[i1];
    	//	double mult = 0.0;
    		
    		double target=Math.pow(b_var[i][i1],2);
    		
    		
    		
    		
           double bcopt1 = trans1(bcopt,true);
	    		mat.setQuick(k, 0, 1.0);
	    		mat.setQuick(k, 1, ncop);
	    		mat.setQuick(k, 2, Math.pow(ncop,2));
	    		mat.setQuick(k, 3, bcopt);
	    		mat.setQuick(k, 4, bcopt1);
	    		Y.setQuick(k, 0, target);
	    		rows[k] = k;
	    		k++;
    	}
    	}
    	DoubleMatrix2D mat1 = mat;//mat.viewSelection(rows, new int[] {0,1,2});
    	 DoubleMatrix2D Q = new DenseDoubleMatrix2D(mat1.columns(), mat1.columns());
    	 DoubleMatrix2D x0 = new DenseDoubleMatrix2D( mat1.columns(),1);
    	 for(int i=0; i<mat1.columns(); i++){
    		 Q.set(i, i, 1e-3);
    		 x0.set(i, 0, 0);
    	 }
    //	 Q.set(2, 2, 2);
    
    	DoubleMatrix2D res =  RegressParams.solve(mat1, Y, Q, x0);
		
			return res.viewColumn(0);
	}*/


//  final int[] reg_alias;
 // final int[] covar_alias;
 // final int[][] covar_indices;
 // final int[][] nonskew_index;
      //  public final int[] b_alias;
//final public int[] bfrac_alias;
 
/* public int bfrac_alias(int i){
		 return this.bfrac_alias[i];
 }
 public int[] bfrac_indices(int i) {
		 return this.bfrac_indices[i];
	}*/
 

/*
private double avg(double d, double e, double f, double nocop) {
//  double st =  Math.pow(d, 2);
  double st1 =  Math.pow(e, 2);
  double end = Math.pow(f, 2);
  if(nocop==1) return d;
  if(nocop==2) return e;
 
  double res = 
	
	  st1 + ((end - st1)/(Math.log(4) -Math.log(2) ))*(Math.log(nocop)-Math.log(2));
  double res1 = Math.sqrt(Math.max(0.0001,res));
  if(Double.isNaN(res1)){
	  throw new RuntimeException("!!");
  }
  return res1;

}
private boolean skew(double skew_i_b) {
	// TODO Auto-generated method stub
	return false;//skew_i_b!=0;
	//return false;
}*/






/*final DoubleMatrix1D rVarM, rVarMPrior;
 public DoubleMatrix1D getAssumedLRRVar(double[] rVar, double[] r_ratio, double bg){
	  int len = rVar.length-1;
    	DoubleMatrix2D mat = new DenseDoubleMatrix2D(len,3);
    	DoubleMatrix2D Y = new DenseDoubleMatrix2D(len,1);
    	
    	for(int i=1; i<len+1; i++){
    		int k = i-1;
    		double v = trans(r_ratio[i]*bg);
	    		mat.setQuick(k, 0, 1.0);
	    		mat.setQuick(k, 1, v);
	    		mat.setQuick(k, 2, Math.pow(v,2));
	    		Y.setQuick(k, 0, rVar[i]);
    	}
    	 DoubleMatrix2D Q = new DenseDoubleMatrix2D(mat.columns(), mat.columns());
    	 DoubleMatrix2D x0 = new DenseDoubleMatrix2D( mat.columns(),1);
    	 for(int i=0; i<mat.columns(); i++){
    		 Q.set(i, i, 1e-3);
    		 x0.set(i, 0, 0);
    	 }
    	 Q.set(2, 2, 1);
    
    	DoubleMatrix2D res =  RegressParams.solve(mat, Y, Q, x0);
    	return res.viewColumn(0);
    }*/


/* double backgroundCount() {
return Constants.backgroundCount();
}*/
/** Assumes is ordered from smallest to largest 
* returns greater than or equals
* 
private static int indexOf(double x_, double[] x) {

for(int i=0; i<x.length; i++){
if(x_ == x[i] || (Double.isNaN(x_) && Double.isNaN(x[i]))) return i;
if(i>0 && x[i] < x[i-1]) {
	throw new RuntimeException("should be sorted in ascending order");
}
if(x_<=x[i]) return i;
}
return -1;
}
private int indexOf(int x_, int[] x) {
for(int i=0; i<x.length; i++){
	if(x_==x[i]) return i;
}
return -1;
}


//final   public Shape[] shape;

// public int[][]bfrac_indices;
//int[][] reg_indices;*/

/*    private static  int[][] transform( int[] bfrac_alias2) {
SortedSet<Integer> vals = new TreeSet<Integer>();
for(int i=0; i<bfrac_alias2.length; i++){
	   vals.add(bfrac_alias2[i]);
}


int[][]res =  new int[vals.size()][];
int[] cnt = new int[vals.size()];
for(int i=0; i<bfrac_alias2.length; i++){
	int al = vals.headSet(bfrac_alias2[i]).size();
	cnt[al]++;
}
for(int i=0; i<vals.size(); i++){
	res[i] = new int[cnt[i]];
}
Arrays.fill(cnt,0);
for(int i=0; i<bfrac_alias2.length; i++){
	int al = vals.headSet(bfrac_alias2[i]).size();
	
	res[al][cnt[al]] = i;
	cnt[al]++;

}
return res;
}*/

/*public void initColors() {
if(this.ca_r!=null) return;
SortedMap<Double, String> m = new TreeMap<Double, String>();
List<Color> l=  new ArrayList<Color>();

for(int i=0; i<this.r[0].length; i++){
   SkewNormal sn =  (SkewNormal)r[0][i];
   m.put(sn.getMean(), sn.name());
}
List<Double> ml = new ArrayList<Double>(m.keySet());
double mid = ((double) ml.size()-1)/2.0;
int floor = (int) Math.floor(mid);
double avg =0;
if(ml.size()<=1) avg =0;
else if(Math.abs(floor-mid)<0.001){
	avg=ml.get(floor);
}
else{
	avg=(ml.get(floor)+ml.get(floor+1))/2.0;
}
Color black, green, red;
black = green=red = null;



for(Iterator<Double> it = m.keySet().iterator(); it.hasNext();){
    Double key = it.next();
    if(key <avg-0.05){
        if(red!=null){
            red = modify(red);
           
        }
        else red = Color.RED;
        l.add(red);
    }
    else if(key>avg+0.05){
        if(green!=null){
            green = modify(green);
           
        }
        else green = Color.GREEN;
        l.add(green);
    }
    else{
        if(black!=null){
            black = modify(black);
           
        }
        else black = Color.black;
        l.add(black);
    }
}
this.ca_r = new ColorAdapter(l.toArray(new Color[0]));
for(Iterator<String> it = m.values().iterator(); it.hasNext();){
    ca_r.getColor(it.next());
}
}*/