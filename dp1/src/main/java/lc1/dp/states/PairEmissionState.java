package lc1.dp.states;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import lc1.dp.core.DoublePool;
import lc1.dp.data.representation.Emiss;
import lc1.dp.emissionspace.CompoundEmissionStateSpace;
import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.stats.PseudoDistribution;
import lc1.util.Constants;

public  class PairEmissionState extends CompoundState {

    private final EmissionState[] dist;
    
    public  void append(EmissionState emissionState){
       for(int k=0; k<dist.length; k++){
    	   dist[k].append(((PairEmissionState)emissionState).dist[k]);
       }
         }
    protected final CompoundEmissionStateSpace emStSp;
   // final public  PairMarkovModel parentModel;
    final boolean decomp;
   /* public CompoundMarkovModel getParentModel(){
        return parentModel;
    }*/
    public Object clone() {
       PairEmissionState res =  new PairEmissionState(copy(this.dist) ,  emStSp, decomp);
       //res.calculateIndex();
       return res;
    }
    
    @Override
    public Integer getFixedInteger(int i) {
        // TODO Auto-generated method stub
        return null;
    }
    final int[] indices_tmp;
    @Override
 public Integer calculateIndex(int i){
        EmissionState[] states = getMemberStates(false);
     
     //   Integer[] result = new Integer[this.noSnps()];
     //   Integer[] internal_indices = new Integer[states.length];//[this.noSnps()];
        for(int j=0; j<states.length; j++){
            if(   states[j] instanceof CompoundState){
                Integer ind = ((CompoundState)states[j]).calculateIndex(i);
                if(ind==null) return null;
                indices_tmp[j] = ind;
               // if(internal_indices[j]==null){
               //   throw new RuntimeException("!!");
              //  }
            }
            else{
             //   for(int i=0; i<internal_indices[j].length; i++){
                    Integer ind = states[j].getFixedInteger(i);
                    if(ind==null) return null;
                    indices_tmp[j] = ind;
               // }
            }
        }
       // for(int i=0; i<this.noSnps(); i++){
         //   int[] indices = new int[states.length];
       // for(int jk=0; jk<states.length; jk++){
          
        //   Integer index_i = internal_indices[jk];
         //  if(index_i==null){
         //      return null;
         //  }
         //  indices[jk] =  index_i;
       // }
         Integer result = ((CompoundEmissionStateSpace)this.emStSp).getIndex(indices_tmp);
         if(result==null){
         /*Map<int[], Integer> m =  this.emStSp.membersIndexToIndex;
         for(Iterator<Entry<int[], Integer>> it = m.entrySet().iterator(); it.hasNext();){
             Entry<int[], Integer> nxt = it.next();
             System.err.println(nxt.getKey()[0]+" "+nxt.getKey()[1]+" "+nxt.getKey().equals(indices));
         }*/
             throw new RuntimeException("is null "+indices_tmp[0]+" "+indices_tmp[1]+" "+indices_tmp.length+" ");
         }
       // }
      //  int[] res =  emStSp.getMemberIndices(result);
       //  if(result > 30){
        //     throw new RuntimeException("!!");
         //}
        return result;
    }
    
    
    @Override
    public double[] getEmiss(int i) {
        throw new RuntimeException("should not call this!");
    }
    
    public String toString(int i){
        StringBuffer sb = new StringBuffer(dist[0].toString(i));
        for(int j=1; j<dist.length; j++){
            sb.append(dist[j].toString(i));
        }
        return sb.toString();
  }
    
    protected static List<EmissionState> copy(EmissionState[] em){
        EmissionState[] res = new EmissionState[em.length];
        for(int k=0; k<res.length; k++){
            res[k] =(EmissionState)em[k].clone();
        }
        return Arrays.asList(res);
    
    }
    
    
    
    public int mostLikely(int pos){
        int[] mostLikely = new int[dist.length];
        for(int j=0; j<mostLikely.length; j++){
            mostLikely[j] = (this.getInnerState(j, false)).mostLikely(pos);
        }
        return this.emStSp.getIndex(mostLikely);
    }
    
    static String getConcatenation(Iterator<EmissionState> states){
        StringBuffer sb = new StringBuffer();
        while(states.hasNext()){
            EmissionState nxt = states.next();
          //  if(states.get(i).adv!=states.get(0).adv) throw new RuntimeException("must have same advance");
            sb.append(nxt.getIndex());
            if(states.hasNext()) sb.append("_");
        }
        return sb.toString();
    }
    
   String name1;
   public static Integer getNoCop(List<EmissionState> dist){
       int noCop = 0;
       for(int i=0; i<dist.size(); i++){
    	   Integer nc = dist.get(i).noCop();
    	   if(nc==null) return null;
           noCop+=nc.intValue();
       }
       return noCop;
   }
   
   /*public static Integer[] getNoCop(List<EmissionState> dist){
	   Integer[] res = new Integer[Constants.inputDir.length];
	   for(int i=0; i<res.length; i++){
		   res[i] = getNoCop(dist,i);
	   }
      
       return res;
   }*/
    public String getEmissionName(){
        return name1;
    }
    final int noCopies;
    
    final Integer noCop;
   
   
  /*  public Integer noCop(int di){
        return noCop[di];
    }*/
    public Integer noCop(){
    	return noCop;
    }
    
   // static ProbRManager pm = new ProbRManager();
  
   /* @Override
public int mod(int j, int di){
	//NOTE - this will work based on ordering of emstsp for B allele counts =0, not sure about B allele counts > 0
	//USE ONLY FOR READ PAIR AT THIS STAGE
	
	return //offset==null ? j :
		    j+offset[di];
}*/
//public int[] offset =null;
    
    public PairEmissionState(List<EmissionState> st1,boolean decompose ){
    	this(st1, getEmissionStateSpace(st1,decompose), decompose);
    	//this.offset = makeOffset(st1);
    	
    	
    }
    
   /* public int[] makeOffset(List<EmissionState> st1){
    	int len = Constants.inputDir.length;
    	
    	
		this.offset=  new int[len];
		for(int i=0; i<st1.size(); i++){
    			for(int k=0; k<len; k++){
    				offset[k]+=st1.get(i).mod(0, k);
    			}
    	}
		return offset;
    }*/
    
    private static CompoundEmissionStateSpace getEmissionStateSpace(List<EmissionState> st1,
			boolean decompose) {
    	  int[] numCopies = new int[st1.size()];
          boolean useFixed = true;
           for(int i=0; i<st1.size(); i++){
        	   EmissionStateSpace emstsp = ((EmissionState)st1.get(i)).getEmissionStateSpace();
        	  numCopies[i] = emstsp.copyNumber.get(0);
        	  if(emstsp.copyNumber.size()>1){
        		  useFixed = false;
        	  }
           }
           CompoundEmissionStateSpace  emstsp = !useFixed ? Emiss.emiss.spaceByPloidy[st1.size()-1]:
         // CompoundEmissionStateSpace emstsp =//Emiss.getSpaceForNoCopies(Constants.backgroundCount());
         	Emiss.getEmissionStateSpace(numCopies);
         return emstsp;
	}

//** assume they have the same emission space */
   public  PairEmissionState(List<EmissionState> dist, CompoundEmissionStateSpace emStSp, boolean decompose){
        super(getConcatenation(dist.iterator()), 1);
        noCop = getNoCop(dist);
        if(emStSp.getMembers().length!=dist.size()) 
        	throw new RuntimeException("!!");
     //  if(dist.size()==2 && dist.get(0).noCop()+dist.get(1).noCop()==0){
    	//   Logger.global.info("h");
     //  }
        this.emStSp =  emStSp;
        this.decomp = decompose;
    //    this.parentModel = parent;
        this.noCopies = dist.size();
        this.dist = dist.toArray(new EmissionState[0]);
        SortedSet<EmissionState> distSet = new TreeSet<EmissionState>(dist);
        this.name1 = getConcatenation(distSet.iterator());
        this.indices_tmp = new int[noCopies];
    	//this.offset = makeOffset(dist);
    //this just sets the state space
      
         
        }
    
    public EmissionState getInnerState(int j, boolean real){
        return this.dist[j];
    }
    /*protected EmissionStateSpace getEmissionStateSpace(int i){
        return this.getInnerState(i, true).getEmissionStateSpace();
       
    }*/
    
    /*
    public double score(ComparableArray key, int i,  boolean recursive, boolean decompose){
        if(decompose){
            double score =0;
               Comparable[] l = this.emStSp.getHapls(key);;
                for(Iterator<Comparable> poss = Arrays.asList(l).iterator(); poss.hasNext();
                 ){
                    List li =  ((ComparableArray)poss.next()).elements();
                    Iterator<Comparable> obj =li.iterator();
                    double sc = this.getInnerState(0, false).score(getEmissionStateSpace(0).get(obj.next()), i);
                      //      getParentModel(0).getEmissionStateSpaceIndex(obj.next()),i);
                    for(int j=1; obj.hasNext(); j++){
                        sc*=this.getInnerState(j, false).score(getEmissionStateSpace(j).get(obj.next()), i);
                    }
                    score+=sc;
                }
            return score;///(double)l.size();
        }
        else{
            double sc = 1;
            for(int j=0; j<key.size(); j++){
                EmissionState innerSt = getInnerState(j, false);
                if(recursive && innerSt instanceof CompoundState){
                    sc*=((CompoundState)innerSt).score((ComparableArray)key.get(j), i,  recursive, decompose);
                }
                else{
                    sc*=innerSt.score(this.getEmissionStateSpace(j).get(key.get(j)), i);
                }
            }
            return sc;
        }
    }*/
    
   public double score(int key, int i,  boolean recursive, boolean decompose){
        if(decompose){
            
               int[] l = this.emStSp.getHaps(key);
              double sc =0;
           //   double[] w = this.emStSp.getHapWeights(key);
                for(int k=0; k<l.length;k++
                 ){
                    int li = l[k];// ((ComparableArray)poss.next()).elements(); //haplotype index
                   sc+= score(li, i, recursive, false);
                }
            return sc;///(double)l.size();
        }
        else{
            double sc = 1;
           int[] indices =  this.emStSp.getMemberIndices(key);
           int[] memInd = emStSp.getMemberTypeIndices(key);
          
            for(int j=0; j<indices.length; j++){
	                EmissionState innerSt = getInnerState(j, false);
	                EmissionStateSpace emstsp = emStSp.getMembers()[memInd[j]];
	            	if( innerSt.getEmissionStateSpace().size()==emstsp.size()){
	                if(recursive && innerSt instanceof CompoundState){
	                    sc*=((CompoundState)innerSt).score(indices[j], i,  recursive, decompose);
	                }
	                else{
	                    sc*=innerSt.score(indices[j], i);
	                }
	            	}
            //	}
	            	else{
	            		sc = 0.0;
	            	}
            	//}
            	//else sc = 0;
            }
            return sc;
        }
    }
   
  public  static DoublePool pool = new DoublePool();
   public void addCount(int key,  double value, int i, boolean decompose) {
       try{
       if(decompose){
    	  // this.score(key, i, false, decompose);
               int[] l =this.emStSp.getHaps(key); 
               if(l.length>1){
                   double[] prob =pool.getObj(l.length);
//                               new double[l.length];
                   double sum=0;
                   for( int j=0; j<prob.length; j++){
                       prob[j] = this.score(l[j], i, true, false);//*w[j];
                       sum+=prob[j];
                   }
                   if(sum==0){
                       Arrays.fill(prob, 1.0 / (double)prob.length );
                       sum = 1.0;
                   }
                   for(int j=0;j<prob.length; j++){
                       double prob_j = value* (prob[j]/sum);
                     if(prob_j>0)  this.addCount(l[j], prob_j, i, false);
                   }
                   pool.returnObj(prob);
               }
               else{
                   this.addCount(l[0], value, i,false);
               }
      }
      else{
          int[] indices =  this.emStSp.getMemberIndices(key);
          for(int j=0; j<indices.length; j++){
        	  EmissionState inner1 = getInnerState(j, false);
                inner1.addCount(indices[j], value, i);
          }
      }
       }catch(Exception exc){
           exc.printStackTrace();
           System.exit(0);
       }
  }
   public void addCountDT(double key,  int phen_index, double value, int i) {
       for(int k=0; k<this.dist.length; k++){
           dist[k].addCountDT(key,phen_index, value/dist.length, i);
       }
  }

    
   @Override
   public double scoreEmiss(Double[] object_index, int i1){
       double sc = 1.0;
       for(int k=0; k<object_index.length; k++){
           if(object_index[k]==null) continue;
           throw new RuntimeException("!! ");
//           sc*=emissionsDT[k][i1].probability(object_index[k]);
       }
       //int i= emissions.length==1 ? 0 : i1;
      return sc;
       //return probs[object_index];
   }
   
    @Override
    public void addCount(int key1, double value, int i) {
        if(value==0) return;
       // ComparableArray key = (ComparableArray)this.getEmissionStateSpace().get(key1);
        this.addCount(key1, value, i, decomp);
    }
    
   
    
    public double score(int key1, int i){
       // Comparable key =  this.getEmissionStateSpace().get(key1);
        double sc = score(key1, i,  false, decomp);
      //  System.err.println(key+"  cf "+this.dist[0].getBestIndex(i)+" "+this.dist[1].getBestIndex(i)+" "+sc);
        return sc;
    }
    
    /*public void setRandom(double u, boolean restart){
       for(int i=0; i<dist.length; i++){
           dist[i].setRandom(u, restart);
       }
    }*/
    
    public boolean transferCountsToProbs( double pseudo) {
        throw new RuntimeException("!!! "+this.getClass());
    }

    

    @Override
    public void initialiseCounts() {
        for(int j=0; j<this.dist.length; j++){
            dist[j].initialiseCounts();
        }
    }

  
  
public int sample(int i){
   int[] sample = new int[this.dist.length];
    for(int j=0; j<this.dist.length; j++){
        sample[j] = this.getInnerState(j, false).sample(i);
    }
    return ((CompoundEmissionStateSpace)this.emStSp).getIndex(sample);
}
  

    
    public void print(PrintWriter pw, String prefix, List<Integer>cols){
      
        for(int i=0; i<dist.length; i++){
            dist[i].print(pw, prefix+" "+dist[i].getName(), cols);
        }
    }
   
    
   
    
    @Override
    public void validate() throws Exception {
        for(int k=0; k<dist.length; k++){
            this.dist[k].validate();
        }
        this.lengthDistrib.validate();
        for(int i=0; i<this.noSnps(); i++){
            double[] sum=new double[emStSp.size()];;
           /// List haploList = emStSp.haploList;
            for(int j=0; j<this.emStSp.size(); j++){
                sum[j]=this.score(j,i);
            }
            if(Math.abs(Constants.sum(sum)-1.0)>0.1) {
            	 for(int j=0; j<this.emStSp.size(); j++){
                     sum[j]=this.score(j,i);
                 }
                throw new RuntimeException("!! "+sum);
            }
        }
    }

    @Override
    public EmissionState[] getMemberStates(boolean real) {
       return dist;
    }
  
   
    
  /*  public void initialise(StateIndices dat) {
        if(this.dist[0]!=this.dist[1]) throw new RuntimeException("!!");
        for(int i=0; i<this.noSnps(); i++){
            double[] comp =  dat.getDistribution(i);
           EmissionStateSpace emStSp =getEmissionStateSpace();
            double[] init = new double[Emiss.stateSpace.length];
            double pseudo = 0.0;
            Arrays.fill(init, pseudo);
           
           outer: for(int k=0; k<comp.length; k++){
               ComparableArray key = (ComparableArray) emStSp.get(k);
               for(int j=0; j<key.size(); j++){
                   init[this.getEmissionStateSpace(j).get(key.get(j))]+=comp[k];
                 //       getParentModel(j).getEmissionStateSpaceIndex(key.get(j))]+= comp[k];
               }
                  
               }
            double sum=0;
            for(int k1=0; k1<init.length; k1++){
                sum+=init[k1] ;
            }
           
            for(int k1=0; k1<init.length; k1++){
                init[k1] = init[k1]/sum;
            }
           System.arraycopy(init, 0,((HaplotypeEmissionState)this.dist[0]).emissions[i].probs, 0, init.length);
           }
         
        
    }*/

    @Override
    public EmissionStateSpace getEmissionStateSpace() {
        return emStSp;
    }

    @Override
    public int getParamIndex() {
      int max=0;
      for(int i=0; i<this.dist.length; i++){
          int in = dist[i].getParamIndex();
          if(in>max) max = in;
      }
      return max;
    }

	@Override
	public void modifyDirectCounts(double d) {
		if(true) throw new RuntimeException("!!");
	 for(int i=0; i<dist.length; i++){
		PseudoDistribution[] dists = ((HaplotypeEmissionState) dist[i]).emissions;
		for(int i1=0; i1<dists.length; i1++){
			dists[i1].multiplyCounts(d);
		}
	 }
		
	}
   
    


   
    
    
    
   // static final double fraction_null = 1.0;
}
