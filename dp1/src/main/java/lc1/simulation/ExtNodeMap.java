package lc1.simulation;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import org.biojava.bio.symbol.Location;

import pal.alignment.Alignment;
import pal.tree.Node;
import pal.tree.Tree;
import cern.jet.random.Poisson;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;

public class ExtNodeMap extends TreeMap<BigInteger,ExtNode> {
    static Poisson poisson = new Poisson(0, new MersenneTwister());
    static Uniform uniform = new Uniform(new MersenneTwister());
   public ExtNodeMap(MSInput ms){
    this(ms.trees.toArray(new Tree[0]), 
            ms.startL.toArray(new Integer[0]),
            ms.end.toArray(new Integer[0])
           );
 
   }
   
  
   
   public void check(Map subset){
	   ExtNodeMap overlap = new ExtNodeMap();
	  ExtNodeMap all = new ExtNodeMap();
	  ExtNodeMap wrong = new ExtNodeMap();
	  all.putAll(this);
	   for(Iterator<Map.Entry<BigInteger, ExtNode>>it = subset.entrySet().iterator(); it.hasNext(); ){
		   Map.Entry<BigInteger, ExtNode> nxt = it.next();
		  ExtNode l =  this.get(nxt.getKey());
		  if(l!=null){
		  Location loc = l.loc;
		  Location loc1 = nxt.getValue().loc;
				  if(loc.overlaps(loc1)){
					  overlap.put(nxt.getKey(), null);
				  }
				  else{
					  all.remove(nxt.getKey());
					  System.err.println("not found location for "+nxt.getKey()+": " +loc+" vs "+loc1);
				  }
		  }
		  else{
			  wrong.put(nxt.getKey(), null);
		  }
	   }
	   System.err.println("detected "+overlap.size()+" of "+all.size()+" "+((double) overlap.size()/(double)all.size())+" with "+wrong.size()+" wrong");
	   
   }
    
   public ExtNodeMap(Alignment al, Integer[] sites, ExtNodeMap truth){
	   super(new Comparator<BigInteger>(){

           public int compare(BigInteger o1, BigInteger o2) {
              return -1*o1.compareTo(o2);
           }
           
       });
	   this.truth = truth;
	   b= new char[al.getSequenceCount()];
	   List<BigInteger> currentNodes = new ArrayList<BigInteger>();
	   List<Boolean> inf = new ArrayList<Boolean>();
       for(int i=0; i<al.getSiteCount(); i++){
           getExtNode(al,i, sites[i], currentNodes,inf);
       }
   }
   
   
   
   /** start, end are the start end coords for each tree */
   public  ExtNodeMap(Tree[] tree,Integer[] start, Integer[] end){
        super(new Comparator<BigInteger>(){

            public int compare(BigInteger o1, BigInteger o2) {
               return -1*o1.compareTo(o2);
            }
            
        });
        for(int i=0; i<tree.length; i++){
            getExtNode(tree[i].getRoot(), start[i], end[i]);
        }
        /*ProbabilityDistribution[] dist = new ProbabilityDistribution[geomean.length];
        for(int i=0; i<dist.length; i++){
            dist[i] =  new LognormalDistribution(Math.log(geomean[i]), Math.log(geodev[i]));
        }
        for(Iterator<Map.Entry<BigInteger,List<ExtNode>>>it = this.entrySet().iterator(); it.hasNext();){
            Map.Entry<BigInteger,List<ExtNode>> nxt = it.next();
            BigInteger key = nxt.getKey();
            List<Integer> extNodes = this.decompose(key);
            List<ExtNode> ss = nxt.getValue();
            Collections.sort(ss);
         //   Logger.global.info("bef "+ss);
            double bl_prv =0;
            for(int i=0; i<ss.size(); i++){
                ExtNode ss_i = ss.get(i);
                int[] noEvents = new int[mutationRate.length];
                double bl = 
                    ss_i.bl() - bl_prv;
                bl_prv = ss_i.bl();
                ss_i.setbl(bl);
                boolean sim =false;
               for(int j=0; j<noEvents.length; j++){
                        noEvents[j] = poisson.nextInt(mutationRate[j]*bl);
                        sim = sim || noEvents[j] >0;
                }
                if(sim){
                    ss_i.simulate(noEvents,  uniform, dist, length);
                }
                
            }
           // Logger.global.info("aft "+ss);
        }*/
    }
    
    public ExtNodeMap() {
	// TODO Auto-generated constructor stub
}

	public static boolean isAncestralTo(BigInteger anc, BigInteger ext){
        return ext.compareTo(ext.and(anc))==0;
    }
   /* public Location[] getMutations(BigInteger anc){
        Location[] res = new Location[3];
        for(Iterator<ExtNode> extnodes = this.get(anc).iterator(); extnodes.hasNext();){
            Location[] mutations = extnodes.next().mutations;
            for(int k=0; k<mutations.length; k++){
                if(mutations[k]!=null){
                    if(res[k]==null){
                        res[k] = mutations[k];
                    }
                    else{
                 //       Logger.global.info("bef "+mutations[k]);
                        res[k] = mutations[k].union(res[k]);
                    //    Logger.global.info("aft "+mutations[k]);
                    }
                }
            }
        }
        return res;
    }
    public Location[] getMutations(Node external){
        Location[] res = new Location[3];
        BigInteger ext = (BigInteger)external.getIdentifier().getAttribute("big");
        for(Iterator<BigInteger> it = getAllAncestral(ext).iterator(); it.hasNext();){
            BigInteger anc = it.next();
            Location[] mutations = getMutations(anc);
            for(int k=0; k<mutations.length; k++){
                    if(mutations[k]!=null){
                        if(res[k]==null){
                            res[k] = mutations[k];
                        }
                        else{
                     //       Logger.global.info("bef "+mutations[k]);
                            res[k] = mutations[k].union(res[k]);
                        //    Logger.global.info("aft "+mutations[k]);
                        }
                    }
            }
        }
        return res;
    }
   public List<BigInteger> getAllAncestral(BigInteger ext){
       List<BigInteger> l = new ArrayList<BigInteger>();
       for(Iterator<BigInteger> it = this.keySet().iterator(); it.hasNext();){
           BigInteger nxt = it.next();
           if(isAncestralTo(nxt, ext)){
               l.add(nxt);
           }
       }
       return l;
   }
 /*   public List<Integer[]>[] extract(List<ExtNode> l, double bl, int start, int end){
        List<Integer[]>[] res = new List[3];
        for(int i=0; i<res.length; i++){
            res[i] =   new ArrayList<Integer[]>();
        }
        double sum=0;
        for(int i=0; i<l.size(); i++){
            ExtNode en = l.get(i);
            List<Integer[]> [] additional = en.extract(start, end);
            for(int j=0; j<res.length; j++){
                if(additional[j]!=null) res[j].addAll(additional[j]);
            }
            sum+=en.bl();
            if(sum>bl) break;
        }
        for(int i=0; i<res.length; i++){
           Collections.sort(res[i], ExtNode.comp);
        }
        return res;
    }*/
    
    
    static BigInteger two = new BigInteger("2");
  static BigInteger getBigInteger(int i){
      return two.pow(i-1);
  }
     BigInteger getExtNode(Node n, int start, int end){
        if(n.getChildCount()>2) throw new RuntimeException("!!");
        BigInteger res;
        if(n.isLeaf()){
            res = getBigInteger(Integer.parseInt( n.getIdentifier().getName()));
          
        }
        else{
            res =  getExtNode(n.getChild(0), start, end).add(getExtNode(n.getChild(1), start, end));
           
        }
        n.getIdentifier().setName(res.toString());
        if(!n.isRoot()){
            ExtNode cf = this.get(res);
            if(cf==null) this.put(res, cf =new ExtNode(n.getBranchLength(), start, end) );
         //   ExtNode r = 
          
            else cf.extend(n.getBranchLength(),start, end);
               
            
        }
        return res;
    }
     
   char[] b = null;
   
   ExtNodeMap truth;
   BigInteger getExtNode(Alignment al , int col, int pos, List<BigInteger> currentNodes, List<Boolean> inf){
         
    	 int cnt=0;
    	 for(int k=0; k<b.length; k++){
    		 b[k] =al.getData(
    				// k
    				 b.length - k-1
    				 , col);
    		 if(b[k]=='1') cnt++;
    	 }
    	 BigInteger res = new BigInteger(new String(b),2);
    	 if(cnt>1 && cnt  <b.length){
       
        
       
             ExtNode cf = this.get(res);
             if(cf==null) this.put(res, cf = new ExtNode(0, pos, pos));
             else{
                     cf.extend(0,pos, pos);
             }
        SortedSet<BigInteger> inc = new TreeSet<BigInteger>();
       List<Integer> toRemove = new ArrayList<Integer>();
         for(int k=0; k<currentNodes.size(); k++){
        	 BigInteger b1 = currentNodes.get(k);
        	 int compare = b1.compareTo(res);
        	 if(compare!=0 &&
        		 !( compare >0 ? compatible(b1,res) :compatible(res,b1) )){
        		 inc.add(b1);
        		// toRemove.add(k);
        	 }
         }
         //now resolve inconsist
         SortedSet<BigInteger> inc1 =new TreeSet<BigInteger>();
         if(inc.size()>0){
        	
        	 BigInteger small = inc.first();
        	
        	 for( Iterator<BigInteger> it = inc.iterator();it.hasNext();){
        		 BigInteger key = it.next();
        		 if(key.and(small).intValue()==0){
        			 key.and(res);
        			 inc1.add(key);//, inc.get(nxt));
        			 it.remove();
        		 }
        	 }
        	
        	 

         }
         {
        	 Set<BigInteger> toadd = new HashSet<BigInteger>();
     		for(int k=0; k<currentNodes.size(); k++){
     			BigInteger curr = currentNodes.get(k);
     			boolean first = inc.contains(curr);
     			boolean second = inc1.contains(curr);
     			BigInteger or = curr.or(res);
     			BigInteger and = curr.and(res);
     			if((first || second) && false){
	        			boolean tr_or = truth.containsKey(or);// && truth.get(or).loc.contains(pos);
	        			boolean tr_and = truth.containsKey(and);// && truth.get(and).loc.contains(pos);
	        			if(tr_or && tr_and){
	        				tr_or = truth.containsKey(or) && truth.get(or).loc.contains(pos);
	        			}
	        			BigInteger repl = null;// tr_or ? or : and;
	        		    if(tr_or) repl = or;
	        		    else if (tr_and) repl = and;
	        		    if(repl!=null){
	        			currentNodes.set(k, repl);
	        			inf.set(k, true);	
	        			ExtNode en = this.get(repl);
	        				if( en==null){
	        					inferred.put(repl, en = new ExtNode(0,pos,pos));
	        				}
	        		    }
	        		    else{
	        		    	toRemove.add(k);
	        		    }
     			}
     			else if(!first && ! second && inf.get(k)){
     				int bitcount=0;
     				boolean tr_or = truth.containsKey(or);
     				if(!currentNodes.contains(or) && (bitcount=or.bitCount())>1 && bitcount < this.b.length && tr_or){
     					toadd.add(or);
     					
     				}
     			/*	if(!currentNodes.contains(and) && (bitcount = and.bitCount())>1 && bitcount < b.length){
     					toadd.add(and);
     				}*/
     			}
     				
     		}
     	        	 Collections.sort(toRemove);
       	 for(int k=toRemove.size()-1; k>=0; k--){
       		 int j = toRemove.get(k);
       		 currentNodes.remove(j);
       		 inf.remove(j);
       	 }    		
     		if(true){
     		for(Iterator<BigInteger> it1 = toadd.iterator(); it1.hasNext();){
     			 BigInteger and = it1.next();
     			currentNodes.add(and);
     			inf.add(false);
					if(!this.containsKey(and)){
						inferred.put(and, new ExtNode(0,pos,pos));
					}
     		 }
     		}
     	 }
         
           if(!currentNodes.contains(res)){
        	   currentNodes.add(res); 
        	   inf.add(false);
           }
    	 }
         return res;
     }
   
   Map<BigInteger, ExtNode> inferred = new TreeMap<BigInteger, ExtNode>();
   /** assume first is bigger */
     private boolean compatible(BigInteger big, BigInteger small) {
    	BigInteger or = big.or(small);
		if(or.subtract(big).intValue()==0){
			return true;
		}
		else if(big.and(small).intValue()==0){
			return true;
		}
		else 
			return false;
	}



	BigInteger getExtNode(Tree tree, int start, int end){
        return getExtNode(tree.getRoot(), start, end);
    }

    public static void main(String[] args){
        System.err.println(decompose(new BigInteger("19")));
    }
    
static      List<Integer>decompose(BigInteger i){
        int len = i.bitLength();
        List<Integer> res = new ArrayList<Integer>();
        for(int k=0; k<len; k++){
            if(i.testBit(k)){
                res.add(k+1);
            }
        }
        return res;
    }
  /*  public void print(PrintWriter pw){
        String[] name = new String[] {"snp", "del", "cnv"};
        for(Iterator<BigInteger> it = this.keySet().iterator(); it.hasNext();){
            BigInteger nxt = it.next();
            Location[] mutations = getMutations(nxt);
            if(allNull(mutations)) continue;
            List<Integer> nodes = decompose(nxt);
            pw.println("Nodes "+nodes);
            inner: for(int i=0; i<name.length; i++){
                if(mutations[i]==null) continue;
                else{
                    pw.println(name[i]+": "+mutations[i]);
                }
            }
        }
    }
    public Location getLocations(int j){
        Location loc=null;
        for(Iterator<BigInteger> it = this.keySet().iterator(); it.hasNext();){
            BigInteger nxt = it.next();
            Location mutations = getMutations(nxt)[j];
            if(mutations==null) continue;
            else if(loc==null) loc = mutations;
            else{
                loc = loc.union(mutations);
            }
        }
        return loc;
    }*/
     boolean allNull(Location[] mutations) {
       for(int i=0; i<mutations.length; i++){
           if(mutations[i]!=null) return false;
       }
       return true;
    }
}
