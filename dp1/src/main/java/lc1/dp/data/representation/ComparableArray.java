package lc1.dp.data.representation;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import lc1.util.Constants;


/** equals is inconsistent with ordering */
public class ComparableArray  implements AbstractEmiss, Serializable {
    final private List<Comparable> elements;
  
   // private Comparable[] compString;
    
//public Comparable[] compString(){
  //  return compString;
//}
    public Comparable getReal(int j){
        return this.get(j);
    }
 /*   public boolean isEqualGenotype( ComparableArray comp2){
        if(size()!=comp2.size()) return false;
        else return this.compareTo(comp2)==0;

    }*/
    
    public int size(){
        return elements.size();
    }
    public Comparable get(int i){
        return elements.get(i);
    }
    public ComparableArray(boolean b){
        throw new RuntimeException("!!");
    }
    
  /*  private void writeObject(java.io.ObjectOutputStream out)
    throws IOException{
        out.writeChar('(');
        for(int i=0; i<this.size(); i++){
            out.writeObject(this.get(i));
        }
        out.writeChar(')');
    }
private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException{
    
}*/
    
    public ComparableArray copy(){
       return  new ComparableArray(this.elements);
   /*     for(int i=0; i<this.size(); i++){
            Comparable compa = this.get(i);
            if(compa instanceof ComparableArray){
                res.add(((ComparableArray)compa).copy());
            }
            else res.add(compa);
        }
        return res;*/
    }
    
   final  boolean order_known;
private transient Comparator<String> reverse = new Comparator<String>(){

	@Override
	public int compare(String arg0, String arg1) {
		int l1 = arg0.length();
		int l2 = arg1.length();
		if(l1!=l2) return l1<l2  ? -1 : 1;
		else return arg0.compareTo(arg1);
	}
	
};
    
  /*  public boolean equals(Object o){
       // System.err.println(o.getClass());
        Class[] clasz = o.getClass().getInterfaces();
        List co = (List) o;
        if(co.size()!=this.size()) return false;
        for(int i=0; i<co.size(); i++){
            if(!co.get(i).equals(this.get(i))) return false;
        }
        return true;
    }*/
  /*  @Override 
    public void add(int index, Comparable comp){
        super.add(index, comp);
        recalcString();
    }
    
    @Override 
    public void remove(Object comp){
        super.remove( comp);
        recalcString();
    }
    @Override 
    public void add(Object comp){
        super.remove(o);
        super.add( comp);
        recalcString();
    }*/
   
    
  /*  private Comparator getComparator(){
        if(order_known){
            if(this.get(0) instanceof ComparableArray){
                return ORDER_UNORDER;
            }
            else return this.ORDER;
        }
        else return this.IGNORE_ORDER;
    }*/
    
    public void set(int i,Comparable emiss) {
       elements.set(i, emiss);
        
    }
  /*  public int compareTo(Object o) {
        ComparableArray compa = (ComparableArray)o;
   
        if(this.size()!=compa.size()) throw new RuntimeException("sizes must be the same");
        //if(this.compString.length!=compa.compString.length) {
        //    return this.compString.length<compa.compString.length  ? -1 : 1;
      //  }
        for(int i=0; i<compString.length; i++){
            int res = this.compString[i].compareTo(compa.compString[i]);
            if(res!=0) return res;
        }
       return 0;
        /* if(order_known){
           for(int i=0; i<this.size(); i++){
                int res = this.get(i).compareTo(compa.get(i));
                if(res!=0) return res;
            }
            return 0;
        }
        else{
            ComparableArray ex1 = ComparableArrayHelper.perm.get(this).get(0);
            ComparableArray ex2 = ComparableArrayHelper.perm.get(compa).get(0);
            for(int i=0; i<ex1.size(); i++){
                int res = ex1.get(i).compareTo(ex2.get(i));
                if(res!=0) return res;
            }
            return 0;
            
        }
  //      return getComparator().compare(this, o);
     }*/
    
  
    
    /*private static class IgnoreOrderComparator implements Comparator<ComparableArray>{
        public int compare(ComparableArray c1, ComparableArray c2) {
           ComparableArray ex1 = ComparableArrayHelper.perm.get(c1).get(0);
           ComparableArray ex2 = ComparableArrayHelper.perm.get(c2).get(0);
           return ORDER.compare(ex1, ex2);
        }
    };*/
   
   // public static Comparator IGNORE_ORDER = new IgnoreOrderComparator();
   
    
   /* private  static Comparator getComparator(final Comparator inner){
        return new Comparator<ComparableArray>(){
            public int compare(ComparableArray c1, ComparableArray c2) {
                if(c1.size()!=c2.size()){
                    throw new RuntimeException("sizes not equal!!!");
                }
                Iterator<Comparable> o1 = c1.iterator();
                Iterator<Comparable>  o2 = c2.iterator();
              while(o1.hasNext()){
                  Comparable n1 = o1.next();
                  Comparable n2 = o2.next();
                //  if(n1 instanceof ComparableArray){
                  int res =  inner.compare(n1,n2);
                  if(res!=0) return res;
               }
               return 0;
            }
        };
    }
    
    public static Comparator ORDER_UNORDER = getComparator(IGNORE_ORDER);
    */
    
 public  static ComparableArray make(Comparable comp){
       return new ComparableArray(new Comparable[] {comp});
//       compA.add(comp);
 //      return compA;
   }
 
 
 public  boolean containsNull(){
     for(int i=0; i<size(); i++){
         if(get(i) instanceof ComparableArray && ((ComparableArray)get(i)).containsNull()) return true;
         else if(get(i).equals(Emiss.N())) return true;
     }
     return false;
 }
 
 public  int countNull(){
     int count=0;
     for(int i=0; i<size(); i++){
         if(get(i) instanceof ComparableArray) count+=((ComparableArray)get(i)).countNull() ;
         else if(get(i).equals(Emiss.N())) count++;
     }
     return count;
 }

 
public void addObjectsRecursive(List<Comparable> l){
    
        for(Iterator<Comparable> it = this.iterator(); it.hasNext();){
            Comparable nxt =  it.next();
            if(nxt instanceof ComparableArray) ((ComparableArray) nxt).addObjectsRecursive(l);
            else{
               /* if(nxt instanceof Emiss){
                    ComparableArray compa  = ((Emiss)nxt).expand();
                   for(int i=0; i<compa.size(); i++){
                        l.add(compa.get(i));
                    }
                }
                else{*/
                    l.add(nxt);
               // }
            }
        }
 }
public Iterator<Comparable> iterator() {
  return elements.iterator();
}
public static ComparableArray make(Comparable ma, Comparable ma2) {
     return new ComparableArray(new Comparable[] {ma, ma2});
    /* ComparableArray compA = new ComparableArray(false);
     compA.add(ma); compA.add(ma2);
     return compA;*/
   }
 
 public ComparableArray(Comparable[] objects){
     this.elements = new ArrayList<Comparable>(objects.length);
     for(int i=0; i<objects.length; i++){
        elements.add(objects[i]);
     }
   this.order_known = objects[0] instanceof ComparableArray;
  // this.recalcString();
 }
   /* public ComparableArray(boolean order_known){
        super(2);
        this.order_known = order_known;
    }*/

    public ComparableArray(List name) {
       this((Comparable[])name.toArray(new Comparable[0]));
    }
    
  
public ComparableArray(ComparableArray obj) {
        this(obj.elements);
    }
public void flatten(List<Comparable> l){
    for(Iterator<Comparable> it = this.iterator(); it.hasNext();){
        Comparable obj = it.next();
        if(obj instanceof ComparableArray){
            ((ComparableArray)obj).flatten(l);
        }
        else l.add(obj);
    }
}

    public int noCopies(boolean expandEmiss) {
        int no_cop=0;
        for(Iterator it = this.iterator(); it.hasNext();){
            Object obj = it.next();
            if(obj instanceof ComparableArray){
                no_cop+= ((ComparableArray)obj).noCopies(expandEmiss);
            }
            else if (obj instanceof IntegerEmiss){
                no_cop+=1;
            }
            else {
                no_cop+= expandEmiss ? ((Emiss)obj).noCopies() : 1;
            }
        }
        return no_cop;
    }



   

    
   /* public void shuffle() {
        if(this.size()>0 && this.get(0) instanceof ComparableArray){
            for(Iterator it = this.iterator(); it.hasNext();){
                ((ComparableArray)it.next()).shuffle();
            }
        }
        else{
            Collections.shuffle(elements);
        }
    }*/
    
    /** return all possible re-arrangements of key  - i.e. [0,1] -> {[1,0],[0,1]}
     * note: have to fix this - does not do recursive if order not known
     * 
    public  List<ComparableArray> decompose() {
        if(this.order_known){
                return decomposeRecursive();
        }
        else{
            return ComparableArrayHelper.perm.get(this);
        }
     }*/
    
    public boolean needsPhasing() {
       if(this.get(0) instanceof ComparableArray){
           for(int i=0; i<this.size(); i++){
               if(((ComparableArray)this.get(i)).needsPhasing()) return true;//!=this.get(0)) return true;
           }
           return false;
       }
       else{
           for(int i=1; i<this.size(); i++){
               if(this.get(i)!=this.get(0)) return true;
           }
           return false;
       }
    }
   
    
 /*   public List<ComparableArray> decompose(boolean recursive){
        if(recursive){
            if(!this.order_known) throw new RuntimeException("!!");
            return decomposeRecursive();
        }
        else{
            if(this.order_known){
                ArrayList<ComparableArray> a = new ArrayList<ComparableArray>();
                a.add(this);
                return a;
            }
            else{
                return ComparableArrayHelper.perm.get(this);
            }
        }
    }*/
   
    /*private List<ComparableArray> decomposeRecursive() {
        if(!this.order_known) throw new RuntimeException("!!");
        final List<ComparableArray> res = new ArrayList<ComparableArray>();
        CopyEnumerator cnp = new CopyEnumerator(size()){
            @Override
            public void doInner(Comparable[] list) {
              res.add(new ComparableArray(Arrays.asList(list)));
            }

            @Override
            public boolean exclude(Object obj, Object previous) {
                return false;
            }

            @Override
            public Iterator getPossibilities(int depth) {
               Comparable compa = get(depth);
               if(compa instanceof ComparableArray){
                   return ((ComparableArray)compa).decompose().iterator();
               }
               else return  Arrays.asList(new Comparable[] {compa}).iterator();
            }
            
        };
      
        cnp.run();
        return res;
    }*/


    public boolean orderKnown() {
       // if(order_known==true) throw new RuntimeException("!!");
       return order_known;
    }

    public boolean isNested() {
      return this.get(0) instanceof ComparableArray;
    }

    public int homoCount() {
        if(this.isNested()){
            int res =0;
            for(int i=0; i<this.size(); i++){
                res += ((ComparableArray)this.get(i)).homoCount();
            }
            return res;
        }
        else{
            for(int i=1; i<size(); i++){
                if(this.get(i)!=this.get(0)){
                   return 0;
                }
            }
            return 1;
        }
    }

    public void setAsBoolean(Comparable all) {
        for(int k=0;k<size(); k++){
            Comparable nxt =get(k);
            if(nxt instanceof ComparableArray){
                ((ComparableArray)nxt).setAsBoolean(all);
            }
            else{
                if(!nxt.equals(Emiss.N())){
                    if(nxt.equals(all)){
                        nxt = Emiss.b();
                    }
                    else{
                        nxt = Emiss.a();
                    }
                }
                elements.set(k, nxt);
            }
        }
        
    }

    public void convertHemizygousHomozygous() {
        if(this.isNested()){
            for(int j=0; j<size(); j++){
               ((ComparableArray)get(j)).convertHemizygousHomozygous();
            }
        }
        else{
            int countFalse =0;
            int countTrue=0;
            int countNull = countNull();
            if(countNull>0 && countNull<size()){
                for(int j=0; j<size(); j++){
                    if(get(j).equals(Emiss.b())) countTrue++;
                    else if(get(j)==Emiss.a()) countFalse++;
                    else if(get(j)!=Emiss.N()) throw new RuntimeException("!!");
                }
                Emiss val = countTrue>=countFalse ?Emiss.b() : Emiss.a();
                for(int j=0; j<size(); j++){
                    if(get(j)==Emiss.N()) set(j,val);
                }
            }
        }
        
    }

   

    public void addError(double error) {
        if(this.isNested()){
            for(int j=0; j<size(); j++){
               ((ComparableArray)get(j)).addError(error);
            }
        }
        else{
            for(int j=0; j<size(); j++){
                double ra=Constants.rand.nextDouble(); 
                 if(ra < error &&get(j)!=Emiss.N()){
                     set(j, (get(j)==Emiss.b()) ? Emiss.a() : Emiss.b());
                 }
            }
        }
        
    }
    
 public String toString(){
        return this.elements.toString();
    }

     public String toStringShort(){
         StringBuffer sb = new StringBuffer();
         for(int i=0; i<this.size(); i++){
             Object obj  = this.get(i);
             sb.append(obj instanceof Emiss ? ((Emiss)obj).toStringShort() : obj);
         }
         return sb.toString();
     }

     public String toStringPrint(){
         StringBuffer sb = new StringBuffer();
         for(int i=0; i<this.size(); i++){
             Object obj  = this.get(i);
             sb.append(obj instanceof Emiss ? ((Emiss)obj).toStringPrint() : obj);
         }
         return sb.toString();
         
     }
    public int copyNumber() {
       int no = 0;
       for(int i=0; i<this.size(); i++){
          no+= ((Emiss) this.get(i)).noCopies();
       }
       return no;
    }

    public int noB() {
        int no = 0;
        for(int i=0; i<this.size(); i++){
        	if(this.get(i) instanceof AbstractEmiss){
           no+= ((AbstractEmiss) this.get(i)).noB();
        	}
        }
        return no;
        
    }
    
    public int noB(int i1) {
        int no = 0;
        for(int i=0; i<this.size(); i++){
        	if(this.get(i) instanceof AbstractEmiss){
           no+= ((AbstractEmiss) this.get(i)).noB(i1);
        	}
        }
        return no;
        
    }
    public void add(Comparable o) {
        elements.add(o);
   //   recalcString();
        
    }
    
    public String getHaplotypeString(){
       return getString(elements);
    }
    public String getHaploPairString(){
     
          List<Comparable> orig1 = new ArrayList<Comparable>(elements);
         
          Collections.sort(orig1);
         return getString(orig1);
        
    } 
    
    public static String getString(List<Comparable> orig1){
        StringBuffer sb = new StringBuffer();
        
        for(int i=0; i<orig1.size(); i++){
            Comparable comp = orig1.get(i);
            
            if(comp instanceof ComparableArray){
             //   throw new RuntimeException("!!");
                sb.append(((ComparableArray)comp).getGenotypeString());
            }
            else if(comp instanceof IntegerEmiss){
                sb.append( Integer.toString(((IntegerEmiss)comp).v, Constants.radix()));
            }
            else if(comp instanceof Integer){
                sb.append( Integer.toString((Integer)comp, Constants.radix()));
            }
            else{
                sb.append(((Emiss)comp).toStringPrint());
            }
            if(i<orig1.size()-1) sb.append(",");
           // sb.append(this.members[i].getGenotypeString(orig.get(i)));
        }
        String hapString = sb.toString();
        char[] ch = hapString.replaceAll("_", "").toCharArray();
       // Arrays.sort(ch);
        String res = new String(ch);
        return res;
    }
    
    public String getGenotypeString(){
        StringBuffer sb = new StringBuffer();
      // List<String> str = new ArrayList<String>();
        for(int i=0; i<size(); i++){
            Comparable comp = elements.get(i);
          
            	 if(comp instanceof Integer){
                     sb.append( Integer.toString(((Integer)comp), Constants.radix()));
                 }
           
            else
            sb.append(((AbstractEmiss)comp).toStringShort());
           
           // sb.append(this.members[i].getGenotypeString(orig.get(i)));
        }
        String hapString = sb.toString();
        
        char[] ch = hapString.replaceAll("_", "").replaceAll(",", "").toCharArray();
        Arrays.sort(ch);
//        str.add(new String(ch));
  //      Collections.sort(str, reverse);
    //    String str1 = str.toString();
     //   str1 = str1.substring(1,str1.length()-1);
        return new String(ch);//.replaceAll(",", "");
    }
    
    
   
  // private void recalcString() {
   // 
   // }
    public List elements() {
       return elements;
    }
    public Comparable remove(int i) {
       Comparable res = elements.remove(i);
     //  this.recalcString();
       return res;
    }
    public void add(int i, Comparable newFirstPosition) {
       elements.add(i, newFirstPosition);
   //    this.recalcString();
        
    }
    public boolean contains(Comparable comparable) {
        return elements.contains(comparable);
    }
    
     
   
    public boolean equals(Object obj){
    	if(!(obj instanceof ComparableArray)) return false;
        return this.elements.equals(((ComparableArray)obj).elements);
      /*  if(array.elements.size()==this.elements.size()){
            List<Comparable> els = elements;
            List<Comparable> els1 = array.elements;
          for(int i=0; i<els.size(); i++){
              if(!els.get(i).equals(els1.get(i))){
                  return false;
              }
          }
          return true;
        }
        else return false;*/
    }

public int hashCode(){
    return elements.hashCode();
}

public int compareTo(Object o) {
	if(o instanceof Emiss || o instanceof IntegerEmiss) return 1;
    ComparableArray compA = (ComparableArray) o;
    if(compA.size()!=this.size()) return this.size() < compA.size() ? -1 : 1;
    else{
        for(int i=0; i<size(); i++){
           int res = this.get(i).compareTo(compA.get(i));
           if(res!=0) return res;
        }
    }
    return 0;
}

@Override
public int noCopies() {
	int nocop=0;
	for(int i=0; i<this.elements.size(); i++){
		nocop+=((AbstractEmiss)	this.elements.get(i)).noCopies();
	}
	return nocop;
}

public int numLevels() {
    int max=0;
	for(int i=0; i<this.elements.size(); i++){
		Comparable compa = elements.get(i);
		if(compa instanceof ComparableArray){
			int v = ((ComparableArray)compa).numLevels();
			if(v>max){
				max = v;
			}
		}
		
	}
	return max+1;
}

public boolean het() {
	Comparable el = elements.get(0);
	for(int i=1; i<this.elements.size(); i++){
		if(el.compareTo(this.elements().get(i))!=0){
			return true;
		}
	}
	return false;
}


}
