package lc1.dp.data.representation;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import lc1.dp.emissionspace.EmissionStateSpace;
import lc1.dp.states.PhasedDataState;
import lc1.util.Constants;



public  abstract class SimpleScorableObject implements ScorableObject, SSOData {
    
    public static PhasedDataState make(String string, List<String> asList,
            EmissionStateSpace emStSp, short dataIndex){
        return new PhasedDataState(string, asList, emStSp, dataIndex);
        //return new PhasedIntegerGenotypeData(string, asList, class1);
    }
    
    
    public static PhasedDataState make(String string, int noSnps, EmissionStateSpace emStSp, short data_index) {
        // TODO Auto-generated method stub
        return new PhasedDataState(string, noSnps, emStSp, data_index);
    }
    List<Comparable> l1 ;
    String id;
    final Class clazz;
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#getName()
     */
    public String getName(){
        return id;
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#setName(java.lang.String)
     */
    public void setName(String name) {
        this.id = name;
        
    }

    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#clone()
     */
    public abstract ScorableObject clone();
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#remove(int)
     */
    public final void remove(int i){
        this.l1.remove(i);
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#getStringRep(int)
     */
    public  String getStringRep(int i){
        return this.getElement(i).toString();
    }
    
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#restrictTo(java.lang.Integer, java.lang.Integer)
     */
    public void restrictTo(Integer integer, Integer integer2) {
        this.l1 = l1.subList(integer, integer2+1);
    }

    
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#restrictSites(int)
     */
    public void restrictSites(int i) {
       if(i<l1.size()){
           this.l1 = l1.subList(0,i);
       }
        
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#restrictSites(int, int)
     */
    public void restrictSites(int min, int max) {
            this.l1 = l1.subList(min, max);
         
     }
    
public SimpleScorableObject(String id, int noSites, Class clazz){
    this.l1 = new ArrayList(noSites);
    this.id = id;
    this.clazz = clazz;
}
public Class clazz(){
    return this.clazz;
}

public static void printIdLine(String idLine, PrintWriter pw, int len){
    char[] ch = new char[Math.max(idLine.length(),len)];
    Arrays.fill(ch, ' ');
    System.arraycopy(idLine.toCharArray(), 0, ch, 0, idLine.length());
  int jmp =10;
  int jmp2 = 100;
    for(int i=0; i<ch.length; i+=jmp){
        if(Math.IEEEremainder(i, jmp2)==0 && idLine.length()+i < ch.length){
            System.arraycopy(idLine.toCharArray(), 0, ch, i, idLine.length());
        }
        if(i>=idLine.length()){
            ch[i] = '|';
        }
    }
    pw.println(new String(ch));
}
public SimpleScorableObject(SSOData data){
    this(data.getName(), data.length(), data.clazz());
    for(int i=0; i<data.length(); i++){
        this.l1.add(copyElement(data.getElement(i)));
    }
}
/* (non-Javadoc)
 * @see lc1.dp.data.representation.SSOData#copyElement(java.lang.Comparable)
 */
public abstract Comparable copyElement(Comparable element);

/* (non-Javadoc)
 * @see lc1.dp.data.representation.SSOData#addMissingData(double)
 */
public  int addMissingData(double perc){
    int cnt =0;
    for(int i=0; i<this.length(); i++){
            if(Constants.rand.nextDouble()<perc){
                this.set(i, null);
                cnt++;
            }
            else if(Constants.rand.nextDouble()<perc){
               throw new RuntimeException();
            }
    }
    return cnt;
}



/* (non-Javadoc)
 * @see lc1.dp.data.representation.SSOData#allNull()
 */
public final  boolean allNull(){
    for(int i=0; i<l1.size(); i++){
        if(l1.get(i)!=null) return false;
    }
    return true;
}




    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#length()
     */
    public final int length() {
       return l1.size();
    }

    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#getElement(int)
     */
    public  Comparable getElement(int i) {
        return l1.get(i);
    }
    
    
/* (non-Javadoc)
 * @see lc1.dp.data.representation.SSOData#printElement(java.io.PrintWriter, java.lang.Object, boolean)
 */
public void printElement(PrintWriter pw, Object el, boolean expand){
    if(el==null) pw.print(" ");
    else{ 
        pw.print( el instanceof Emiss ? 
                (expand ? ((Emiss)el).toStringShort() : ((Emiss)el).toStringPrint()) : 
        el instanceof Integer ? ((Integer)el).toString((Integer)el, Constants.radix()) : el.toString());
    }
}

    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#print(boolean, java.io.PrintWriter)
     */
    public void print(boolean idline, PrintWriter pw){
        if(idline) pw.println("# id "+this.getName());
        for(int i=0; i<this.length(); i++){
           printElement(pw, this.l1.get(i), false);
        }
        pw.println();
    }
   
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#addPoint(java.lang.Comparable)
     */
    public  void addPoint(int i, Comparable i1){
        l1.add(i1);
      //  if(l1.size()>87) throw new RuntimeException("!!");
    }
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#set(int, java.lang.Comparable)
     */
    public  void set(int i, Comparable obj) {
       // if(obj!=null && !clazz.isInstance(obj)) throw new ClassCastException(""+obj.getClass()+" "+clazz);
       l1.set(i, obj);
        
    }
  /*  public void print(PrintWriter pw) {
        for(int i=0; i<l1.size(); i++){
            Object num =  l1.get(i);
            pw.print(num==null ? '?' : num);
        }
        
    }*/
    
    /* (non-Javadoc)
     * @see lc1.dp.data.representation.SSOData#toString()
     */
    public  String toString(){
        if(this.getElement(0) instanceof Emiss){
            StringBuffer sb = new StringBuffer();
            for(int i=0; i<this.length(); i++){
                sb.append(this.getStringRep(i));
            }
            return sb.toString();
        }
        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        print(pw,true, false, null, null, null, null);
        return sw.getBuffer().toString();
    }


    public static PIGData make(PhasedDataState[] unit, boolean b, String join, EmissionStateSpace emstsp) {
       // return new PhasedDataState(unit, b, join);
       return new PhasedDataState(unit, b, join, emstsp);
    }


    public static PhasedDataState make(PIGData next) {
       return new PhasedDataState((PhasedDataState)next);
    }


    public static PhasedDataState make(PIGData[] statesD, List<int[]> list) {
       
        return new PhasedDataState(statesD, list);
        // return new PhasedIntegerGenotypeData(statesD, list);
    }
  
    public static Comparable switchAlleles(Comparable comp){
        if(comp instanceof Emiss){
           return  Emiss.switchElement((Emiss)comp);
         }
         else{
             ComparableArray comp1 = ((ComparableArray) comp).copy();
             for(int j=0; j<comp1.size(); j++){
                
            	 comp1.set(j, switchAlleles(comp1.get(j)));
             }
             return comp1;
         }
    }

}
