package lc1.dp.data.representation;


public interface SSOData extends ScorableObject{

    public abstract String getName();

    public abstract void setName(String name);

  //  public abstract ScorableObject clone();

//    public abstract void remove(int i);

    public abstract String getStringRep(int i);

   // public abstract void restrictTo(Integer integer, Integer integer2);

    public abstract void restrictSites(int i);

    public abstract void restrictSites(int min, int max);


  //  public abstract int addMissingData(double perc);

  // public abstract boolean allNull();

    public abstract int length();

    public abstract Comparable getElement(int i);

    public abstract String printElement( Object el, boolean expand);

   // public abstract void print(boolean idline, PrintWriter pw);

    public abstract void addPoint(int i,Comparable i1);

    public abstract void set(int i, Comparable obj);

    /*  public void print(PrintWriter pw) {
          for(int i=0; i<l1.size(); i++){
              Object num =  l1.get(i);
              pw.print(num==null ? '?' : num);
          }
          
      }*/

    public abstract String toString();

    public abstract Class clazz();

}