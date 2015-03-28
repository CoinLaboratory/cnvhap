/**
 * 
 */
package lc1.util;

import java.util.Arrays;
import java.util.Iterator;

public class PseudoIterator implements Iterator<double[]>{ 
    final boolean count;
    static String[] optional_extensions = new String[] {"limit", "start", "index_start", "update_freq", "frac"};
        int index =0; 
        
         public double[] start;
         public double[] index_start;
         public double[] limit;
         public double[] update_freq;
         public double[] frac;
        
        
        public void set(String extension, double[] st) {
            try{
              //  Class clazz = this.getClass();
          //   Field[] fileds =    this.getClass().getFields();
            this.getClass().getField(extension).set(this, st);
            }catch(Exception exc){
                exc.printStackTrace();
            }
            
        }
        
        
      public PseudoIterator(boolean count) {
    	 this(Constants.noPseudoCountParams,count);
        
        }

      
      public PseudoIterator(int len, boolean count) {
    	  start = new double[len];
          limit = new double[len];
          this.count = count;
          index_start = new double[len];
          update_freq = new double[len];
          frac = new double[len];
          Arrays.fill(limit, 0.0);
          Arrays.fill(start, 1e-5);
          Arrays.fill(index_start, 0);
          Arrays.fill(frac, 1.0);
          Arrays.fill(update_freq, 0);
          current = new double[start.length];
          res = new double[start.length];
            // TODO Auto-generated constructor stub
	}


	public void refresh() {
         System.arraycopy(start, 0,current,0, start.length);
          index =0;
      }
      
        //  int phase = ph >= start.length ? start.length-1 : ph;
        double[] current;// =  new double[start.length];
      double[] res ;//= new double[start.length];
        public boolean hasNext() {
            return true;
        }
        
        public double[] next() {
            
            for(int i=0; i<current.length; i++){
            	if(Double.isNaN(start[i]) || Double.isNaN(limit[i])) res[i] = Double.NaN;
            	else{
               if(index == index_start[i]){
                       current[i] = start[i];
              
               }
               else if(index > index_start[i]) current[i] = Math.abs( current[i]*frac[i]);
               else current[i] =count ? 0 : 1e10;
               if(!count && current[i] < limit[i] ) current[i] = limit[i];
               if(count && current[i] > limit[i]) current[i] = limit[i];
               if(Math.abs(Math.IEEEremainder(index, update_freq[i]))>0.0001){
                  res[i] =count ? 0 : 1e6;
               }
               else res[i] = current[i];
            	}
            }
            index++;
           return res;
        }

        public void remove() {
            // TODO Auto-generated method stub
            
        }


      

       
    }