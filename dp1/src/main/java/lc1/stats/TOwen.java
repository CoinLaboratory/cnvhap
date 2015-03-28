package lc1.stats;
import cern.colt.function.DoubleDoubleFunction;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
public class TOwen {
 
        DoubleMatrix1D h;
        double a;
        Algebra alg = new Algebra();
        int jmax, cutpoint;
        TOwen(double[] h, double a, int jmax, int cutpoint){
            this.h = new DenseDoubleMatrix1D(h);
            this.a =a ;
            this.jmax = jmax ;
            this.cutpoint = cutpoint;
        }
        
        public static DoubleMatrix2D multOuter(DoubleMatrix1D a, DoubleMatrix1D b, DoubleDoubleFunction f){
            double[][]res = new double[a.size()][b.size()];
            for(int i=0; i<res.length; i++){
                for(int j=0; j<res[i].length; j++){
                    res[i][j] = f.apply(a.get(i), b.get(j));
                }
            }
            return new DenseDoubleMatrix2D(res);
        }
        
  class Fui implements DoubleDoubleFunction{
      public double apply(double h, double i){
      return Math.pow(h, 2*i)/ ((Math.pow(2,i))*cern.jet.stat.Gamma.gamma(i+1)) ;
      }
  }
  double[] makevec(int st, int end){
     double[] res = new double[end-st];
     for(int i=0; i<res.length; i++){
         res[i] = st+i;
     }
     return res;
  }
  class LowProcedure implements cern.colt.function.DoubleProcedure{
      boolean low = true;
      public LowProcedure(boolean low){
          this.low = low;
      }
    public boolean apply(double arg0) {
        return arg0 < cutpoint;
    }
      
  }
  
 /* double[] getRowSums(DoubleMatrix2D){
      
  }
  
     double tInt()
       {
       Double seriesL, seriesH;
       seriesL = null; seriesH = null;
       DoubleMatrix1D  i = new DenseDoubleMatrix1D(makevec(0, jmax));
       DoubleMatrix1D hL = this.h.viewSelection(new LowProcedure(true));
       DoubleMatrix1D hH = this.h.viewSelection(new LowProcedure(false));
       double  L =hL.size();
         if (L>0) {
           DoubleMatrix2D b  = multOuter(hL,i, new Fui());
           cumb <- apply(b,1,cumsum)
           b1   <- exp(-0.5*hL^2)*t(cumb)
           matr <- matrix(1,jmax+1,L)-t(b1)
           jk   <- rep(c(1,-1),jmax)[1:(jmax+1)]/(2*i+1)
           matr <- t(matr*jk) %*%  a^(2*i+1)
           seriesL  <- (atan(a)-as.vector(matr))/(2*pi)
         }
         if (length(hH) >0) 
           seriesH <- atan(a)*exp(-0.5*(hH^2)*a/atan(a))*
                      (1+0.00868*(hH^4)*a^4)/(2*pi)
         series <- c(seriesL,seriesH)
         id <- c((1:length(h))[low],(1:length(h))[!low]) 
         series[id] <- series  # re-sets in original order
         series
      }
      if(!is.vector(a) | length(a)>1) stop("a must be a vector of length 1")
      if(!is.vector(h)) stop("h must be a vector")
      aa <- abs(a)    
      ah <- abs(h)
      if(is.na(aa)) stop("parameter 'a' is NA") 
      if(aa==Inf) return(0.5*pnorm(-ah))
      if(aa==0)   return(rep(0,length(h)))
      na  <- is.na(h)
      inf <- (ah==Inf)
      ah  <- replace(ah,(na|inf),0)
      if(aa<=1)
        owen <- T.int(ah,aa,jmax,cut.point)
      else
        owen<-0.5*pnorm(ah)+pnorm(aa*ah)*(0.5-pnorm(ah))- 
                   T.int(aa*ah,(1/aa),jmax,cut.point)
      owen <- replace(owen,na,NA)
      owen <- replace(owen,inf,0)
      return(owen*sign(a))
    }*/
}
