package data;

public class RUtils {
	public static double[][] getIMatrix(int d, double[] phi, int cats, int num1, double[][] x){
		double[][] I = new double[d][d];
		for(int j =0; j< d; j++){
			for(int k =j; k<d; k++){
				int num=0;
				for(int i =0; i<cats; i++){
			//		num = num+ (i-1)*(i-1)*x[,j]*x[,k]*Math.exp(phi[i]);
				}
			//	I[j,k] = sum(num1.get(j))
			}
		}
		return I;
	}
	/*getIMatrix<-function(d,phi,cats,num1,x){
		 I = matrix(ncol=d,nrow=d)
		 for( j in 1:d ){
		    for( k in j:d ){
		      num = 0
		      for( i in 1:cats ) num = num + (i-1)*(i-1)*x[,j]*x[,k]*exp(phi[i])
		      I[j,k] = sum(num1[[j]]*num1[[k]])
		      I[j,k] = I[j,k] - sum(num)
		      I[k,j] = I[j,k]
		    }
		 }
		 I
		}
		getUNum1Matrix<-function(d,y,x,phi,num1,cats){
		  U = matrix(ncol=1,nrow=d)
		  num1 = list(length=d)
		  for( j in 1:d ){U[j] = sum(y*x[,j]);    num1[[j]]=0; for( i in 1:cats ) num1[[j]] = num1[[j]] + x[,j] * (i-1) * exp(phi[i]); U[j,1] = U[j,1] - sum(num1[[j]])}
		  list(U,num1)
		}

	*/
	
}
