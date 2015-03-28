package lc1.dp.core;



public class Term extends AbstractTerm{
   
       
        public int j;
       public int getBestPath(){
           return j;
       }
       protected Term(int j, int i, double sc){
           super(i,sc);
           this.j = j;
       }
       
      public void setj(int j){
          this.j = j;
      }
       
      public Term(int i, double sc, int modelLength){
          super(i, sc);
      }
       
        
      
    }



        
        

