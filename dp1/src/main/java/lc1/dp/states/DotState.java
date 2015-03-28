/*
 * Created on 17-Aug-2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package lc1.dp.states;

import lc1.stats.SimpleDistribution;


/**
 * @author lc
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class DotState extends State {
    static final long serialVersionUID = 1;
    short num;
    public int getClassId(){
        return 0;
    }
	public DotState(String name, SimpleDistribution start, SimpleDistribution end){
		super(name, 0);
        this.distStart = start;
        this.distEnd = end;
        
	}
    public DotState(DotState st1){
        super(st1);
        this.distStart= new SimpleDistribution(st1.distStart);
        this.distEnd= new SimpleDistribution(st1.distEnd);
    }
    
    
    public String toString(){
        return this.name+" "+this.num;
    }
    public DotState(String name){
        super(name, 0);
        this.distStart = SimpleDistribution.noOffset;
        this.distEnd = SimpleDistribution.noOffset;
        
    }
    public Object clone(){
        return new DotState(this);
    }
    
   public void validate() throws Exception{}  
    public DotState(String name, short num){
        this(name,SimpleDistribution.noOffset, SimpleDistribution.noOffset);
        this.num = num;
    }

 final SimpleDistribution distStart;
 final SimpleDistribution distEnd;

 public int hashCode(){
     return name.hashCode()+num;
 }
public SimpleDistribution adv(State s){
    if(s==null) return distStart;
    else return distEnd;
}

/*public double score(Object element, boolean logspace) {
    return logspace ? 0 : 1;
}*/
}
