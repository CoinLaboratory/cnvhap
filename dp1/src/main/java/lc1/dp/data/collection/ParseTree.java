package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.Sequence;


public class ParseTree {
  PhylogenyNode root;
  SortedMap<String, PhylogenyNode> m = new TreeMap<String,PhylogenyNode>();
  
  private PhylogenyNode makeNode(String str) throws Exception {
		PhylogenyNode n =  new PhylogenyNode(str);
		Sequence s = new Sequence();
		DomainArchitecture d = new DomainArchitecture();
		s.setDomainArchitecture(d);
		n.getNodeData().setSequence(s);
		return n;
	}
  
  ParseTree(File yt, List<String> samples){
		try{
			BufferedReader br1 = new BufferedReader(new FileReader(yt));
			
			String st = "";
			
			m.put("root", root = makeNode("root"));
			while((st = br1.readLine())!=null){
				String[] str = st.split("\\s+");
				PhylogenyNode par = m.get(str[1]);
				PhylogenyNode chi = m.get(str[0]);
				if(chi==null) m.put(str[0], chi = makeNode(str[0]));
				setChildParent(chi,par);
			}
			for(int i=0; i<samples.size(); i++){
				addNode(samples.get(i));
			}
	//	PhylogenyNode no = 	this.m.get("D1a1").getParent().getParent().getParent().getParent().getParent();
		//prune(this.root);
		}catch(Exception exc){
			exc.printStackTrace();
		}
  }
  
 private void prune(PhylogenyNode n){
	 if(n.isExternal()) return;
	 int cnt =  n.getNumberOfDescendants();
	 if(cnt ==1){
		
		PhylogenyNode chi = n.getChildNode(0);
		PhylogenyNode par =  n.getParent();
		this.setChildParent(chi, par);
		String nme = n.getNodeName();
		PhylogenyNode removed = this.m.remove(nme);
		n.removeChildNode(0);
		//if(removed==null){
			//throw new RuntimeException("!!");
		//}
		//else {
		//	System.err.println("removed "+nme);
		//}
	 }else{
		for(int i=0; i<cnt; i++){
			prune(n.getChildNode(i));
		}
	 }
	 
  }

 private void setChildParent(PhylogenyNode c, PhylogenyNode p){
	 c.setParent(p);
	 p.addAsChild(c);
	 /*String nme = p.getNodeName();
	 if(nme.equals("Q1a1") ){//&& c.getNodeName().equals("MNOPS")){
		 System.err.println("h");
	 }*/
 }
 
private PhylogenyNode addNode(String string) throws Exception{
	PhylogenyNode n = this.m.get(string);
	if(n==null){
		m.put(string, n = makeNode(string));
		if(string.length()==1){
			System.err.println(" should be root");
		}
		setChildParent(n,addNode(string.substring(0,string.length()-1)));
	}
	return n;
}

/** adds all children */
public void expand(List<String> n, String st) {
	//for(Iterator<String> it = m.keySet().iterator(); it.hasNext();){
	 String key = st;//it.next();
	 if(st.endsWith("*")) key = st.substring(0,st.length()-1);
	// if(key.startsWith(st)){
	PhylogenyNode pn = m.get(key);
	List<String> ln = pn.getAllExternalDescendantsNames();
	for(int i=0; i<ln.size(); i++){
		if(!n.contains(ln.get(i)))n.add(ln.get(i));
	}
	//}
	//}
}

public boolean isAncestral(String nme, String ke) {
	PhylogenyNode pn = this.m.get(nme);
	PhylogenyNode pn1 = this.m.get(ke);
	//if(ke.equals("P") && nme.equals("K")){
	//	System.err.println("h");
	//}
	if(pn==root) return true;
	while(pn1!=this.root){
		if(pn1.equals(pn)){
			return true;
		}
		pn1 = pn1.getParent();
		//System.err.println(pn1.getNodeName());
	}
	return false;
}
}
