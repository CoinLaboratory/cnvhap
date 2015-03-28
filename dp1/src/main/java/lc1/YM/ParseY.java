package lc1.YM;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;


import org.forester.phylogeny.Phylogeny;
import org.forester.phylogeny.PhylogenyNode;
import org.forester.phylogeny.data.DomainArchitecture;
import org.forester.phylogeny.data.NodeData;
import org.forester.phylogeny.data.ProteinDomain;
import org.forester.phylogeny.data.Sequence;



public class ParseY {
private Phylogeny phyl;

public static void main(String[] args){
	try{
		File f = new File("chrY.csv");
		File yf = new File("ytree.txt");
		ParseY y = new ParseY(f,yf);
		y.makeTree();
		System.err.println("h");
		y.makeAlignment(false);
		File out = new File("chrY.haplogroups");
		y.print(out);
		File dir = new File("builds");
		y.printBranchBuilds(dir);
	}catch(Exception exc){
		exc.printStackTrace();
	}
	
}

short[][] alignment;

String[] extn;

private List<String> getExtNodeNames(){
	Iterator<PhylogenyNode> it = this.m.values().iterator();
	List<String> l = new ArrayList<String>();
	while(it.hasNext()){
		PhylogenyNode n = it.next();
		if(n.isExternal()){
			String st = n.getNodeName();
			l.add(n.getNodeName());
			if(st.equals("BT")){
				System.err.println("h");
			}
			if(n.getDescendants().size()>0){
				throw new RuntimeException("!!");
			}
		}
	}
	return l;
}

private void makeAlignment(boolean all) {
		 extn = 
			getExtNodeNames().toArray(new String[0]);
		if(all){
			List<String> nme = new ArrayList<String>(this.m.keySet());
			nme.remove("root");
			extn=  nme.toArray(new String[0]);
		}
		PhylogenyNode[] pn = new PhylogenyNode[extn.length];
		for(int k=0; k<pn.length; k++){
			pn[k] = this.m.get(extn[k]);
			if(k==35){
				PhylogenyNode pnk = pn[k];
				System.err.println("h");
			}
		}
		Integer[] pos = (new ArrayList<Integer>(snps.keySet())).toArray(new Integer[0]);
		short[][] b = new short[pos.length][extn.length];
		for(int i=0; i<pos.length; i++){
			String group = snps.get(pos[i]);
			PhylogenyNode groupn = this.m.get(group);
			
			for(int j=0; j<extn.length; j++){
				boolean anc =isAncestral(groupn,pn[j]); 
				boolean dec = isAncestral(pn[j],groupn) & ! pn[j].equals(groupn); 
				if(dec){
					PhylogenyNode pnj = pn[j];
					System.err.println("is dec "+groupn+" of "+pn[j]);
				}
				//System.err.println(pn[j]+"\t"+groupn+"\t"+anc);
				b[i][j] = dec ? -1 : (anc ? (short)1 : (short)0);
			}
		}
		alignment = b;
//		List<PhylogenyNode> s = new ArrayList<PhylogenyNode>(this.phyl.getExternalNodes());
	// TODO Auto-generated method stub
}

private boolean isAncestral(PhylogenyNode groupn, PhylogenyNode pn) {
	return groupn.getAllExternalDescendants().contains(pn);
}
 String header="";
public void print(File out) throws Exception{
	PrintWriter pw = new PrintWriter(new FileWriter(out));
	//String[] extn = phyl.getAllExternalNodeNames();
	Integer[] pos = (new ArrayList<Integer>(snps.keySet())).toArray(new Integer[0]);
	pw.print(header.replace('\t', ' '));
	for(int k=0; k<extn.length; k++){
		pw.print(" "+extn[k]+"-geno");
	}
	pw.println();
	for(int i=0; i<pos.length; i++){
		//String[] str = rsids.get(pos[i]).split("\\s");
		//if(!compl(str[1]).equals(str[2])){
		pw.print(this.rsids.get(pos[i]).replace('\t', ' ').replace("chr", ""));
		for(int k=0; k<extn.length; k++){
			pw.print(" "+(alignment[i][k]<0 ? '-' : (this.alignment[i][k]==0 ? 'A' : 'B')));
		}
		pw.println();
		//}
	}
	pw.close();
}
private String compl(String s){
	if(s.equals("A")) return "T";
	if(s.equals("C")) return "G";
	if(s.equals("T")) return "A";
	if(s.equals("G")) return "C";
	else return "";
}

static class OrderedPrint{
	PrintWriter pw;
	OrderedPrint(File f) throws Exception{
		pw = new PrintWriter(new FileWriter(f));
	}
	SortedMap<Integer, String> m1 = new TreeMap<Integer, String>();
	public void println(String st){
		m1.put(Integer.parseInt(st.split("\t")[1]),st);
	}
	public void close(){
	Iterator<String> m2 = m1.values().iterator();
	while(m2.hasNext()){
		String nxt = m2.next();
		pw.println(nxt);
	}
	pw.close();
	}
}

public void printBranchBuilds(File dir) throws Exception{
	dir.mkdir();
	this.printBranchBuild(dir,this.phyl.getRoot(), new ArrayList<OrderedPrint>());
}
public void printBranchBuild(File dir, PhylogenyNode n, List<OrderedPrint> pws1) throws Exception{
	NodeData nd = n.getNodeData();
	boolean recursive = true;
	OrderedPrint pw = new OrderedPrint(new File(dir,n.getNodeName()));
	List<OrderedPrint> pws =  new ArrayList<OrderedPrint>();
	if(recursive)pws.addAll(pws1);
	pws.add(pw);
		if(nd!=null) {
		List<Sequence> s = nd.getSequences();
		if(s!=null && s.size()>0) {
		DomainArchitecture arch = s.get(0).getDomainArchitecture();
		if(arch!=null){
		Iterator< ProteinDomain> m = arch.getDomains().values().iterator();
		
		while(m.hasNext()){
			ProteinDomain dom = m.next();
			int st = dom.getFrom();
			String string = this.rsids.get(st) ;
			int ind0 =string.indexOf('\t',0)+1;
			int ind1=  string.indexOf('\t',ind0);
			Integer pos = Integer.parseInt(string.substring(ind0,ind1));
			Integer pos1 = pos+20;
			//int end = dom.getFrom()+20;
			//String id = "Y_"+st;
			// "Y\t"+st+"\t"+end+"\t"+id
			String left = string.substring(0,ind0);
			String right = string.substring(ind1);
			String middle = pos+"\t"+pos1;
			for(int k=0; k<pws.size(); k++){
				pws.get(k).println(left+middle+right);
			}
			
		}
		
		
		
		}
		}
		}
		for(int i=0; i<n.getNumberOfDescendants(); i++){
			printBranchBuild(dir,n.getChildNode(i), pws);
		}
	pw.close();
}
	SortedMap<String, PhylogenyNode> m = new TreeMap<String, PhylogenyNode>();
	SortedMap<Integer, String> snps = new TreeMap<Integer, String>();
	Map<Integer, String> rsids = new HashMap<Integer, String>();
	PhylogenyNode root;
	
	ParseY(File f, File yt){
		try{
			PrintWriter excluded = new PrintWriter(new FileWriter(new File("excl.txt")));
			BufferedReader br = new BufferedReader(new FileReader(f));
			BufferedReader br1 = new BufferedReader(new FileReader(yt));
			
			String st = "";
			
			m.put("root", root = makeNode("root"));
			while((st = br1.readLine())!=null){
				String[] str = st.split("\\s+");
				PhylogenyNode par = m.get(str[1]);
				PhylogenyNode chi = m.get(str[0]);
				if(chi==null) m.put(str[0], chi = makeNode(str[0]));
				chi.setParent(par);
				par.addAsChild(chi);
			}
			//SNP	Hapgroup	Comments	RefSNP ID	Y-position (Build 36.3)
			//2611	O3a1c  comm rs2075181	7606726
			//
			//
			//chrom   start   end     rsid    identifier      altid   haplogroup      alleleA alleleB anc_eq_ref
			this.header = br.readLine();
			List<String> header = Arrays.asList(this.header.split("\t"));
			//SortedSet<Integer> rest = new TreeSet<Integer>();
			//for(int i=0; i<header.size(); i++) rest.add(i);
			int haplo_index = header.indexOf("haplogroup");
			int pos_index = header.indexOf("start");
			
		//	int rs_id =header.indexOf(arg0)
			//int rs_index = header.indexOf("altid");
			int allA = 5; int allB = 6;
			while((st = br.readLine())!=null){
				if(st.indexOf('#')>=0 || st.indexOf("NA")>=0) {
					continue;
				}
				String[] str = st.split("\t");
				for(int k=0; k<str.length; k++){
					str[k] = str[k].replaceAll("\\s+", "");
				}
				if(true){
					try{
					if(str[haplo_index].indexOf('-')<0 
							//allB<str.length && check(str[allA]) && check(str[allB]) &&
							){
						System.err.println(str[allA]+"->"+str[allB]);
				    String hstr = str[haplo_index];
					PhylogenyNode n = m.get(hstr);
					String pos = str[pos_index].split(";")[0];
					int pos_ = Integer.parseInt(pos.split("\\.")[0]);
					snps.put(pos_, str[haplo_index]);
					rsids.put(pos_, st);
					if(n==null) m.put(str[haplo_index], n = makeNode(str[haplo_index]));
					//else System.err.println("already there "+n.getNodeName());
					update(n, str,pos_);
					}
					else{
						excluded.print(st);
					}
					}catch(Exception exc){
						System.err.println(exc.getMessage());
					}
				}
			}
			excluded.close();
		}catch(Exception exc){
			exc.printStackTrace();
		}
	}
	
	
	
	private boolean check(String string) {
		return string.length()==1 && "ACGT".indexOf(string)>=0;
	}

	public void makeTree()throws Exception{
		List<String> l = new ArrayList<String>(m.keySet());
		l.remove("root");
		this.makeTree(l, root);
		phyl =  new Phylogeny();
		phyl.setRoot(root);
		
	}
	public void makeTree(List<String> l, PhylogenyNode root) {
		//List<String> l = new ArrayList<String>(m.keySet());
		List<List<String>> l1 = new ArrayList<List<String>>();
		List<PhylogenyNode> ln = new ArrayList<PhylogenyNode>();
	
		
		
		while(l.size()>0){
			
			Iterator<String> it = l.iterator();
			String st = it.next();
			
			//System.err.println("h"+l.size()+" "+st);
			PhylogenyNode parent = this.m.get(st);
			PhylogenyNode gp = parent.getParent();
			List<String> l_new = new ArrayList<String>();
			if(gp==null){
				parent.setParent(root);
				root.addAsChild(parent);
			}
			it.remove();
			for(; it.hasNext();){
				String st1 = it.next();
				if(st1.startsWith(st)){
					
				
					l_new.add(st1);
					it.remove();
				}
			}
			if(l_new.size()>0){
				l1.add(l_new);
				ln.add(parent);
			}
			//System.err.println("h"+l.size()+" "+l_new);
		}
		for(int i=0; i<l1.size(); i++){
			makeTree(l1.get(i), ln.get(i));
		}
	}

	private void transferToChild(PhylogenyNode parent, PhylogenyNode child) {
		Iterator<ProteinDomain> it = parent.getNodeData().getSequence().getDomainArchitecture().getDomains().values().iterator();
		DomainArchitecture da = child.getNodeData().getSequence().getDomainArchitecture();
		while(it.hasNext()){
			da.addDomain(it.next());
		}
		
	}

	private void update(PhylogenyNode n, String[] str, int st) {
		//n.s
	
		ProteinDomain d = new ProteinDomain(str[3], st,st+1);
		DomainArchitecture da = n.getNodeData().getSequence().getDomainArchitecture();
		da.addDomain(d);
		//System.err.println(n.getNodeName()+" "+da.getNumberOfDomains());
	}

	private PhylogenyNode makeNode(String str) throws Exception {
		PhylogenyNode n =  new PhylogenyNode(str);
		Sequence s = new Sequence();
		DomainArchitecture d = new DomainArchitecture();
		s.setDomainArchitecture(d);
		n.getNodeData().setSequence(s);
		return n;
	}
	
	
}
