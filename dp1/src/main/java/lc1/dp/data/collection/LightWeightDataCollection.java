package lc1.dp.data.collection;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;

import lc1.dp.data.representation.ComparableArray;
import lc1.dp.data.representation.Emiss;
import lc1.dp.data.representation.SimpleScorableObject;
import lc1.dp.states.EmissionState;
import lc1.dp.states.HaplotypeEmissionState;
import lc1.util.ApacheCompressor;
import lc1.util.Constants;

import org.apache.commons.compress.archivers.zip.ZipFile;

/**
 * like likelihood datacollection, but only reads a single snp ata time for mem
 * efficiency
 * 
 * 
 */
public class LightWeightDataCollection extends LikelihoodDataCollection {

	private static final boolean single = false;

	 public static  int midpoint = 0;
	 public static int last = 2*midpoint+1;
	
	public LightWeightDataCollection(File f, short index, int no_copies,
			int[][] mid, File bf, int buff,Collection<String> snpidrest) throws Exception {
		super(f, index, no_copies, mid, bf,snpidrest);
		zf = new ZipFile(f);
		this.alleleA.clear();
		this.alleleB.clear();
		this.strand.clear();
		tmp = new String[this.indiv.size()];
		len = tmp.length;

	}
	public LightWeightDataCollection(File f1, short s, int i, int[][] js,
			List<Integer> locs, List<String> snpid, int buff,Collection<String> snpidrest) throws Exception{
		this(f1, s, i, js,  null, buff,snpidrest);
		this.loc = locs;
		this.snpid = snpid;
	}
	public LightWeightDataCollection(String string, int buff) {
		this.name = name;
		
	}
	public void rename(Map<String, String> convert) {
		throw new RuntimeException("!!");
	}

	@Override
	public void readBuildFile(ZipFile zf, String prefix, BufferedReader br, String chrom,
			final int[][] fromTo, List<Integer> loc, List<String> chr,
			List<String> snpid, List<Character> majorAllele,
			List<Character> minorAllele, List<Boolean> forward, int loc_index,
			int chr_index, int snp_index, int strand_index, int bin,int[] maf_index,
			Collection<String> snp_ids_to_rest)
			throws Exception {
		List<String> l = new ArrayList<String>();

		String st = "";
		outer: for (int i = 0; (st = br.readLine()) != null; i++) {
			String[] str = st.split("\\s+");

			if (chrom.endsWith("all") || str[chr_index].equals(chrom)) {
				int no = Integer.parseInt(str[loc_index]);
				for (int k = 0; k < fromTo.length; k++) {

					if (no >= fromTo[k][0] && no <= fromTo[k][1] && (snp_ids_to_rest==null || snp_ids_to_rest.contains(str[snp_index]))) {
						
						this.process(str, i, no, loc_index, maf_index, chr_index, strand_index, snp_index, l, chr, majorAllele, minorAllele, forward,bin)
						;//this.process(str, i);
						/*l.add(str[loc_index]);
						loc.add(no);
						chr.add(str[chr_index].substring(3));
						snpid.add(str[snp_index]);
						if (maf_index != null && maf_index[0] >= 0
								&& str.length > maf_index[0]) {
							majorAllele.add(str[maf_index[0]].charAt(0));
							minorAllele.add(str[maf_index[1]].charAt(0));
						}*/
						if (single)
							break outer;
						else
							continue outer;
						// System.err.println(no);
					}

				}
			}
		}
		br.close();
		// return readZip(zf, l);

	}

	@Override
	public void createDataStructure(List<String> indiv, List<Integer> ploidy, List<Integer> sampid) {
		for (int i1 = 0; i1 < sampid.size(); i1++) {
			int i = sampid.get(i1);
			String key = indiv.get(i);

			HaplotypeEmissionState value = createEmissionState(key, ploidy
					.get(i));
			value.setNoCop(ploidy.get(i));
			dataL.put(key, value);
			data.put(key, SimpleScorableObject.make(key, 2*midpoint+1, value
					.getEmissionStateSpace(), this.index));

		}
	}
	public void dropIndiv(String[] toDel) {
		if(toDel.length>0) throw new RuntimeException("!!");
	}

	@Override
	public List<Emiss> getHaplotypes(int pos) {
		this.updateIndex(pos);
		return super.getHaplotypes(0);

	}

	@Override
	public void append(int pos, StringBuffer[] sb) {
		this.updateIndex(pos);
		for (int i = 0; i < indiv.size(); i++) {
			ComparableArray comp = (ComparableArray) (data.get(indiv.get(i))
					.getElement(0));
			int no_copies = comp.size();
			for (int k = 0; k < no_copies; k++) {
				sb[i * no_copies + k].append(comp.get(k).toString());
			}
		}
	}

	public void getCompa(int pos, Comparable[] genotypes) {
		this.updateIndex(pos);
		for (int i = 0; i < indiv.size(); i++) {
			genotypes[i] = (data.get(indiv.get(i)).getElement(0));

		}
	}

	@Override
	protected Boolean process(String snpd, int i, ZipFile zf,
			List<Integer> ploidy, List<Integer> samps, double[] missing) {
		return null;
	}

	@Override
	public HaplotypeEmissionState createEmissionState(String key, int no_copies) {
		// if(stSp[1].size()==stSp1[1].size())
		return SimpleScorableObject.make(key, this.last, null, this.index);
		// return new Illumina1NoBg(key, stSp[1],stSp1[1], trans[1],
		// getR(key),getB(key), this.length, index) ;

	}

	int currPosIndex = -1;

	public synchronized String getPhenInfo(String string, int pos_index,
			int phenIndex, int type) {
		updateIndex(pos_index);
		if (pos_index != currPosIndex) {
			this.currentPosScIndex = -1;
			currPosIndex = pos_index;
		}
		return super.getPhenInfo(string, 0, phenIndex, type);
	}

	public String getInfo(String string, int pos_index) {
		updateIndex(pos_index);
		if (string.indexOf("maf") >= 0) {
			return super.getInfo(string, 0);
		} else {
			return super.getInfo(string, pos_index);
		}
	}

	int[] linesToRead = null;

	@Override
	public int restricToAlias(Collection<String> alias) {
		// boolean b = alias.contains("22086");
		// super.restricToAlias(alias)
		List<String> keys = new ArrayList<String>(indiv);
		if (linesToRead != null)
			throw new RuntimeException("can only apply once");
		List<Integer> linesToRead1 = new ArrayList<Integer>();
		for (int i = 0; i < keys.size(); i++) {
			linesToRead1.add(i);
		}
		for (Iterator<String> it = keys.iterator(); it.hasNext();) {
			String key = it.next();

			if (!alias.contains(key)) {
				Integer index = keys.indexOf(key);
				linesToRead1.remove(index); // /note - is removing the object
											// not the index
				dataL.remove(key);
				recSites.remove(key);
				viterbi.remove(key);
				data.remove(key);
				indiv.remove(key);
			}
		}
		int sze = dataL.size();
		System.err.println("sze is " + sze);
		this.linesToRead = new int[linesToRead1.size()];
		for (int i = 0; i < linesToRead.length; i++) {
			linesToRead[i] = linesToRead1.get(i);
		}
		len1 = linesToRead.length;
		return sze;
	}

	String[] tmp;
	int len, len1;

	public boolean canProcess(String snp_id){
		if(snpid==null) return false;
		ZipEntry ent = zf.getEntry(snp_id);
		if(ent==null) return false;
		else{
			
			return true;
		}
	}
	
	public boolean process1(String snp_id, int i, ZipFile zf, int ploidy, double[] missing)
			throws Exception {
		
		Arrays.fill(missing,0);
		// List<String> l = null;
		if (!ApacheCompressor.getIndiv(zf, snp_id, null, tmp))
			return false;
		//if(allnull(tmp)) throw new RuntimeException("!!!");
		if (this.linesToRead == null && tmp.length != indiv.size()) {
			throw new RuntimeException("!!" + tmp.length + " " + indiv.size());
		}
		if (linesToRead == null) {
			for (int j = 0; j < len; j++) {
				String stri = tmp[j];
				String[] st = stri.trim().split("\\s+");

				this.process(indiv.get(j), header, st, i, ploidy,missing);

			}
		} else {
			for (int j = 0; j < linesToRead.length; j++) {
				String stri = tmp[linesToRead[j]];
				String[] st = stri.trim().split("\\s+");

				this.process(indiv.get(j), header, st, i, ploidy,missing);

			}
		}
		return true;
	}

	private boolean allnull(String[] tmp2) {
		// TODO Auto-generated method stub
		return false;
	}

	String currId = "";
double[] missing = new double[2];
	public boolean updateIndex(String rsid) {
		if(rsid==null) return false;
		if (!currId.equals(rsid)) {
			currId = rsid;
			try {
				shuffle();
				//System.err.println("reading "+rsid);
				if (!process1(rsid, last-1, zf, Constants.maxPloidy(), missing))
					return false;
				// super.calculateMaf(false,false);
			} catch (Exception exc) {
				System.err.println("could not process " + rsid);
				exc.printStackTrace();
			}
		}
		return true;
	}

	private void shuffle() {
	for(Iterator<EmissionState> it = dataL.values().iterator(); it.hasNext();){
		HaplotypeEmissionState hes = (HaplotypeEmissionState) it.next();
		for(int i=1; i<hes.emissions.length; i++){
			hes.emissions[i-1] = hes.emissions[i];
		}
		hes.emissions[hes.emissions.length-1] = null;
	}
		
	}
	public String getInfo(String tag, String key, int i, boolean style)
			throws Exception {
		return super.getInfo(tag, key, 0, style);
	}

	@Override
	public void updateIndex(int i) {
		if(true) throw new RuntimeException("!!");
		if (i != currentIndex) {
			currentIndex = i;
			updateIndex(snpid.get(i));

		}
	}

	int currentIndex = 0;
	ZipFile zf;

	public Iterator<String> getRS() throws Exception {

		final BufferedReader br = ApacheCompressor.getBufferedReader(zf, "SNPS");
		return new Iterator<String>() {
			String st = br.readLine();

			public boolean hasNext() {
				if (st == null) {
					try {
						br.close();
					} catch (Exception exc) {
						exc.printStackTrace();
					}
					return false;
				}
				return true;
			}

			public String next() {
				String[] str = st.split("\\s+");
				try {
					st = br.readLine();

				} catch (Exception exc) {
					exc.printStackTrace();
				}
				return str[3];
			}

			public void remove() {
				// TODO Auto-generated method stub

			}

		};
	}

	public String getCompressedString(String key, int i, boolean b,boolean b2) {
		updateIndex(i);
		return super.getCompressedString(key, 0, b,b2);
	}
	/*
	 * public EmissionState makeMafState(EmissionStateSpace emStSp1){ if(emStSp1
	 * instanceof CompoundEmissionStateSpace){ CompoundEmissionStateSpace emStSp
	 * = (CompoundEmissionStateSpace) emStSp1; EmissionStateSpace[] ems =
	 * emStSp.getMembers(); EmissionState[] st = new EmissionState[ems.length];
	 * for(int i=0; i<st.length; i++){ st[i] = new
	 * HaplotypeEmissionState("maf_"+i, 1, ems[i].size(), ems[i], null, null); }
	 * 
	 * this.maf = new AlleleCopyPairEmissionState(Arrays.asList(st), emStSp,
	 * false, null);// emStSp.size(), emStSp); } else{ return new
	 * HaplotypeEmissionState("maf_", 1, emStSp1.size(), emStSp1, null, null); }
	 * maf.initialiseCounts(); return maf; }
	 */

}
