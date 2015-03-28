package lc1.dp.data.collection;

import java.util.List;

public interface DataCInterface {

	public String headSNP();
	public String head_snp();
    
    public String name();
	public List<Integer>loc();
	public List<String> snpid();
	public String getCompressedString(String key, int i, boolean b, boolean b1);
	public String chrom();
	public List<Character> alleleA();
	public List<Character> alleleB();
	public List<String> indiv();
}
