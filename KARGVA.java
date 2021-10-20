import java.io.*;
import java.io.File;
import java.lang.*;
import java.math.*;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

class AMRVariantGene
{
	public HashMap<String,Integer> kmerFreque;
	public HashMap<String,Float> kmerMapped;
	public HashMap<String,Integer> kmerIsSNP;
	public AMRVariantGene()
	{
		kmerFreque = new HashMap<String,Integer>();
		kmerMapped = new HashMap<String,Float>();
		kmerIsSNP = new HashMap<String,Integer>();
	}
}

public class KARGVA
{
	public static String checkAndAmendAminoAcidSequence(String s)
	{
		StringBuffer k = new StringBuffer();
		for (int i=0; i<s.length(); i++)
		{
			char c=s.charAt(i);
			if (c=='$' || c=='A' || c=='C' || c=='D' || c=='E' || c=='F' || c=='G' || c=='H' || c=='I' || c=='K' || c=='L' || c=='M' || c=='N' || c=='P' || c=='Q' || c=='R' || c=='S' || c=='T' || c=='V' || c=='W' || c=='Y') {k.append(c);}
			else
				if (c=='a' || c=='c' || c=='d' || c=='e' || c=='f' || c=='g' || c=='h' || c=='i' || c=='k' || c=='l' || c=='m' || c=='n' || c=='p' || c=='q' || c=='r' || c=='s' || c=='t' || c=='v' || c=='w' || c=='y') {k.append(Character.toUpperCase(c));}
				else
					if (c=='*') {k.append('$');}
					else
						if (c!=' ') {k.append('?');}
		}
		return k.toString();		
	}
	
	public static String checkAndAmendRead(String s)
	{
		StringBuffer k = new StringBuffer();
		for (int i=0; i<s.length(); i++)
		{
			char c=s.charAt(i);
			if (c=='A' || c=='a') {k.append('A');}
			else
				if (c=='C' || c=='c') {k.append('C');}
				else
					if (c=='G' || c=='g') {k.append('G');}
					else
						if (c=='T' || c=='t' || c=='U' || c=='u') {k.append('T');}
							else 
								if (c!=' ') {k.append('N');}
		}
		return k.toString();		
	}
	
	public static String reverseComplement(String s)
	{
		char[] reverse = new char[s.length()];
		for (int i=0; i<s.length(); i++) 
		{
			char c = s.charAt(i);
			if (c=='A') {reverse[(reverse.length-1)-i]='T';}
			else
				if (c=='C') {reverse[(reverse.length-1)-i]='G';}
				else
					if (c=='G') {reverse[(reverse.length-1)-i]='C';}
					else
						if (c=='T') {reverse[(reverse.length-1)-i]='A';}
							else 
								if (c=='N') {reverse[(reverse.length-1)-i]='N';}
		}
		return String.valueOf(reverse);
	}
	
	public static String nucleotidesToAminoAcids(String s, HashMap<String,String> h)
	{
		StringBuffer k = new StringBuffer();
		for (int i=0; i<s.length()-2; i=i+3)
		{
			String tri = s.substring(i,i+3); 
			String tra = h.get(tri); 
			if (tra==null) {k.append("?");} else {k.append(tra);}
		}
		return k.toString();
	}
	
	public static String randomString(int n)
	{
		StringBuffer k = new StringBuffer();
		for (int i=0; i<n; i++)
		{
			double d = Math.random();
			if (d<0.000001d) k.append('N');
				else 
				{
					d = Math.random();
					if (d<0.25d) k.append('A');
						else if (d<0.5d) k.append('C');
							else if (d<0.75d) k.append('G');
								else k.append('T');
				}
		}
		return k.toString();
	}
	
	public static Comparator<HashMap.Entry<String,Float>> sortHashMapByValueFloat = new Comparator<HashMap.Entry<String,Float>>()
	{
		@Override
		public int compare(Map.Entry<String,Float> e1, Map.Entry<String,Float> e2)
		{
			Float f1 = e1.getValue();
			Float f2 = e2.getValue();
			return f2.compareTo(f1);
		}
	};
	
	public static void main(String[] args) throws Exception
	{
		long time0 = System.currentTimeMillis();
		long startTime = System.currentTimeMillis();
		long endTime = System.currentTimeMillis();
		long elapsedTime = endTime - startTime;
		final int DEFAULT_BUFFER_SIZE=16384;
		float allram = (float)(Runtime.getRuntime().maxMemory());
		float usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());
		
		String [] nucleotide_triplets ={"AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"};
		String [] amino_acids = {"K","N","K","N","T","T","T","T","R","S","R","S","I","I","M","I","Q","H","Q","H","P","P","P","P","R","R","R","R","L","L","L","L","E","D","E","D","A","A","A","A","G","G","G","G","V","V","V","V","$","Y","$","Y","S","S","S","S","$","C","W","C","L","F","L","F"};
		HashMap<String,String> nuc_ami = new HashMap<String,String>();
		for (int t=0; t<nucleotide_triplets.length; t++) {nuc_ami.put(nucleotide_triplets[t],amino_acids[t]);}
		
		int k = 9;
		int numT = 12500;
		String dbfile="kargva_db_v5.fasta";
		String readfile="";
		boolean reportMultipleHits = true;
		
		for (int t=0; t<args.length; t++)
		{
			if (args[t].startsWith("d:")) dbfile=args[t].split(":")[1];
			if (args[t].endsWith(".fastq") || args[t].endsWith(".gz")) readfile=args[t];
			if (args[t].startsWith("f:")) readfile=args[t].split(":")[1];
			if (args[t].startsWith("k:")) k=Integer.parseInt(args[t].split(":")[1]);
			if (args[t].startsWith("i:")) numT=Integer.parseInt(args[t].split(":")[1]);
			if (args[t].equals("m:n") || args[t].equals("m:no")) reportMultipleHits = false;
			if (args[t].equals("m:y") || args[t].equals("m:yes")) reportMultipleHits = true;
		}
		if (k<4) {System.out.println("Minimum value of k must be 4"); k=4;}
		if (readfile.equals("")) {System.out.println("Please specify a read file"); System.exit(0);}
		
		System.out.println("Reading antibiotic resistance gene variant database, creating k-mer mapping (k="+k+")");
		startTime = System.currentTimeMillis();
		
		HashMap<String,ArrayList<String>> kmerGeneMapping = new HashMap<String,ArrayList<String>>();
		HashMap<String,ArrayList<String>> kmerGeneSNPMapping = new HashMap<String,ArrayList<String>>();
		HashMap<String,AMRVariantGene> geneKmerMapping = new HashMap<String,AMRVariantGene>();
		BufferedReader r = new BufferedReader(new FileReader(dbfile));
		String header = r.readLine();
		long i=0;
		while(true)
		{
			if (!header.startsWith(">")) {System.out.println("Wrong fasta format"); System.exit(0);}
			if (header==null) break;
			String sequence = r.readLine();
			if (sequence==null) break;
			String nextl = r.readLine();
			if (nextl==null) break;
			while(nextl!=null && !nextl.startsWith(">")) {sequence=sequence+nextl; nextl=r.readLine();}
			if (sequence.length()>=k*3)
			{
				AMRVariantGene amrhg = new AMRVariantGene();
				sequence = checkAndAmendAminoAcidSequence(sequence);
				String [] mandPos = header.split("\\|");
				mandPos = mandPos[1].split(";");
				ArrayList<Integer> mandatoryPositions = new ArrayList<Integer>();
				for (int l=0; l<mandPos.length; l++)
				{
					String pp=mandPos[l].replaceAll("[^0-9]","");
					int mpos=Integer.parseInt(pp)-1;
					mandatoryPositions.add(mpos);
					//System.out.println(header.substring(0,header.indexOf("|")));System.out.println(sequence.substring(Math.max(0,mpos-25),Math.min(sequence.length(),mpos+25)));
				}
				for (int g=0; g<sequence.length()-k+1; g++)
				{
					String fk = sequence.substring(g,g+k);
					ArrayList<String> al = kmerGeneMapping.get(fk);
					if (al==null)
					{
						al = new ArrayList<String>();
						al.add(header);
						kmerGeneMapping.put(fk,al);
					}
					else
					{
						if (!al.contains(header)) al.add(header);
						kmerGeneMapping.put(fk,al);
					}
					if (amrhg.kmerFreque.get(fk)==null) {amrhg.kmerFreque.put(fk,1);} else {amrhg.kmerFreque.put(fk,amrhg.kmerFreque.get(fk)+1);}
					for (int l=0; l<mandatoryPositions.size(); l++)
					{
						int mpos = mandatoryPositions.get(l);
						if (g<=mpos && mpos<=(g+k)) //mpos==(g+k)/2
						{
							ArrayList<String> alm = kmerGeneSNPMapping.get(fk);
							if (alm==null)
							{
								alm = new ArrayList<String>();
								alm.add(header);
								kmerGeneSNPMapping.put(fk,alm);
							}
							else
							{
								if (!alm.contains(header)) alm.add(header);
								kmerGeneSNPMapping.put(fk,alm);
							}
							amrhg.kmerIsSNP.put(fk,1);
						}
					}
				}
				geneKmerMapping.put(header,amrhg);
			}
			header=nextl;
			if (nextl==null) break;
			i++;
			if (i%1000==0)
			{
				System.gc();
				allram = (float)(Runtime.getRuntime().maxMemory());
				usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());	
				System.out.println("\t"+i+" genes processed; used RAM = "+100*usedram/allram+"%");
			}
		}
		
		r.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.println(i+" genes read and k-mers mapped in "+elapsedTime/1000+" seconds");
		
		System.out.print("Estimating background/random k-mer match distribution");
		startTime = System.currentTimeMillis();
		if(readfile.endsWith(".gz"))
		{
			r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
		}
		else
		{
			r=new BufferedReader(new FileReader(readfile),DEFAULT_BUFFER_SIZE);
		}
		i=0;
		double avg=0f;
		String line;
		while((line=r.readLine())!=null || i<numT)
		{
			line=r.readLine();
			String fwd = line;
			if (fwd==null) break;
			avg=avg+(double)(fwd.length());
			r.readLine();
			r.readLine();
			i++;
		}
		avg=avg/(double)(i);
		System.out.println(" (average read length is "+Math.round(avg)+" bases)");
		if ( (avg/3d)<k ) {System.out.println("Average read length too short for the chosen k"); System.exit(0);}
		int [] matchDist = new int [numT];
		int [] matchDistNotSNP = new int [numT];
		System.out.print("\t");
		for (int y=0; y<numT; y++)
		{
			String fk = randomString((int)(avg));
			String rk = reverseComplement(fk);
			String f1 = nucleotidesToAminoAcids(fk,nuc_ami); String f2 = nucleotidesToAminoAcids(fk.substring(1),nuc_ami); String f3 = nucleotidesToAminoAcids(fk.substring(2),nuc_ami); 
			String f4 = nucleotidesToAminoAcids(rk,nuc_ami); String f5 = nucleotidesToAminoAcids(rk.substring(1),nuc_ami); String f6 = nucleotidesToAminoAcids(rk.substring(2),nuc_ami); 
			matchDist[y] = 0;
			matchDistNotSNP[y] = 0;
			int hits = 0;
			int hitsNotSNP = 0;
			for (int g=0; g<f1.length()-k+1; g++)
			{
				String q = f1.substring(g,g+k);
				if (kmerGeneSNPMapping.get(q)!=null) {hits++;}
				if (kmerGeneMapping.get(q)!=null) {hitsNotSNP++;}
			}
			if (matchDist[y] < hits) matchDist[y] = hits;
			if (matchDistNotSNP[y] < hitsNotSNP) matchDistNotSNP[y] = hitsNotSNP;
			hits = 0; hitsNotSNP = 0;
			for (int g=0; g<f2.length()-k+1; g++)
			{
				String q = f2.substring(g,g+k);
				if (kmerGeneSNPMapping.get(q)!=null) {hits++;}
				if (kmerGeneMapping.get(q)!=null) {hitsNotSNP++;}
			}
			if (matchDist[y] < hits) matchDist[y] = hits;
			if (matchDistNotSNP[y] < hitsNotSNP) matchDistNotSNP[y] = hitsNotSNP;
			hits = 0; hitsNotSNP = 0;
			for (int g=0; g<f3.length()-k+1; g++)
			{
				String q = f3.substring(g,g+k);
				if (kmerGeneSNPMapping.get(q)!=null) {hits++;}
				if (kmerGeneMapping.get(q)!=null) {hitsNotSNP++;}
			}
			if (matchDist[y] < hits) matchDist[y] = hits;
			if (matchDistNotSNP[y] < hitsNotSNP) matchDistNotSNP[y] = hitsNotSNP;
			hits = 0; hitsNotSNP = 0;
			for (int g=0; g<f4.length()-k+1; g++)
			{
				String q = f4.substring(g,g+k);
				if (kmerGeneSNPMapping.get(q)!=null) {hits++;}
				if (kmerGeneMapping.get(q)!=null) {hitsNotSNP++;}
			}
			if (matchDist[y] < hits) matchDist[y] = hits;
			if (matchDistNotSNP[y] < hitsNotSNP) matchDistNotSNP[y] = hitsNotSNP;
			hits = 0; hitsNotSNP = 0;
			for (int g=0; g<f5.length()-k+1; g++)
			{
				String q = f5.substring(g,g+k);
				if (kmerGeneSNPMapping.get(q)!=null) {hits++;}
				if (kmerGeneMapping.get(q)!=null) {hitsNotSNP++;}
			}
			if (matchDist[y] < hits) matchDist[y] = hits;
			if (matchDistNotSNP[y] < hitsNotSNP) matchDistNotSNP[y] = hitsNotSNP;
			hits = 0; hitsNotSNP = 0;
			for (int g=0; g<f6.length()-k+1; g++)
			{
				String q = f6.substring(g,g+k);
				if (kmerGeneMapping.get(q)!=null) {hits++;}
				if (kmerGeneMapping.get(q)!=null) {hitsNotSNP++;}
			}
			if (matchDist[y] < hits) matchDist[y] = hits;
			if (matchDistNotSNP[y] < hitsNotSNP) matchDistNotSNP[y] = hitsNotSNP;
			if (y%(numT/5)==0) System.out.print(y+"..");
		}
		System.out.println();
		Arrays.sort(matchDist);
		Arrays.sort(matchDistNotSNP);
		int pvalthres=matchDist[99*numT/100];
		int pvalthresNotSNP=matchDistNotSNP[99*numT/100];
		System.out.println("99th percentile of random variant/gene k-mers match distribution is "+pvalthres+"/"+pvalthresNotSNP+" (max is "+matchDist[numT-1]+"/"+matchDistNotSNP[numT-1]+")");
		r.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.println("Empirical distribution for "+numT+" random reads estimated in "+elapsedTime/1000+" seconds");
		
		System.out.println("Reading file and mapping genes");
		startTime = System.currentTimeMillis();
		
		String readOutFile = readfile.substring(0,readfile.indexOf("."))+"_KARGVA_mappedReads.csv";
		FileWriter rfilewriter = new FileWriter(readOutFile);
		BufferedWriter rwriter = new BufferedWriter(rfilewriter);
		rwriter.write("Idx,");
		rwriter.write("GeneAnnotation,");
		rwriter.write("GeneScore/KmerSNPsHits/KmersHitsOnGene/KmersHitsOnAllGenes/KmersTotal");
		if (reportMultipleHits) rwriter.write(",...");
		rwriter.write("\r\n");
		
		if(readfile.endsWith(".gz"))
		{
			r=new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(readfile),DEFAULT_BUFFER_SIZE)),DEFAULT_BUFFER_SIZE);
		}
		else
		{
			r=new BufferedReader(new FileReader(readfile),DEFAULT_BUFFER_SIZE);
		}
		i=0;
		while((line=r.readLine())!=null)
		{
			header = line;
			line=r.readLine();
			String fwd = line;
			i++;
			if (line==null) break;
			r.readLine();
			r.readLine();
			fwd = checkAndAmendRead(fwd);
			if ((fwd.length()/3)>k)
			{
				String rwd = reverseComplement(fwd);
				String f1 = nucleotidesToAminoAcids(fwd,nuc_ami); String f2 = nucleotidesToAminoAcids(fwd.substring(1),nuc_ami); String f3 = nucleotidesToAminoAcids(fwd.substring(2),nuc_ami); 
				String f4 = nucleotidesToAminoAcids(rwd,nuc_ami); String f5 = nucleotidesToAminoAcids(rwd.substring(1),nuc_ami); String f6 = nucleotidesToAminoAcids(rwd.substring(2),nuc_ami); 
				String [] fs = {f1,f2,f3,f4,f5,f6};
				ArrayList<String> kmerHitsMandatoryBest = new ArrayList<String>();
				ArrayList<String> kmerHitsBest = new ArrayList<String>();
				HashMap<String,Float> geneHitsMandatoryWeightedBest = new HashMap<String,Float>();
				HashMap<String,Float> geneHitsWeightedBest = new HashMap<String,Float>();
				HashMap<String,Integer> geneHitsMandatoryUnweightedBest = new HashMap<String,Integer>();
				HashMap<String,Integer> geneHitsUnweightedBest = new HashMap<String,Integer>();
				ArrayList<String> kmerHitsNoSNP = new ArrayList<String>();
				HashMap<String,Float> geneHitsNoSNP = new HashMap<String,Float>();
				int bestMandatory = 0;
				for (int l=0; l<fs.length; l++)
				{
					ArrayList<String> kmerHitsMandatory = new ArrayList<String>();
					ArrayList<String> kmerHits = new ArrayList<String>();
					HashMap<String,Float> geneHitsMandatoryWeighted = new HashMap<String,Float>();
					HashMap<String,Float> geneHitsWeighted = new HashMap<String,Float>();
					HashMap<String,Integer> geneHitsMandatoryUnweighted = new HashMap<String,Integer>();
					HashMap<String,Integer> geneHitsUnweighted = new HashMap<String,Integer>();
					for (int g=0; g<fs[l].length()-k+1; g++)
					{
						String fk = fs[l].substring(g,g+k);
						ArrayList<String> kmerGenesMandatory = kmerGeneSNPMapping.get(fk);
						if (kmerGenesMandatory!=null)
						{
							kmerHitsMandatory.add(fk);
							for (int y=0; y<kmerGenesMandatory.size(); y++)
							{
								String key = kmerGenesMandatory.get(y);
								float frac = 1f/(float)(kmerGenesMandatory.size());
								if (geneHitsMandatoryWeighted.get(key)==null) {geneHitsMandatoryWeighted.put(key,frac);} else {geneHitsMandatoryWeighted.put(key,geneHitsMandatoryWeighted.get(key)+frac);}
								if (geneHitsMandatoryUnweighted.get(key)==null) {geneHitsMandatoryUnweighted.put(key,1);} else {geneHitsMandatoryUnweighted.put(key,geneHitsMandatoryUnweighted.get(key)+1);}
							}
						}
						ArrayList<String> kmerGenes = kmerGeneMapping.get(fk);
						if (kmerGenes!=null)
						{
							kmerHits.add(fk);
							for (int y=0; y<kmerGenes.size(); y++)
							{
								String key = kmerGenes.get(y);
								float frac = 1f/(float)(kmerGenes.size());
								if (geneHitsWeighted.get(key)==null) {geneHitsWeighted.put(key,frac);} else {geneHitsWeighted.put(key,geneHitsWeighted.get(key)+frac);}
								if (geneHitsUnweighted.get(key)==null) {geneHitsUnweighted.put(key,1);} else {geneHitsUnweighted.put(key,geneHitsUnweighted.get(key)+1);}
							}
						}
					}
					if (kmerHitsMandatory.size()>kmerHitsMandatoryBest.size()) {kmerHitsMandatoryBest=kmerHitsMandatory; geneHitsMandatoryWeightedBest=geneHitsMandatoryWeighted; geneHitsMandatoryUnweightedBest=geneHitsMandatoryUnweighted; bestMandatory=l; kmerHitsBest=kmerHits; geneHitsWeightedBest=geneHitsWeighted; geneHitsUnweightedBest=geneHitsUnweighted;}
					if (kmerHitsMandatory.size()==0 && kmerHitsNoSNP.size()<kmerHits.size()) {geneHitsNoSNP=geneHitsWeighted; kmerHitsNoSNP=kmerHits;}
				}
				if (kmerHitsMandatoryBest.size()>pvalthres)
				{
					HashMap<String,Float> scores = new HashMap<String,Float>();
					Set<String> geneKeys = geneHitsMandatoryWeightedBest.keySet();
					for (String gk : geneKeys) 
					{
						float sc1 = geneHitsWeightedBest.get(gk)/(float)(kmerHitsBest.size());
						float sc2 = geneHitsMandatoryWeightedBest.get(gk)/(float)(kmerHitsMandatoryBest.size());
						float sc = sc1 + sc2 - (sc1*sc2);
						scores.put(gk,sc);
					}
					ArrayList<HashMap.Entry<String,Float>> genehitsarr = new ArrayList<HashMap.Entry<String,Float>>();
					for (HashMap.Entry<String,Float> e: scores.entrySet()) {genehitsarr.add(e);}
					Collections.sort(genehitsarr,sortHashMapByValueFloat);
					rwriter.write(header+",");
					float ratio = genehitsarr.get(0).getValue();
					int  stp = 0;
					for (int y=0; y<genehitsarr.size(); y++)
					{
						stp=y+1;
						rwriter.write(genehitsarr.get(y).getKey()+",");
						float fr = genehitsarr.get(y).getValue();
						float fp = (float)Math.round(fr*100)/100;
						rwriter.write(fp+"/"+geneHitsMandatoryUnweightedBest.get(genehitsarr.get(y).getKey())+"/"+geneHitsUnweightedBest.get(genehitsarr.get(y).getKey())+"/"+kmerHitsBest.size()+"/"+(fs[bestMandatory].length()-k+1));
						if (y>19 || fr/ratio<0.05 || !reportMultipleHits) {break;}
						rwriter.write(",");
					}
					rwriter.write("\r\n");
					if (!reportMultipleHits) {stp=1;}
					for (int y=0; y<stp; y++)
					{
						AMRVariantGene genehit = geneKmerMapping.get(genehitsarr.get(y).getKey());
						for (int c=0; c<kmerHitsBest.size(); c++)
						{
							String kh = kmerHitsBest.get(c);
							if (genehit.kmerFreque.get(kh)!=null)
							{
								if (genehit.kmerMapped.get(kh)==null) {genehit.kmerMapped.put(kh,1f);}
								else {genehit.kmerMapped.put(kh,genehit.kmerMapped.get(kh)+1f);}
							}
						}
					}
				}
				else
				{
					rwriter.write(header+",");
					rwriter.write("?,");
					rwriter.write("?/?/?/?/?");
					rwriter.write("\r\n");
				}
				if (kmerHitsMandatoryBest.size()<=pvalthres && kmerHitsNoSNP.size()>pvalthresNotSNP)
				{
					ArrayList<HashMap.Entry<String,Float>> genehitsarr = new ArrayList<HashMap.Entry<String,Float>>();
					for (HashMap.Entry<String,Float> e: geneHitsNoSNP.entrySet()) {genehitsarr.add(e);}
					Collections.sort(genehitsarr,sortHashMapByValueFloat);
					float ratio = genehitsarr.get(0).getValue();
					int  stp = 0;
					for (int y=0; y<genehitsarr.size(); y++)
					{
						stp=y+1;
						float fr = genehitsarr.get(y).getValue();
						if (y>19 || fr/ratio<0.05 || !reportMultipleHits) {break;}
					}
					if (!reportMultipleHits) {stp=1;}
					for (int y=0; y<stp; y++)
					{
						AMRVariantGene genehit = geneKmerMapping.get(genehitsarr.get(y).getKey());
						for (int c=0; c<kmerHitsNoSNP.size(); c++)
						{
							String kh = kmerHitsNoSNP.get(c);
							if (genehit.kmerFreque.get(kh)!=null)
							{
								if (genehit.kmerMapped.get(kh)==null) {genehit.kmerMapped.put(kh,1f);}
								else {genehit.kmerMapped.put(kh,genehit.kmerMapped.get(kh)+1f);}
							}
						}
					}	
				}
			}
			if (i%100000==0)
			{
				System.gc();
				allram = (float)(Runtime.getRuntime().maxMemory());
				usedram = (float)(Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory());	
				endTime = System.currentTimeMillis();
				elapsedTime = endTime - startTime;
				System.out.print(i+" reads processed; used RAM = "+100*usedram/allram+"%; time = "+elapsedTime/1000+" s \r\n");
			}
		}
		r.close();
		rwriter.close();

		FileWriter filewriter = new FileWriter(readfile.substring(0,readfile.indexOf("."))+"_KARGVA_mappedGenes.csv");
		BufferedWriter writer = new BufferedWriter(filewriter);
		writer.write("GeneIdx,KmerSNPHits,PercentGeneCovered,AverageKMerDepth\r\n");
		Collection<String> keysc = geneKmerMapping.keySet();
		ArrayList<String> keys = new ArrayList<String>(keysc);
		Collections.sort(keys);
		for (String key : keys)
		{
			AMRVariantGene ag = geneKmerMapping.get(key);
			Set<String> actualKmers = ag.kmerFreque.keySet();
			int totKmers = actualKmers.size();
			int totSNP = ag.kmerIsSNP.size();
			double kmerSNPHits = 0;
			double percCovered = 0;
			double kmerDepth = 0;
			for (String fk : actualKmers)
			{
				if (ag.kmerMapped.get(fk)!=null)
				{
					if (ag.kmerIsSNP.get(fk)!=null) {kmerSNPHits=kmerSNPHits+1d;}
					percCovered = percCovered + 1d;
					double dd = 0;
					if (ag.kmerMapped.get(fk)!=null) dd+=(double)(ag.kmerMapped.get(fk));
					dd=dd/(double)(ag.kmerFreque.get(fk));
					kmerDepth = kmerDepth + dd;
				}
			}
			kmerSNPHits = kmerSNPHits;
			percCovered = percCovered/(double)(totKmers);
			kmerDepth = kmerDepth/(double)(totKmers);
			if (percCovered>0.001f && kmerSNPHits>0)
			{
				writer.write(key+",");
				writer.write((int)(kmerSNPHits)+"/"+(int)(totSNP)+",");
				writer.write(100*percCovered+"%,");
				writer.write(kmerDepth+"\r\n");
			}
		}
		writer.close();
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - startTime;
		System.out.print("Reads and genes mapped in = "+elapsedTime/1000+" s\r\n");
		
		endTime = System.currentTimeMillis();
		elapsedTime = endTime - time0;
		System.out.print("Total time employed  = "+elapsedTime/1000+" s\r\n");
	}
}