
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.molgenis.genotype.*;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.RandomAccessGenotypeData;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.variantFilter.VariantCombinedFilter;
import org.molgenis.genotype.variantFilter.VariantFilter;
import org.molgenis.genotype.variantFilter.VariantFilterableGenotypeDataDecorator;
import org.molgenis.genotype.variantFilter.VariantIdIncludeFilter;
import umcg.genetica.io.text.TextFile;


/**
 *
 * @author Patrick Deelen
 */
public class Compare2TriTyperDatasets {
	RandomAccessGenotypeData dataset1;
	RandomAccessGenotypeData dataset2;
	HashMap<String, Integer> dataset1SharedSampleMap;
	HashMap<String, Integer> dataset2SharedSampleMap;
	ArrayList<SharedSample> sharedSamplesList;
	Map<String, String> gte;

	public Compare2TriTyperDatasets(){}

	public Compare2TriTyperDatasets(String f1, String f2){

		try {
			dataset1 = new TriTyperGenotypeData(f1);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 1: " + f1);
		}

		try {
			System.out.println("Second");
			dataset2 = new TriTyperGenotypeData(f2);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 2: " + f2);

		}

	}

	public Compare2TriTyperDatasets(String f1, String f2, HashSet<String> SNPsToIncludeFirstDataset){

		try {
			VariantFilter variantFilter = new VariantIdIncludeFilter(SNPsToIncludeFirstDataset);
			dataset1 =  RandomAccessGenotypeDataReaderFormats.TRITYPER.createFilteredGenotypeData(f1, 1000, variantFilter, null);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 1: " + f1);

		}

		try {
			System.out.println("Second");
			dataset2 = new TriTyperGenotypeData(f2);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 2: " + f2);

		}

	}

	/**
	 * Compare using only SNPs from eQTL mapping results
	 * @param f1
	 * @param f2
	 * @param eQTLfile
	 */
	public Compare2TriTyperDatasets(String f1, String f2, String eQTLfile){
		HashSet eqtlSnps = null;
		try {
			TextFile eqtls = new TextFile(eQTLfile, false);
			eqtlSnps = new HashSet(eqtls.readAsSet(1, TextFile.tab));
			eqtls.close();
		} catch (IOException e) {
			System.out.println("Error reading eQTL file " + eQTLfile);
		}
		try {
			dataset1 = new TriTyperGenotypeData(f1);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 1: " + f1);
		}

		try {
			System.out.println("Second");
			VariantFilter variantFilter = new VariantIdIncludeFilter(eqtlSnps);
			dataset2 =  RandomAccessGenotypeDataReaderFormats.TRITYPER.createFilteredGenotypeData(f2, 1000, variantFilter, null);
		} catch (IOException ex) {
			System.err.println("Error reading dataset 2: " + f2);

		}
	}

	/**
	 * Get samples present in the genome-to-expression couplings file and the TriTyper datasets
	 * @param gteFile - first column contains samples from the first dataset, second column contains samples from the second dataset
	 * @throws IOException
	 */
	private void getSharedSamples(String gteFile) throws IOException {
		TextFile gte_coupling = new TextFile(gteFile, false);
		gte = gte_coupling.readAsHashMap(0, 1);
		gte_coupling.close();

		ArrayList<String> samples1 = new ArrayList<String>();
		ArrayList<String> samples2 = new ArrayList<String>();
		for (Sample s : dataset1.getSamples())
			samples1.add(s.getId());
		for (Sample s : dataset2.getSamples())
			samples2.add(s.getId());


		//check if the gte column order is correct. If it's the other way round, make a correct gte map
		for(String sample : samples1){
			if (gte.containsKey(sample)){
				break;
			}
			if (gte.containsValue(sample)){
				gte_coupling = new TextFile(gteFile, false);
				gte = gte_coupling.readAsHashMap(1, 0);
				gte_coupling.close();
				break;
			}
		}



		//Getting shared samples. Add both sample id from the first dataset and from the second dataset
		HashSet<String> sharedSamples = new HashSet<String>();

		for (String sample : samples1){
			String translation = gte.get(sample);
			if (translation != null){
				if (samples2.contains(translation)){
					sharedSamples.add(sample);
					sharedSamples.add(translation);
				}
			}
		}

		//Make a map mapping sample ids to sample position numbers in TriTyper
		System.out.println("\nNumber of shared samples: " + sharedSamples.size()/2 + "\n");

		dataset1SharedSampleMap = new HashMap<String, Integer>();
		dataset2SharedSampleMap = new HashMap<String, Integer>();

		int i = 0;
		for(Sample sample : dataset1.getSamples()){
			String id = sample.getId();
			if (sharedSamples.contains(id))
				dataset1SharedSampleMap.put(id, i);
			++i;
		}

		i = 0;
		for(Sample sample : dataset2.getSamples()){
			String id = sample.getId();
			if (sharedSamples.contains(id))
				dataset2SharedSampleMap.put(id, i);
			++i;
		}

		//List of sample ids from the first dataset
		sharedSamplesList = new ArrayList<SharedSample>();
		for (String id : dataset1SharedSampleMap.keySet()){
			sharedSamplesList.add(new SharedSample(id));
		}

	}

	private boolean isHeterogygous(Alleles alleles){
		char[] alls = alleles.getAllelesAsChars();
		if (alls[0] == alls[1])
			return false;
		return true;
	}

	private boolean passQC(GeneticVariant snp){
		if ((snp.getCallRate() > 0.5) && (snp.getMinorAlleleFrequency() > 0.05))
			return true;
		return false;
	}


	public void compare(){
		for(GeneticVariant dataset1Variant : dataset1){
			GeneticVariant dataset2Variant = dataset2.getSnpVariantByPos(dataset1Variant.getSequenceName(), dataset1Variant.getStartPos());

			if (dataset2Variant == null){
				//skipping non shared variants
				continue;
			}


			if ((! passQC(dataset1Variant)) || (! passQC(dataset2Variant))){
				//skipping SNPs not passing QC
				continue;
			}

			List<Alleles> dataset1VariantAlleles = dataset1Variant.getSampleVariants();
			List<Alleles> dataset2VariantAlleles = dataset2Variant.getSampleVariants();


			int nSNPs = 0;
			int numSharedSamples = 0, numDiscordantGenotypes = 0;

			for(SharedSample sharedSample : sharedSamplesList){

				//Get alleles for this shared sample
				Alleles x = dataset1VariantAlleles.get(dataset1SharedSampleMap.get(sharedSample.id));
				Alleles y = dataset2VariantAlleles.get(dataset2SharedSampleMap.get(gte.get(sharedSample.id)));

				if ((! isCalled(x) || (! isCalled(y)))){
					continue;
				}

				sharedSample.numShared++;
				numSharedSamples++;
				//sharedPerSample[sample]++;

				//Compare if same alleles. In this case AT == TA
				if(x.sameAlleles(y)){
					sharedSample.numConcordant++;
					//samePerSample[sample]++;
				}
				else{
					numDiscordantGenotypes++;
				}

				if (isHeterogygous(x))
					sharedSample.numHeterozygous++;

				if (isHeterogygous(y))
					sharedSample.numHeterozygousD2++;

				//sample++;
				nSNPs++;

			}

			/*
			float fractionDiscordant = (float)numDiscordantGenotypes/numSharedSamples;
			if (fractionDiscordant > 0.6)
				System.out.println(dataset1Variant.getSequenceName() + ":" + dataset1Variant.getStartPos() + "\t" + fractionDiscordant);
			*/

			//System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + corr + "\t" + countIdentical + "\t" + countMismatch + "\t" + hetero);
			//System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + countIdentical + "\t" + countMismatch);

			if (nSNPs % 1000 == 0){
				System.out.println("Processed " + nSNPs + " shared SNPs");
			}
		}

	int num = 0;
	for (SharedSample sample : sharedSamplesList){
		System.out.println(sample.toFullString());
		num++;
	}
	}

	/**
	 * Gets and integer for each genotype:
	 *   if both samples homozygous for same allele: {0,0}
	 *   if both samples homozygous for different alleles: {0,2}
	 *   if both samples are heterozygous: {1,1}
	 * @param alleles1
	 * @param alleles2
	 * @return
	 */
	private int[] getIntGenotype(char[] alleles1, char alleles2[]){
		int[] genotype = new int[] {-1, -1};
		//if any of the samples is not called, exit:
		if (((alleles2[0] == '0') || (alleles2[1] == '0')) || ((alleles1[0] == '0') || (alleles1[1] == '0')))
			return genotype;

		if ((alleles1[0] == alleles1[1]) && (alleles1[0] == alleles2[0]) && (alleles2[0] == alleles2[1])){ //if both samples are homozygous by the same allele
			genotype[0] = 0;
			genotype[1] = 0;
			return genotype;
		}

		HashSet<Character> all_alleles = new HashSet<Character>();
		all_alleles.add(alleles1[0]);
		all_alleles.add(alleles1[1]);
		all_alleles.add(alleles2[0]);
		all_alleles.add(alleles2[1]);
		ArrayList<Character> allele_list = new ArrayList<Character>(all_alleles);
		if (all_alleles.size() == 2){
			if (alleles1[0] == allele_list.get(0).charValue()){
				genotype[0] = (alleles1[0] != alleles1[1]) ? 1 : 0;
				genotype[1] = (alleles2[0] != alleles2[1]) ? 1 : 2;
			}
			else if (alleles1[0] == allele_list.get(1).charValue()){
				genotype[0] = (alleles1[0] != alleles1[1]) ? 1 : 2;
				genotype[1] = (alleles2[0] != alleles2[1]) ? 1 : 0;
			}
		}
		return genotype;
	}

	/**
	 * Calculates correlation for a SNP between genotypes from the two datasets
	 * @param dataset1VariantAlleles
	 * @param dataset2VariantAlleles
	 * @return
	 */
	private double calculateCorrelationSharedSamples(List<Alleles> dataset1VariantAlleles, List<Alleles> dataset2VariantAlleles){
		int[] alleles1;
		int[] alleles2;
		ArrayList<Integer> allelesList1 = new ArrayList<Integer>();
		ArrayList<Integer> allelesList2= new ArrayList<Integer>();
		int cnt = 0;
		boolean empty = true;
		for(SharedSample sharedSample : sharedSamplesList){

			//Get alleles for this shared sample
			Alleles x = dataset1VariantAlleles.get(dataset1SharedSampleMap.get(sharedSample.id));
			Alleles y = dataset2VariantAlleles.get(dataset2SharedSampleMap.get(gte.get(sharedSample.id)));

			if ((! isCalled(x) || (! isCalled(y)))){
				continue;
			}
			empty = false;

			//get integer genotypes as in Lude's method
			int[] genotypes = getIntGenotype(x.getAllelesAsChars(), y.getAllelesAsChars());
			allelesList1.add(genotypes[0]);
			allelesList2.add(genotypes[1]);

			cnt++;
		}

		if (! empty){
			alleles1 = ArrayUtils.toPrimitive(allelesList1.toArray(new Integer[cnt]));
			alleles2 = ArrayUtils.toPrimitive(allelesList2.toArray(new Integer[cnt]));
			//System.out.println("\t" + JSci.maths.ArrayMath.correlation(alleles1, alleles2));
			return JSci.maths.ArrayMath.correlation(alleles1, alleles2);
		}
		else
			return -1;

	}

	private boolean isCalled(Alleles al){
		if ((al.getAllelesAsChars()[0] == '0') || (al.getAllelesAsChars()[1] == '0'))
			return false;
		return true;
	}


	/**
	 * Compares genotype concordance only for highly correlated SNPs:
	 * to get outliers only based on SNPs with high concordance for most of the samples
 	 */

	public void compareHighlyCorrelatedSNPs(){
		//System.out.println("variantId\tcountIdentical\tcountMismatch");
		double corrThreshold = 0.9;
		int numShared = sharedSamplesList.size();
		int[] sharedPerSample = new int[numShared];
		int[] samePerSample = new int[numShared];
		int numPass = 0, numPassCorr = 0;

		System.out.println("Comparing SNPs with correlation higher than " + corrThreshold);

		for(GeneticVariant dataset1Variant : dataset1){
			GeneticVariant dataset2Variant = dataset2.getSnpVariantByPos(dataset1Variant.getSequenceName(), dataset1Variant.getStartPos());

			if (dataset2Variant == null){
				//skipping non shared variants
				continue;
			}


			/*if ((! passQC(dataset1Variant)) || (! passQC(dataset2Variant))){
				//skipping SNPs not passing QC'
				continue;
			}*/

			if (! passQC(dataset1Variant)) {
				//skipping SNPs not passing QC only in the first dataset
				continue;
			}
			numPass++;
			List<Alleles> dataset1VariantAlleles = dataset1Variant.getSampleVariants();
			List<Alleles> dataset2VariantAlleles = dataset2Variant.getSampleVariants();

			int countIdentical = 0;
			int countMismatch = 0;
			int sample = 0;
			int hetero = 0;
			double corr = calculateCorrelationSharedSamples(dataset1VariantAlleles, dataset2VariantAlleles);
			if (corr > 0.9){
				numPassCorr++;
				for(SharedSample sharedSample : sharedSamplesList){

					//Get alleles for this shared sample
					Alleles x = dataset1VariantAlleles.get(dataset1SharedSampleMap.get(sharedSample.id));
					Alleles y = dataset2VariantAlleles.get(dataset2SharedSampleMap.get(gte.get(sharedSample.id)));

					if ((! isCalled(x) || (! isCalled(y)))){
						continue;
					}

					sharedSample.numShared++;
					//sharedPerSample[sample]++;

					//Compare if same alleles. In this case AT == TA
					if(x.sameAlleles(y)){
						sharedSample.numConcordant++;
						countIdentical++;
					}
					else{
						countMismatch++;
					}
					if (isHeterogygous(x)){
						sharedSample.numHeterozygous++;
						hetero++;
					}
					if (isHeterogygous(y)){
						sharedSample.numHeterozygousD2++;

					}
					sample++;

				}
				//System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + corr + "\t" + countIdentical + "\t" + countMismatch + "\t" + hetero);
				//System.out.println(dataset1Variant.getPrimaryVariantId() + "\t" + countIdentical + "\t" + countMismatch);
			if (numPassCorr % 1000 == 0){
				System.out.println("Processed " + numPassCorr + " shared SNPs");
			}
			}
		}

		for (SharedSample sample : sharedSamplesList){
			System.out.println(sample.toFullString());

		}
		System.out.println("number of SNPs " + numPass);
	}

	private void getAllSNPPassQC(){
		for(GeneticVariant var : dataset2){
			if (passQC(var)){
				//System.out.println(var.getPrimaryVariantId());
				GeneticVariant var2 = dataset1.getSnpVariantByPos(var.getSequenceName(), var.getStartPos());
				if (var2 != null){
					System.out.println(var.getSequenceName() +" " + var.getStartPos());
				}
			}
		}
	}
	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws IOException {

		/*String f2 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/LL/genotypes/SNVMix-TriTyper/",
				f1="/Users/dashazhernakova/Documents/UMCG/data/BBMRI/LL/genotypes/TriTyper_pos/";

		Compare2TriTyperDatasets compare = new Compare2TriTyperDatasets(f1, f2);
		compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/LL/expression_table/gte_coupling.txt");
		//String f1 = args[1], f2 = args[2];
		*/

		/*String f1 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/CODAM-imputed-20130828-trityper/",
				f2 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/SNVMix-TriTyper/";

		HashSet<String> includedSNPs = new HashSet<String>();
		includedSNPs.add("rs6517457");
		includedSNPs.add("rs1494558");

		Compare2TriTyperDatasets compare = new Compare2TriTyperDatasets(f1, f2, includedSNPs);
		compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/gte_coupling_codam.txt");
		*/
		/*
		String f2 = "/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/SNVMix-TriTyper/",
				f1="/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM/genotypes/CODAM-imputed-20130828-trityper/";

		Compare2TriTyperDatasets compare = new Compare2TriTyperDatasets(f1, f2);
		compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/BBMRI/CODAM//gte_coupling_codam.txt");
		 */
		for (String a : args){
			System.out.println(a);
		}
		Compare2TriTyperDatasets compare;
		if (args.length < 3){
			String f2 = "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/GBR/rna-seq/SNVMix-TriTyper/",
					f1="/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/GBR/dna-seq/TriTyper/";

			//compare = new Compare2TriTyperDatasets(f1, f2);
			compare = new Compare2TriTyperDatasets(f1,f2);
			compare.getSharedSamples("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table_old/GBR/gte_HGtoERR_GBR.txt");
			compare.compare();
		}
		else if (args.length == 3){
			compare = new Compare2TriTyperDatasets(args[0], args[1]);
			compare.getSharedSamples(args[2]);
			compare.compareHighlyCorrelatedSNPs();
			//compare.compare();
		}
		else if (args.length == 4){
			compare = new Compare2TriTyperDatasets(args[0], args[1], args[3]);
			compare.getSharedSamples(args[2]);
			compare.compareHighlyCorrelatedSNPs();
			//compare.compare();
		}
		else
			System.out.println("Wrong number of arguments: " + args.length);
		//compare.getAllSNPPassQC();




	}
}


