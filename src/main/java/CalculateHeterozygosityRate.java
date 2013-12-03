import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.trityper.TriTyperGenotypeData;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.util.trityper.reader.TextFile;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class CalculateHeterozygosityRate {
	TriTyperGenotypeData trityperData;
	ArrayList<TriTyperSample> samples;
	static final int[][] nonpseudoautosomal = new int[][]{{1,60000},{2699521,154931043},{155260561,155270560}};
	public CalculateHeterozygosityRate(String path){
		try {
			trityperData = new TriTyperGenotypeData(path);

			//initialize samples list
			samples = new ArrayList<TriTyperSample>();
			for (Sample s : trityperData.getSamples()){
				TriTyperSample sample = new TriTyperSample(s.getId());
				samples.add(sample);
			}
		} catch (IOException ex) {
			System.err.println("Error reading dataset: " + path);
		}
	}

	private boolean isHeterogygous(Alleles alleles){
		char[] alls = alleles.getAllelesAsChars();
		if (alls[0] == alls[1])
			return false;
		return true;
	}

	public void runNonpseudoautosomal(){

		int start, end;
		for (int i = 0; i < nonpseudoautosomal.length; i++){
			start = nonpseudoautosomal[i][0];
			end = nonpseudoautosomal[i][1];

			for(GeneticVariant variant : trityperData.getVariantsByRange("X", start, end)){

				List<Alleles> variantAlleles = variant.getSampleVariants();

				TriTyperSample sample;
				int sampleCnt = 0;
				for (Sample s : trityperData.getSamples()){
					Alleles alleles = variantAlleles.get(sampleCnt);
					sample = samples.get(sampleCnt);
					sampleCnt++;

					if ((alleles.getAllelesAsChars()[0] == '0') || (alleles.getAllelesAsChars()[0] == '0'))
						continue;
					sample.numSNPs++;
					if (isHeterogygous(alleles)){
						sample.numHeterozygous++;
					}

				}
			}
		}
		for (TriTyperSample s : samples){
			System.out.println(s.id + "\t" + s.numSNPs + "\t" + s.numHeterozygous);
		}
	}

	public void run(String outFileName) throws IOException {
		TextFile out = new TextFile(outFileName, true);
		//System.out.println("sample\tnumSNPs\tnumHetero\tfractionHetero");
		out.write("sample\tnumSNPs\tnumHetero\tfractionHetero\n");
		for(GeneticVariant variant : trityperData){
			List<Alleles> variantAlleles = variant.getSampleVariants();

			TriTyperSample sample;
			int sampleCnt = 0;
			for (Sample s : trityperData.getSamples()){
				Alleles alleles = variantAlleles.get(sampleCnt);
				sample = samples.get(sampleCnt);
				sampleCnt++;

				if ((alleles.getAllelesAsChars()[0] == '0') || (alleles.getAllelesAsChars()[0] == '0'))
					continue;
				sample.numSNPs++;
				if (isHeterogygous(alleles)){
					sample.numHeterozygous++;
				}

			}
		}

		for (TriTyperSample s : samples){
			//System.out.println(s.id + "\t" + s.numSNPs + "\t" + s.numHeterozygous + "\t" + (float)s.numHeterozygous/s.numSNPs);
			out.write(s.id + "\t" + s.numSNPs + "\t" + s.numHeterozygous + "\t" + (float)s.numHeterozygous/s.numSNPs + "\n");
		}
		out.close();
	}
	public static void main(String[] args) throws IOException {
		//CalculateHeterozygosityRate hetRate = new CalculateHeterozygosityRate("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/genotypes/FIN/dna-seq/TriTyper/");
		CalculateHeterozygosityRate hetRate = new CalculateHeterozygosityRate(args[0]);
		hetRate.run(args[1]);


	}
}
