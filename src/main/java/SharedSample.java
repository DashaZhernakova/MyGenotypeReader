import java.util.HashMap;
import java.util.HashSet;

/**
 * Created with IntelliJ IDEA.
 * User: dashazhernakova
 * Date: 09.09.13
 * Time: 11:58
 * To change this template use File | Settings | File Templates.
 */
public class SharedSample extends TriTyperSample{
	public int numShared;
	public int numConcordant;
	public int numHeterozygousD2;
	public HashSet<String> discordantSNPs;

	public SharedSample(String sampleId){
		super(sampleId);
		numShared = 0;
		numConcordant = 0;
		numHeterozygousD2 = 0;
		discordantSNPs = new HashSet<String>();

	}

	public String toString(){

		return id + " (" + numShared + " " + String.format("%.2f", 100*numConcordant/numShared) + "% " + String.format("%.2f", 100*numHeterozygous/numShared) + "% "+ String.format("%.2f", 100*numHeterozygousD2/numShared + "%)");
	}
	public String toFullString(){
		//float conc = numConcordant/numShared;
		//float het = numHeterozygous/numShared;
		return id + "\t" + numShared + "\t" + numConcordant + "\t" + numHeterozygous + "\t" + numHeterozygousD2;

	}

}
