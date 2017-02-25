package ca.ualberta.ece.biomine.predictors.dflpred;

import java.io.File;

public class DFLpred {

	public static void main(String[] args) {
		if (args.length != 2) {
			System.out.println("usage:");
			System.out.println("java -jar DFLpred.jar [input_fasta_file] [output_file]");
			return;
		}
		String input = args[0];
		String output = args[1];
		try {
			new Pred(input, output);
		} catch (Exception e) {
			System.out.println("IUPred was not run successfully. Please make sure the \"myIUPRED\" folder and \"AAind_to_use\" file are "
					+ "in the same direcotry of DFLpred.jar. \"myIUPRED\" was slightly modified from \"IUPred\", which can be complied on "
					+ "any ANSI C compiler e.g. GUN C compiler.");
			Helper.recursiveDelete(new File(output));
		}
	}// main
}// class