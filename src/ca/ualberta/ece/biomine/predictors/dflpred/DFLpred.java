package dflpred;

public class DFLpred {

	public static void main(String[] args) throws Exception {
		if (args.length != 2) {
			System.out.println("usage:");
			System.out.println("java -jar DFLpred.jar [input_fasta_file] [output_file]");
			return;
		}
		String input = args[0];
		String output = args[1];
		new Pred(input, output);
	}// main
}// class