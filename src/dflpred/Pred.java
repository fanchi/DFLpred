package dflpred;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;

public class Pred {
	public Pred(String inputFileName, String outputFileName)
			throws Exception {
		// define file names and paths
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss")
				.format(Calendar.getInstance().getTime());
		String single_fasta_file_path = "./temp" + timeStamp + "/";
		String iupredPath = "./myIUPRED";
		String aaIndFileName = "AAind_to_use";
		new File(single_fasta_file_path).mkdir();
		double[] coeffs = { -1.0929, 6.5992, -4.5844, 3.3256, -0.9237 };
		double thres = 0.18;
		// generate input files for IUPred
		BufferedReader buffer = new BufferedReader(
				new FileReader(inputFileName));
		String currentLine;
		int lineCount = 0;
		ArrayList<String> identifiers = new ArrayList<String>();
		while ((currentLine = buffer.readLine()) != null) {
			if (currentLine.substring(0, 1).equals(">")) {
				identifiers.add(currentLine.substring(1));
				// write to single fasta file
				String singleFastaFileName = single_fasta_file_path
						+ currentLine.substring(1) + ".fasta";
				File singleFastaFile = new File(singleFastaFileName);
				FileOutputStream singleFastaOutputStream = new FileOutputStream(
						singleFastaFile);
				OutputStreamWriter singleFastaFile_osw = new OutputStreamWriter(
						singleFastaOutputStream);
				singleFastaFile_osw.write(currentLine + "\n");
				String seq;
				int currentLineCount = lineCount;
				while ((seq = Helper
						.readLine(inputFileName, ++currentLineCount)) != null
						&& !seq.substring(0, 1).equals(">")) {
					singleFastaFile_osw.write(seq);
				}// while seq lines for current identifier
				singleFastaFile_osw.close();
			}// if identifier line
			lineCount++;
		}// while currentLine in input file
		buffer.close();

		// run IUPred
		if (identifiers.size() == 0) {
			System.out.println("no fasta sequences or no identifiers");
			return;
		}

		// prepare output file
		File outputFile = new File(outputFileName);
		FileOutputStream outputFileStream = new FileOutputStream(outputFile);
		OutputStreamWriter outputFile_osw = new OutputStreamWriter(
				outputFileStream);

		for (String iden : identifiers) {
			outputFile_osw.write(">" + iden + "\n");
			String curIUPL = null, curIUPG = null;
			ProcessBuilder pb = new ProcessBuilder("./myiupred", "."
					+ single_fasta_file_path + iden + ".fasta", "glob");
			pb.directory(new File(iupredPath));
			pb.redirectErrorStream(true);
			Process process = pb.start();
			BufferedReader inStreamReader = new BufferedReader(
					new InputStreamReader(process.getInputStream()));
			String runStatuString = null;
			while ((runStatuString = inStreamReader.readLine()) != null) {
				curIUPG = runStatuString;
			}

			pb = new ProcessBuilder("./myiupred", "." + single_fasta_file_path
					+ iden + ".fasta", "long");
			pb.directory(new File(iupredPath));
			pb.redirectErrorStream(true);
			process = pb.start();
			inStreamReader = new BufferedReader(new InputStreamReader(
					process.getInputStream()));
			while ((runStatuString = inStreamReader.readLine()) != null) {
				curIUPL = runStatuString;
			}

			// get [][]aainds
			buffer = new BufferedReader(new FileReader(aaIndFileName));
			double[][] aainds = new double[2][20];
			int i = 0;
			while ((currentLine = buffer.readLine()) != null) {
				String[] splited = currentLine.split(" ");
				for (int k = 1; k < splited.length; k++) {
					aainds[i][k - 1] = Double.valueOf(splited[k]);
				}
				i++;
			}// while current line in aaind file

			// get sequence for current identifier
			String curSeq = Helper.readLine(single_fasta_file_path + iden
					+ ".fasta", 1);
			// char[] binary = new char[cur_seq.length()];
			// Arrays.fill(binary, '0');
			curSeq = curSeq.toLowerCase();
			char[] ch_cur_seq = curSeq.toCharArray();

			curIUPG = curIUPG.replaceAll(",", "");
			ArrayList<Double> iUPred_L_split = new ArrayList<Double>();
			String[] splited = curIUPL.split(",");
			for (i = 0; i < splited.length; i++)
				iUPred_L_split.add(Double.valueOf(splited[i]));

			String[] scores = new String[curSeq.length()];
			// get features per residue
			for (int resNum = 0; resNum < curSeq.length(); resNum++) {
				Fea_per_resi fea_cur_res = new Fea_per_resi(curSeq, curIUPG,
						iUPred_L_split, aainds, resNum);
				double z = fea_cur_res.num_iupG_0 * coeffs[0]
						+ fea_cur_res.std_iupL_val * coeffs[1]
						+ fea_cur_res.avg_AAInd_AURR980118 * coeffs[2]
						+ fea_cur_res.avg_AAInd_PALJ810114 * coeffs[3]
						+ coeffs[4];
				double score = Helper.sigmoid(z);
				if (score >= thres) {
					// binary[resNum] = '1';
					ch_cur_seq[resNum] = Character
							.toUpperCase(ch_cur_seq[resNum]);
				}
				scores[resNum] = toDecimal(score);
				// outputFile_osw.write(ToDecimal(score) + ",");
			}// for resNum
			String str_cur_seq = String.valueOf(ch_cur_seq);
			// str_binary = String.valueOf(binary);
			outputFile_osw.write(str_cur_seq + "\n");
			for (int resNum = 0; resNum < curSeq.length() - 1; resNum++)
				outputFile_osw.write(scores[resNum] + ",");
			outputFile_osw.write(scores[curSeq.length() - 1]);
			outputFile_osw.write("\n");
			// outputFile_osw.write("\n" + str_binary + "\n");
		}// for every sequence in input file
		outputFile_osw.close();

		// delete single fasta files
		Helper.recursiveDelete(new File(single_fasta_file_path));

	}// DFLpred construction

	private String toDecimal(double val) {
		DecimalFormat df = new DecimalFormat("#.###");
		return df.format(val);
	}// ToDecimal

}// class
