package ca.ualberta.ece.biomine.predictors.dflpred;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;

public class Pred {
	public Pred(String inputFileName, String outputFileName) {
		// define file names and paths
		String timeStamp = new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String singleFastaFilePath = "./temp" + timeStamp + "/";
		String iupredPath = "./myIUPRED";
		String aaIndFileName = "AAind_to_use";
		new File(singleFastaFilePath).mkdir();
		double[] coeffs = { -1.0929, 6.5992, -4.5844, 3.3256, -0.9237 };
		double thres = 0.18;
		ArrayList<String> identifiers = new ArrayList<String>();
		try {
			BufferedReader buffer = new BufferedReader(new FileReader(inputFileName));
			String currentLine;
			int lineCount = 0;
			while ((currentLine = buffer.readLine()) != null) {
				if (currentLine.substring(0, 1).equals(">")) {
					identifiers.add(currentLine.substring(1));
					// write to single fasta file
					String singleFastaFileName = singleFastaFilePath + currentLine.substring(1) + ".fasta";
					File singleFastaFile = new File(singleFastaFileName);
					FileOutputStream singleFastaOutputStream = new FileOutputStream(singleFastaFile);
					OutputStreamWriter singleFastaFileOsw = new OutputStreamWriter(singleFastaOutputStream);
					singleFastaFileOsw.write(currentLine + "\n");
					String seq;
					int currentLineCount = lineCount;
					while ((seq = Helper.readLine(inputFileName, ++ currentLineCount)) != null && !seq.substring(0, 1).equals(">")) {
						singleFastaFileOsw.write(seq);
					}// while seq lines for current identifier
					singleFastaFileOsw.close();
				}// if identifier line
				lineCount ++;
			}// while currentLine in input file
			buffer.close();
		} catch (IOException e) {
			System.out.println("FASTA file was not read correctly");
			Helper.recursiveDelete(new File(singleFastaFilePath));
		}
		// run IUPred
		if (identifiers.size() == 0) {
			System.out.println("no fasta sequences or no identifiers");
			return;
		}
		// prepare output file
		File outputFile = null;
		FileOutputStream outputFileStream = null;
		OutputStreamWriter outputFileOsw = null;
		try {
			outputFile = new File(outputFileName);
			outputFileStream = new FileOutputStream(outputFile);
			outputFileOsw = new OutputStreamWriter(outputFileStream);
		} catch (IOException e) {
			System.out.println("creating output file fail");
			Helper.recursiveDelete(new File(singleFastaFilePath));
		}
		for (String iden : identifiers) {
			try {
				outputFileOsw.write(">" + iden + "\n");
			} catch (IOException e) {
				System.out.println("writing to file error");
				Helper.recursiveDelete(new File(singleFastaFilePath));
			}
			String curIUPL = null, curIUPG = null;
			ProcessBuilder pb = new ProcessBuilder("./myiupred", "." + singleFastaFilePath + iden + ".fasta", "glob");
			pb.directory(new File(iupredPath));
			pb.redirectErrorStream(true);
			try {
				Process process = pb.start();
				BufferedReader inStreamReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
				String runStatuString = null;
				while ((runStatuString = inStreamReader.readLine()) != null) {
					curIUPG = runStatuString;
				}
			} catch (Exception e) {
				System.out.println("run IUPred glob error");
				Helper.recursiveDelete(new File(singleFastaFilePath));
			}
			pb = new ProcessBuilder("./myiupred", "." + singleFastaFilePath + iden + ".fasta", "long");
			pb.directory(new File(iupredPath));
			pb.redirectErrorStream(true);
			try {
				Process process = pb.start();
				BufferedReader inStreamReader = new BufferedReader(new InputStreamReader(process.getInputStream()));
				String runStatuString = null;
				while ((runStatuString = inStreamReader.readLine()) != null) {
					curIUPL = runStatuString;
				}
			} catch (Exception e) {
				System.out.println("run IUPred long error");
				Helper.recursiveDelete(new File(singleFastaFilePath));
			}
			// get [][]aainds
			double[][] aainds = new double[2][20];
			BufferedReader buffer = null;
			try {
				buffer = new BufferedReader(new FileReader(aaIndFileName));
				String currentLine;
				int i = 0;
				while ((currentLine = buffer.readLine()) != null) {
					String[] splited = currentLine.split(" ");
					for (int k = 1; k < splited.length; k ++) {
						aainds[i][k - 1] = Double.valueOf(splited[k]);
					}
					i ++;
				}// while current line in aaind file
			} catch (IOException e) {
				System.out.println("read AA index error");
				Helper.recursiveDelete(new File(singleFastaFilePath));
			} finally {
				try {
					buffer.close();
				} catch (Exception e2) {
					// do nothing
				}
			}
			// get sequence for current identifier
			String curSeq = Helper.readLine(singleFastaFilePath + iden + ".fasta", 1);
			curSeq = curSeq.toLowerCase();
			char[] chCurSeq = curSeq.toCharArray();
			curIUPG = curIUPG.replaceAll(",", "");
			ArrayList<Double> iUPredLSplit = new ArrayList<Double>();
			String[] splited = curIUPL.split(",");
			for (int i = 0; i < splited.length; i ++) {
				iUPredLSplit.add(Double.valueOf(splited[i]));
			}
			String[] scores = new String[curSeq.length()];
			// get features per residue
			for (int resNum = 0; resNum < curSeq.length(); resNum ++) {
				FeaPerResidue feaCurRes = new FeaPerResidue(curSeq, curIUPG, iUPredLSplit, aainds, resNum);
				double z = feaCurRes.numIupG0 * coeffs[0] + feaCurRes.stdIupLVal * coeffs[1] + feaCurRes.avgAAIndAURR980118 * coeffs[2]
						+ feaCurRes.avgAAIndPALJ810114 * coeffs[3] + coeffs[4];
				double score = Helper.sigmoid(z);
				if (score >= thres) {
					// binary[resNum] = '1';
					chCurSeq[resNum] = Character.toUpperCase(chCurSeq[resNum]);
				}
				scores[resNum] = toDecimal(score);
			}// for resNum
			String strCurSeq = String.valueOf(chCurSeq);
			// str_binary = String.valueOf(binary);
			try {
				outputFileOsw.write(strCurSeq + "\n");
				for (int resNum = 0; resNum < curSeq.length() - 1; resNum ++)
					outputFileOsw.write(scores[resNum] + ",");
				outputFileOsw.write(scores[curSeq.length() - 1]);
				outputFileOsw.write("\n");
			} catch (Exception e) {
				System.out.println("writing top result file error");
			} 
		}// for every sequence in input file
		try {
			outputFileOsw.close();
		} catch (Exception e) {
			// do nothing
		}
		Helper.recursiveDelete(new File(singleFastaFilePath));
	}// DFLpred construction

	private String toDecimal(double val) {
		DecimalFormat df = new DecimalFormat("#.###");
		return df.format(val);
	}// ToDecimal

}// class
