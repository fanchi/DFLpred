package dflpred;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class Helper {
	static public String readLine(String strFileName, int line) {
		BufferedReader buffer = null;
		try {
			buffer = new BufferedReader(new FileReader(strFileName));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		int count = 0;
		String currentLine = null;
		try {
			while ((currentLine = buffer.readLine()) != null && count < line) {
				count++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		try {
			buffer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return currentLine;
	}// public static String ReadLine

	static public double sigmoid(double z) {
		return 1.0 / (Math.exp((-1) * z) + 1.0);
	}// static public double sigmoid

	static public void recursiveDelete(File file) {
		if (!file.exists())
			return;
		if (file.isDirectory()) {
			for (File f : file.listFiles()) {
				recursiveDelete(f);
			}
		}
		file.delete();
	}//static public void recursiveDelete

}// class Helper
