package dflpred;
import java.util.ArrayList;
import java.util.List;


public class Fea_per_resi {
	public double num_iupG_0, std_iupL_val, avg_AAInd_AURR980118, avg_AAInd_PALJ810114;
	public Fea_per_resi(String seq, String iupG, ArrayList<Double> iUPredLSplit, double[][] aaInds, int pos) throws Exception{
		//define variables for window
		int winLen = 17, half = (winLen - 1)/2, protLen = seq.length(), winStart = pos - Math.min(pos, half),
				winEnd = pos + Math.min(protLen - 1 - pos, half), actualwinLen = winEnd - winStart + 1;
		double curIUPG0 = 0.0;
		double[][] aaInds_cur_win = new double[aaInds.length][actualwinLen];
		for (int cur = winStart; cur <= winEnd; cur++){
			//due with AA indexes
			char curAmino = seq.charAt(cur);
			switch (curAmino)
			{
			case 'a':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][0]; //for the second dimension 
				}
				break;
			case 'r':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][1]; //for the second dimension 
				}
				break;
			case 'n':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][2]; //for the second dimension 
				}
				break;
			case 'd':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][3]; //for the second dimension 
				}
				break;
			case 'c':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][4]; //for the second dimension 
				}
				break;
			case 'e':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][6]; //for the second dimension 
				}
				break;
			case 'q':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][5]; //for the second dimension 
				}
				break;
			case 'g':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][7]; //for the second dimension 
				}
				break;
			case 'h':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][8]; //for the second dimension 
				}
				break;
			case 'i':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][9]; //for the second dimension 
				}
				break;
			case 'l':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][10]; //for the second dimension 
				}
				break;
			case 'k':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][11]; //for the second dimension 
				}
				break;
			case 'm':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][12]; //for the second dimension 
				}
				break;
			case 'f':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][13]; //for the second dimension 
				}
				break;
			case 'p':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][14]; //for the second dimension 
				}
				break;
			case 's':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][15]; //for the second dimension 
				}
				break;
			case 't':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][16]; //for the second dimension 
				}
				break;
			case 'w':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][17]; //for the second dimension 
				}
				break;
			case 'y':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][18]; //for the second dimension 
				}
				break;
			case 'v':
				for(int i = 0; i < aaInds_cur_win.length; i++){ //for the first dimension
					aaInds_cur_win[i][cur - winStart] = aaInds[i][19]; //for the second dimension 
				}
				break;
			default:
				break;
			}// switch curAmino
			
			//due with IUPG
			char curIUPG = iupG.charAt(cur);
			if(curIUPG == '0')
				curIUPG0 ++;
		}//for cur in a window
		
		this.num_iupG_0 = curIUPG0/actualwinLen;
		
		List<Double> iupL_split_cur_win = iUPredLSplit.subList(winStart, winEnd + 1);
		this.std_iupL_val = std(iupL_split_cur_win);
		
		double[] aaAvg = meanValueOfArray(aaInds_cur_win);
		this.avg_AAInd_AURR980118 = aaAvg[0];
		this.avg_AAInd_PALJ810114 = aaAvg[1];
		
	}//construction
	
	private double std (List<Double> values)
	{
		double mean = meanValueOfList(values);
		double temp = 0.0;
		for(double val:values)
		{
			temp += (val - mean)*(val - mean);
		}
		return Math.sqrt(temp/(values.size()));
	}//double std
	
	private double meanValueOfList (List<Double> values)
	{
		double sum = 0.0;
		for(double val: values)
		{
			sum += val;		
		}
		return sum/values.size();
	}// private mean value of list
	
	private double[] meanValueOfArray(double[][] array){
		double[] val = new double[array.length];
		for(int m = 0; m < val.length; m++){ //for 1-dimension
			for(int n = 0; n < array[0].length; n++){ //for n, 2-dimension
				val[m] += array[m][n];
			}//for n, 2-dimension
			val[m] /= array[0].length;
		}//for m, 1-dimension
		return val;
	}// private mean value of array
	
}//class
