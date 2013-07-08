package edu.utah.math5600.GPS;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class Data {
	
	public double pi, c, r, s;
	Double[] u1, u2, u3, v1, v2, v3, periodicity, altitude, phase;
	
	public Data() {
		u1 = new Double[24]; u2 = new Double[24]; u3 = new Double[24];
		v1 = new Double[24]; v2 = new Double[24]; v3 = new Double[24];
		periodicity = new Double[24]; altitude = new Double[24]; phase = new Double[24];
		
		Scanner dataInput;
		try {
			dataInput = new Scanner(new File("data.dat"));
			ArrayList<Double> dataList = new ArrayList<Double>();
			while(dataInput.hasNextLine()) {
				String dataLine[] = dataInput.nextLine().trim().split(" ");
				String data = dataLine[0].trim();
				dataList.add(Double.parseDouble(data));
			}
			
			pi = dataList.get(0);
			c = dataList.get(1);
			r = dataList.get(2);
			s = dataList.get(3);
			
			int j = 0;
			for(int i = 0; i < 24; i++) {
				u1[i] = dataList.get(4+j);
				u2[i] = dataList.get(5+j);
				u3[i] = dataList.get(6+j);
				v1[i] = dataList.get(7+j);
				v2[i] = dataList.get(8+j);
				v3[i] = dataList.get(9+j);
				periodicity[i] = dataList.get(10+j);
				altitude[i] = dataList.get(11+j);
				phase[i] = dataList.get(12+j);
				j += 9;
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
}
