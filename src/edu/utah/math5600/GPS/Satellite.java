package edu.utah.math5600.GPS;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class Satellite {

	private static double[][] r3;
	private static double[][] x;
	
	public static void main(String[] args) {
		//Get Data from data.dat
		Data.getData();
		
		//Matrices
		r3 = new double[3][3];
		x = new double[3][1];
		
		//Variables for input data
		double tv = 0.0, ad = 0.0, am = 0.0, as = 0.0, bd = 0.0, bm = 0.0, bs = 0.0, h = 0.0;
		int ns = 0, ew = 0;
		double psi, lambda;
		
		//Get vehicle input
		Scanner vehicleInput;
		try {
			vehicleInput = new Scanner(new File("vehicle.dat"));
			
			//While input has a another line to read
			while(vehicleInput.hasNextLine()) {
				String data[] = vehicleInput.nextLine().split(" ");
				
				tv = Double.parseDouble(data[0]);
				ad = Double.parseDouble(data[1]);
				am = Double.parseDouble(data[2]);
				as = Double.parseDouble(data[3]);
				ns = Integer.parseInt(data[4]);
				bd = Double.parseDouble(data[5]);
				bm = Double.parseDouble(data[6]);
				bs = Double.parseDouble(data[7]);
				ew = Integer.parseInt(data[8]);
				h = Double.parseDouble(data[9]);
			}
			
//			System.out.println("This should be h: " + h);
//			
//			System.out.println("This should be bd: " + bd);
		
		
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		psi = getPsi(ns, ad, am, as);
		lambda = getLambda(ew, bd, bm, bs);
		
		//Test 
//		psi = getPsi(1, 40, 45, 55);
//		lambda = getLambda(-1, 111, 50, 58);
//		r3[0][0] = Math.cos((2*Data.pi*tv)/Data.s);
//		r3[0][1] = -Math.sin((2*Data.pi*tv)/Data.s);
//		r3[0][2] = 0.0;
//		r3[1][0] = Math.sin((2*Data.pi*tv)/Data.s);
//		r3[1][1] = Math.cos((2*Data.pi*tv)/Data.s);
//		r3[1][2] = 0.0;
//		r3[2][0] = 0.0;
//		r3[2][1] = 0.0;
//		r3[2][2] = 1.0;
//		
//		x[0][0] = (Data.r + 1372) * Math.cos(psi) * Math.cos(lambda);
//		x[1][0] = (Data.r + 1372) * Math.cos(psi) * Math.sin(lambda);
//		x[2][0] = (Data.r + 1372) * Math.sin(psi);
		
		
		
		//Fill matrices
		r3[0][0] = Math.cos((2*Data.pi*tv)/Data.s);
		r3[0][1] = -Math.sin((2*Data.pi*tv)/Data.s);
		r3[0][2] = 0.0;
		r3[1][0] = Math.sin((2*Data.pi*tv)/Data.s);
		r3[1][1] = Math.cos((2*Data.pi*tv)/Data.s);
		r3[1][2] = 0.0;
		r3[2][0] = 0.0;
		r3[2][1] = 0.0;
		r3[2][2] = 1.0;
		
		x[0][0] = (Data.r + h) * Math.cos(psi) * Math.cos(lambda);
		x[1][0] = (Data.r + h) * Math.cos(psi) * Math.sin(lambda);
		x[2][0] = (Data.r + h) * Math.sin(psi);
		

		System.out.println("This should be u1 of satelite 15: " + Data.u1[14]);
		System.out.println("This should be R: " + Data.r);

	}
	
	
	public static double getPsi(double ns, double ad, double am, double as) {
		return 2 * Data.pi * ns * ( (ad/360) + (am/21600) + (as/1296000) );
	}
	
	public static double getLambda(double ew, double bd, double bm, double bs) {
		return 2 * Data.pi * ew * ( (bd/360) + (bm/21600) + (bs/1296000) );
	}
}
