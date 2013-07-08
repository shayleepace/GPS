package edu.utah.math5600.GPS;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class Satellite {

	private static Double[][] r3;
	private static Double[][] x;
	private static Data d = new Data();
	
	public static void main(String[] args) {
		//Matrices
		r3 = new Double[3][3];
		x = new Double[3][1];
		
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
//		r3[0][0] = Math.cos((2*d.pi*tv)/d.s);
//		r3[0][1] = -Math.sin((2*d.pi*tv)/d.s);
//		r3[0][2] = 0.0;
//		r3[1][0] = Math.sin((2*d.pi*tv)/d.s);
//		r3[1][1] = Math.cos((2*d.pi*tv)/d.s);
//		r3[1][2] = 0.0;
//		r3[2][0] = 0.0;
//		r3[2][1] = 0.0;
//		r3[2][2] = 1.0;
//		
//		x[0][0] = (d.r + 1372) * Math.cos(psi) * Math.cos(lambda);
//		x[1][0] = (d.r + 1372) * Math.cos(psi) * Math.sin(lambda);
//		x[2][0] = (d.r + 1372) * Math.sin(psi);
		
		
		
		//Fill matrices
		r3[0][0] = Math.cos((2*d.pi*tv)/d.s);
		r3[0][1] = -Math.sin((2*d.pi*tv)/d.s);
		r3[0][2] = 0.0;
		r3[1][0] = Math.sin((2*d.pi*tv)/d.s);
		r3[1][1] = Math.cos((2*d.pi*tv)/d.s);
		r3[1][2] = 0.0;
		r3[2][0] = 0.0;
		r3[2][1] = 0.0;
		r3[2][2] = 1.0;
		
		x[0][0] = (d.r + h) * Math.cos(psi) * Math.cos(lambda);
		x[1][0] = (d.r + h) * Math.cos(psi) * Math.sin(lambda);
		x[2][0] = (d.r + h) * Math.sin(psi);
		
		
		//testing matrices
//		double[][] a = new double[3][3];
//		double[][] b = new double[3][1];
//		
//		a[0][0] = 1; 
//		a[0][1] = 2;
//		a[0][2] = 3;
//		a[1][0] = 4;
//		a[1][1] = 5;
//		a[1][2] = 6;
//		a[2][0] = 1;
//		a[2][1] = 1;
//		a[2][2] = 1;
//		
//		b[0][0] = 1;
//		b[1][0] = 1;
//		b[2][0] = 2;
//		double[][] c = multMatrix(a,b);
//		System.out.println("matrix: ");
//		printMatrix(a);

//		System.out.println(c[0][0]);
//		System.out.println(c[1][0]);
//		System.out.println(c[2][0]);
		
		
		System.out.println("This should be u1 of satelite 15: " + d.u1[14]);
		System.out.println("This should be R: " + d.r);

	}
	
	
	public static double getPsi(double ns, double ad, double am, double as) {
		return 2 * d.pi * ns * ( (ad/360) + (am/21600) + (as/1296000) );
	}
	
	public static double getLambda(double ew, double bd, double bm, double bs) {
		return 2 * d.pi * ew * ( (bd/360) + (bm/21600) + (bs/1296000) );
	}
	
	public static double[][] multMatrix(double[][] a, double[][] b) {
		double[][] c = new double[a[0].length][b[1].length];
	    for (int i = 0; i < a[0].length; i++)
	    	for (int j = 0; j < b[1].length; j++)
	    		for (int k = 0; k < a[1].length; k++)
	    			c[i][j] += (a[i][k] * b[k][j]);
	    System.out.println("size: " + c[0].length + " " + c[1].length);
		return c;
	}
	
	public static void printMatrix(double[][] matrix) {
		System.out.println("size: " + matrix[0].length + " " + matrix[1].length);
	    for (int i = 0; i < matrix[0].length; i++) {
	    	for (int j = 0; j < matrix[1].length; j++) {
	    		System.out.print(matrix[i][j] + " ");
	    	}
	    	System.out.print("\n");
	    }
	}

}
