package edu.utah.math5600.GPS;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;


public class Satellite {

	private static double[][] r3;
	private static double[][] x;
	private static double[][] xv;
	
	
	public static void main(String[] args) {
		//Get Data from data.dat
		Data.getData();
		
		//Matrices
		r3 = new double[3][3];
		x = new double[3][1];
		xv = new double[3][1];
		
		//Variables for input data
		double tv = 0.0, ad = 0.0, am = 0.0, as = 0.0, bd = 0.0, bm = 0.0, bs = 0.0, h = 0.0;
		int ns = 0, ew = 0;
		double psi, lambda;
		double xv1 = 0.0, xv2 = 0.0, xv3 = 0.0;
		double ts = 0.0;
		
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
		
		//Computes the position of the vehicle
		xv1 = r3[0][0]*x[0][0]+r3[0][1]*x[1][0]+r3[0][2]*x[2][0];
		xv2 = r3[1][0]*x[0][0]+r3[1][1]*x[1][0]+r3[1][2]*x[2][0];
		xv3 = r3[2][0]*x[0][0]+r3[2][1]*x[1][0]+r3[2][2]*x[2][0];
		
		//Matrix position of vehicle
		xv[0][0] = xv1;
		xv[1][0] = xv2;
		xv[2][0] = xv3;
		
		
		double[] xstv = new double[Data.altitude.length];
		double[] ystv = new double[Data.altitude.length];
		double[] zstv = new double[Data.altitude.length];
		double[] t0 = new double[Data.altitude.length];
		//loop through satellites
		for(int i = 0; i < Data.altitude.length; i++) {
			//Computes position of each satellite at time tv
			xstv[i] = (Data.r+Data.altitude[i])*(Data.u1[i]*Math.cos(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i])+Data.v1[i]*Math.sin(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i]));
			ystv[i] = (Data.r+Data.altitude[i])*(Data.u2[i]*Math.cos(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i])+Data.v2[i]*Math.sin(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i]));
			zstv[i] = (Data.r+Data.altitude[i])*(Data.u3[i]*Math.cos(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i])+Data.v3[i]*Math.sin(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i]));
		
			//Find t_0 for each satellite (start of Newton's Method
			t0[i] = tv - Math.sqrt( Math.pow((xstv[i] - xv1), 2) + Math.pow((ystv[i] - xv2), 2) + Math.pow((zstv[i] - xv3), 2) );
		}
		
		//Determines if the satellite is above the surface of the earth
		double u = 0.0, l = 0.0;
		u = xv1*xstv+xv2*ystv+xv3*zstv;
		l = xv1*xv1+xv2*xv2+xv3*xv3;
		while(u > l) {
			System.out.println("Satellite 1 is above the surface");
		}
		
		System.out.println("This should be u1 of satelite 15: " + Data.u1[14]);
		System.out.println("This should be R: " + Data.r);
		System.out.println("X position of vehicle:" + xv1);
		System.out.println("Y position of vehicle:" + xv2);
		System.out.println("Z position of vehicle:" + xv3);
		System.out.println("X Position of first satellite: " + xstv);
		System.out.println("Y Position of first satellite: " + ystv);
		System.out.println("Z Position of first satellite: " + zstv);
		
		//Computes ts
		ts = tv-((xstv-xv1)*(xstv-xv1)+(ystv-xv2)*(ystv-xv2)+(zstv-xv3)*(zstv-xv3))/Data.c;
		System.out.println("ts: " + ts);
		
		
	}
	

	
	public static double getPsi(double ns, double ad, double am, double as) {
		return 2 * Data.pi * ns * ( (ad/360) + (am/21600) + (as/1296000) );
	}
	
	public static double getLambda(double ew, double bd, double bm, double bs) {
		return 2 * Data.pi * ew * ( (bd/360) + (bm/21600) + (bs/1296000) );
	}
}
