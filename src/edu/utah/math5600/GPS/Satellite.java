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
		double xs1tv = 0.0, ys1tv = 0.0, zs1tv = 0.0, ts = 0.0;
		
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
		
		
		//Computes position of first satellite at time tv
		xs1tv = (Data.r+Data.altitude[0])*(Data.u1[0]*Math.cos(2*Data.pi*tv/Data.periodicity[0]+Data.phase[0])+Data.v1[0]*Math.sin(2*Data.pi*tv/Data.periodicity[0]+Data.phase[0]));
		ys1tv = (Data.r+Data.altitude[0])*(Data.u2[0]*Math.cos(2*Data.pi*tv/Data.periodicity[0]+Data.phase[0])+Data.v2[0]*Math.sin(2*Data.pi*tv/Data.periodicity[0]+Data.phase[0]));
		zs1tv = (Data.r+Data.altitude[0])*(Data.u3[0]*Math.cos(2*Data.pi*tv/Data.periodicity[0]+Data.phase[0])+Data.v3[0]*Math.sin(2*Data.pi*tv/Data.periodicity[0]+Data.phase[0]));
		
		//Determines if the satellite is above the surface of the earth
		double u = 0.0, l = 0.0;
		u = xv1*xs1tv+xv2*ys1tv+xv3*zs1tv;
		l = xv1*xv1+xv2*xv2+xv3*xv3;
		while(u > l) {
			System.out.println("Satellite 1 is above the surface");
		}
		
		System.out.println("This should be u1 of satelite 15: " + Data.u1[14]);
		System.out.println("This should be R: " + Data.r);
		System.out.println("X position of vehicle:" + xv1);
		System.out.println("Y position of vehicle:" + xv2);
		System.out.println("Z position of vehicle:" + xv3);
		System.out.println("X Position of first satellite: " + xs1tv);
		System.out.println("Y Position of first satellite: " + ys1tv);
		System.out.println("Z Position of first satellite: " + zs1tv);
		
		//Computes ts
		ts = tv-((xs1tv-xv1)*(xs1tv-xv1)+(ys1tv-xv2)*(ys1tv-xv2)+(zs1tv-xv3)*(zs1tv-xv3))/Data.c;
		System.out.println("ts: " + ts);
		
		
	}
	

	
	public static double getPsi(double ns, double ad, double am, double as) {
		return 2 * Data.pi * ns * ( (ad/360) + (am/21600) + (as/1296000) );
	}
	
	public static double getLambda(double ew, double bd, double bm, double bs) {
		return 2 * Data.pi * ew * ( (bd/360) + (bm/21600) + (bs/1296000) );
	}
}
