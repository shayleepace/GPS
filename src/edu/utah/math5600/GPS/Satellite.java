package edu.utah.math5600.GPS;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Scanner;

import Jama.Matrix;


public class Satellite {

	private static double[][] r3;
	private static double[][] x;
	private static double[][] xv;
	
	
	public static void main(String[] args) {
		//Get Data from data.dat
		Data.getData();		
		
		File dat = new File("satellite.dat");
		File log = new File("satellite.log");
		
		//create dat file if it doesn't exist
		if(!dat.exists()) {
			try {
				dat.createNewFile();
			} catch(IOException event) {
				event.printStackTrace();
			}
		}
		//create log file if it doesn't exist
		if(!log.exists()) {
			try {
				log.createNewFile();
			} catch(IOException event) {
				event.printStackTrace();
			}
		}
		
		//Matrices
		r3 = new double[3][3];
		x = new double[3][1];
		xv = new double[3][1];
		
		//Variables for vehicle input data
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
			if(vehicleInput.hasNextLine()) {
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
		double[] xst0 = new double[Data.altitude.length];
		double[] yst0 = new double[Data.altitude.length];
		double[] zst0 = new double[Data.altitude.length];
		double[] t0 = new double[Data.altitude.length];
		double[] t1 = new double[Data.altitude.length];
		double[] dft0 = new double[Data.altitude.length];
		double[] ft0 = new double[Data.altitude.length];
		double[] xst1 = new double[Data.altitude.length];
		double[] yst1 = new double[Data.altitude.length];
		double[] zst1 = new double[Data.altitude.length];
		double[] dft1 = new double[Data.altitude.length];
		double[] ft1 = new double[Data.altitude.length];
		double[] t2 = new double[Data.altitude.length];
		double[] xsts = new double[Data.altitude.length];
		double[] ysts = new double[Data.altitude.length];
		double[] zsts = new double[Data.altitude.length];
		
		//loop through satellites
		for(int i = 0; i < Data.altitude.length; i++) {
			//Computes position of each satellite at time tv
			xstv[i] = (Data.r+Data.altitude[i])*(Data.u1[i]*Math.cos(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i])+Data.v1[i]*Math.sin(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i]));
			ystv[i] = (Data.r+Data.altitude[i])*(Data.u2[i]*Math.cos(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i])+Data.v2[i]*Math.sin(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i]));
			zstv[i] = (Data.r+Data.altitude[i])*(Data.u3[i]*Math.cos(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i])+Data.v3[i]*Math.sin(2*Data.pi*tv/Data.periodicity[i]+Data.phase[i]));
		
			
			//Find t_0 for each satellite (start of Newton's Method)
			t0[i] = tv - (Math.sqrt( Math.pow((xstv[i] - xv1), 2) + Math.pow((ystv[i] - xv2), 2) + Math.pow((zstv[i] - xv3), 2))/Data.c) ;
			
			//Computes position of each satellite at time t_0
			xst0[i] = (Data.r+Data.altitude[i])*(Data.u1[i]*Math.cos(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i])+Data.v1[i]*Math.sin(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i]));
			yst0[i] = (Data.r+Data.altitude[i])*(Data.u2[i]*Math.cos(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i])+Data.v2[i]*Math.sin(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i]));
			zst0[i] = (Data.r+Data.altitude[i])*(Data.u3[i]*Math.cos(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i])+Data.v3[i]*Math.sin(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i]));
		
			
			//Find f(t_0) and df(t_0) for each satellite 
			double[][] a = new double[1][3];
			double[][] b = new double[3][1];
			a[0][0] = xst0[i]-xv1;
			a[0][1] = yst0[i]-xv2;
			a[0][2] = zst0[i]-xv3;
			b[0][0] = -Data.u1[i]*Math.sin(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i])+Data.v1[i]*Math.cos(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i]);
			b[1][0] = -Data.u2[i]*Math.sin(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i])+Data.v2[i]*Math.cos(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i]);
			b[2][0] = -Data.u3[i]*Math.sin(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i])+Data.v3[i]*Math.cos(2*Data.pi*t0[i]/Data.periodicity[i]+Data.phase[i]);
			
			Matrix A = new Matrix(a);
			Matrix B = new Matrix(b);
			Matrix C = A.times(B);
			Matrix D = A.transpose().times(A);
			
			//get dft0 and ft0 for each satellite
			dft0[i] = (4*Data.pi*(Data.r+h)/Data.periodicity[i])*(C.getArrayCopy()[0][0])+2*Data.c*Data.c*(tv-t0[i]); 
 			ft0[i] = (D.getArrayCopy()[0][0])-(Data.c*Data.c*(tv-t0[i])*(tv-t0[i]));
 			
 			//Find t_1 for each satellite (First step of Newton's method)
			t1[i] = t0[i] - ft0[i]/dft0[i];
			
			//Computes position of each satellite at time t_1
			xst1[i] = (Data.r+Data.altitude[i])*(Data.u1[i]*Math.cos(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i])+Data.v1[i]*Math.sin(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i]));
			yst1[i] = (Data.r+Data.altitude[i])*(Data.u2[i]*Math.cos(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i])+Data.v2[i]*Math.sin(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i]));
			zst1[i] = (Data.r+Data.altitude[i])*(Data.u3[i]*Math.cos(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i])+Data.v3[i]*Math.sin(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i]));
			
			
			//Find f(t_1) and df(t_1) for each satellite 
			double[][] e = new double[1][3];
			double[][] f = new double[3][1];
			e[0][0] = xst0[i]-xv1;
			e[0][1] = yst0[i]-xv2;
			e[0][2] = zst0[i]-xv3;
			f[0][0] = -Data.u1[i]*Math.sin(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i])+Data.v1[i]*Math.cos(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i]);
			f[1][0] = -Data.u2[i]*Math.sin(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i])+Data.v2[i]*Math.cos(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i]);
			f[2][0] = -Data.u3[i]*Math.sin(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i])+Data.v3[i]*Math.cos(2*Data.pi*t1[i]/Data.periodicity[i]+Data.phase[i]);
			
			Matrix E = new Matrix(e);
			Matrix F = new Matrix(f);
			Matrix G = E.times(F);
			Matrix H = E.transpose().times(E);
			
			//Get dft1 and ft1 for each satellite
			dft1[i] = (4*Data.pi*(Data.r+h)/Data.periodicity[i])*(G.getArrayCopy()[0][0])+2*Data.c*Data.c*(tv-t1[i]);
 			ft1[i] = (H.getArrayCopy()[0][0])-Data.c*Data.c*(tv-t1[i])*(tv-t1[i]);
 			
 			//Find t_2 for each satellite (First step of Newton's method)-- this is ts and needs to be printed out
			t2[i] = t1[i] - ft1[i]/dft1[i];
		
			//Computes position of each satellite at time t_2 (This is xs(ts) and needs to be printed out)
			xsts[i] = (Data.r+Data.altitude[i])*(Data.u1[i]*Math.cos(2*Data.pi*t2[i]/Data.periodicity[i]+Data.phase[i])+Data.v1[i]*Math.sin(2*Data.pi*t2[i]/Data.periodicity[i]+Data.phase[i]));
			ysts[i] = (Data.r+Data.altitude[i])*(Data.u2[i]*Math.cos(2*Data.pi*t2[i]/Data.periodicity[i]+Data.phase[i])+Data.v2[i]*Math.sin(2*Data.pi*t2[i]/Data.periodicity[i]+Data.phase[i]));
			zsts[i] = (Data.r+Data.altitude[i])*(Data.u3[i]*Math.cos(2*Data.pi*t2[i]/Data.periodicity[i]+Data.phase[i])+Data.v3[i]*Math.sin(2*Data.pi*t2[i]/Data.periodicity[i]+Data.phase[i]));
		
			//Prints out results for the satellites	
			System.out.println("ts for satellite " + (i+1) + " is: " + t2[i]);
			System.out.println("xsts for satellite " + (i+1) + " is: " + xsts[i]);
			System.out.println("ysts for satellite " + (i+1) + " is: " + ysts[i]);
			System.out.println("zsts for satellite " + (i+1) + " is: " + zsts[i]);
			
			
			
			
			
			try {
				PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(dat, true)));
				//Output- ts[i], xsts, ysts, zsts, xv1, xv2, xv3 to satellite.dat (receiver)
				if(i == 0) {
					out.println(xv1);
					out.println(xv2);
					out.println(xv3);
				}
				
				out.println(t2[i] + " " + xsts[i] + " " + ysts[i] + " " + zsts[i]);
				out.close();
			} catch (FileNotFoundException event) {
				event.printStackTrace();
			}
			
			try {
				PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(log, true)));
				//Output to satellite.log: All inputs,
				if(i == 0){
					out.println("Position of vehicle = " + "(" + xv1 + ", " + xv2 + ", " + xv3 + ")");
					out.println("pi = " + Data.pi);
					out.println("c = " + Data.c);
					out.println("r = " + Data.r);
					out.println("s = " + Data.s);
				}
				out.println("u1[" + i + "] = " + Data.u1[i]);
				out.println("u2[" + i + "] = " + Data.u2[i]);
				out.println("u3[" + i + "] = " + Data.u3[i]);
				out.println("v1[" + i + "] = " + Data.v1[i]);
				out.println("v2[" + i + "] = " + Data.v2[i]);
				out.println("v3[" + i + "] = " + Data.v3[i]);
				out.println("Periodicity of satellite[" + i + "] = " + Data.periodicity[i]);
				out.println("Altitude of satellite[" + i + "] = " + Data.altitude[i]);
				out.println("Phase of satellite[" + i + "] = " + Data.phase[i]);
				out.println("Time vehicle receives signal =" + tv);
				out.println("Latitude =" + ad + "\u00B0" + am + "'" + as + "\"");
				out.println("NS =" + ns);
				out.println("Longitude =" + bd + "\u00B0" + bm + "'" + bs + "\"");
				out.println("EW =" + ew);
				out.println("h of vehicle = " +  h);
				out.println("Time satellite" + " " + i + " " + "sends a signal" + "=" + t2[i]);
				out.println("Position of satellite" + " " + i +  " " + "=" + "(" + xsts[i] + "," + ysts[i] + "," + zsts[i] + ")");
				out.close();
				
			} catch (FileNotFoundException event) {
				event.printStackTrace();
			}
			
			
			
			
				
		
		} //End of For loop
		
//		
//		

	}
	

	
	public static double getPsi(double ns, double ad, double am, double as) {
		return 2 * Data.pi * ns * ( (ad/360) + (am/21600) + (as/1296000) );
	}
	
	public static double getLambda(double ew, double bd, double bm, double bs) {
		return 2 * Data.pi * ew * ( (bd/360) + (bm/21600) + (bs/1296000) );
	}
}
