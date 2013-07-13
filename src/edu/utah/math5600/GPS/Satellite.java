package edu.utah.math5600.GPS;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import Jama.Matrix;


public class Satellite {

	private static double[][] r3;
	private static double[][] x;
	private static double[][] xv;
	
	
	public static void main(String[] args) {
		//Get Data from data.dat
		Data.getData();		
		
		File log = new File("satellite.log");
		
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
		BufferedReader vehicleInput = new BufferedReader(new InputStreamReader(System.in));
		try {
			String data[] = vehicleInput.readLine().trim().split(" ");
			
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
		
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// Testing vehicle input
		// tv = 102123.0;
		// ad = 40;
		// am =  45; 
		// as = 55.0;
		// ns =  1; 
		// bd = 111; 
		// bm = 50; 
		// bs = 58.0; 
		// ew = -1; 
		// h = 1372.0;
		
//		System.out.println("This should be h: " + h);
//		System.out.println("This should be bd: " + bd);
		
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
			// System.out.println("ts for satellite " + (i+1) + " is: " + t2[i]);
			// System.out.println("xsts for satellite " + (i+1) + " is: " + xsts[i]);
			// System.out.println("ysts for satellite " + (i+1) + " is: " + ysts[i]);
			// System.out.println("zsts for satellite " + (i+1) + " is: " + zsts[i]);
			
			
			//Determines if the satellite is above the surface of the earth
			double u = 0.0, l = 0.0;
			u = xv1*xsts[i]+xv2*ysts[i]+xv3*zsts[i];
			l = xv1*xv1+xv2*xv2+xv3*xv3;
			int count = 0;
			if(u > l) {
				//Output- is, ts, and xs to receiver
				System.out.println(i + " " + t2[i] + " " + xsts[i] + " " + ysts[i] + " " + zsts[i]);
			}
			
		
		} //End of For loop
		
		try {
			PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(log, true)));
			//Output to satellite.log: All inputs,
			out.println("Satellite log, Shaylee Pace and Mitch Norton data.dat:");
			out.println("pi = " + Data.pi);
			out.println("c = " + Data.c);
			out.println("r = " + Data.r);
			out.println("s = " + Data.s);
			
			for(int i = 0; i < Data.altitude.length; i++) {
				out.println("u1[" + i + "] = " + Data.u1[i]);
				out.println("u2[" + i + "] = " + Data.u2[i]);
				out.println("u3[" + i + "] = " + Data.u3[i]);
				out.println("v1[" + i + "] = " + Data.v1[i]);
				out.println("v2[" + i + "] = " + Data.v2[i]);
				out.println("v3[" + i + "] = " + Data.v3[i]);
				out.println("Periodicity of satellite[" + i + "] = " + Data.periodicity[i]);
				out.println("Altitude of satellite[" + i + "] = " + Data.altitude[i]);
				out.println("Phase of satellite[" + i + "] = " + Data.phase[i]);
			}
			
			out.println(" end data.dat\n");
			
			for(int i = 0; i < Data.altitude.length; i++) {
				out.println(i + " " + t2[i] + " " + xsts[i] + " " + ysts[i] + " " + zsts[i]);
			}
			out.close();
			
		} catch (FileNotFoundException event) {
			event.printStackTrace();
		}
	}

	public static double getPsi(double ns, double ad, double am, double as) {
		return 2 * Data.pi * ns * ( (ad/360) + (am/21600) + (as/1296000) );
	}
	
	public static double getLambda(double ew, double bd, double bm, double bs) {
		return 2 * Data.pi * ew * ( (bd/360) + (bm/21600) + (bs/1296000) );
	}
}

//get data from data.dat
class Data {
	public static double pi, c, r, s;
	public static double[] u1, u2, u3, v1, v2, v3, periodicity, altitude, phase;
	
	public static void getData() {
		u1 = new double[24]; u2 = new double[24]; u3 = new double[24];
		v1 = new double[24]; v2 = new double[24]; v3 = new double[24];
		periodicity = new double[24]; altitude = new double[24]; phase = new double[24];
		
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