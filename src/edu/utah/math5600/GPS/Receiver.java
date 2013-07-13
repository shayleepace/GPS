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

import Jama.LUDecomposition;
import Jama.Matrix;

public class Receiver {

	public static void main(String[] args) {
		//Get Data from data.dat
		Data.getData();
		
		File log = new File("receiver.log");
		
		//create log file if it doesn't exist
		if(!log.exists()) {
			try {
				log.createNewFile();
			} catch(IOException event) {
				event.printStackTrace();
			}
		}
		
		//Variables for satellite input data
		double xv1 = 0.0, xv2 = 0.0, xv3 = 0.0;
		double[] t2 = new double[Data.altitude.length];
		double[] xsts = new double[Data.altitude.length];
		double[] ysts = new double[Data.altitude.length];
		double[] zsts = new double[Data.altitude.length];
		
		int[] j = new int[Data.altitude.length];
		int[] k = new int[24];
		int count = 0;
		
		//Get satellite input
		BufferedReader satelliteInput = new BufferedReader(new InputStreamReader(System.in));
		try {
			while(!satelliteInput.readLine().isEmpty()) {
				String data[] = satelliteInput.readLine().trim().split(" ");
				
				k[count] = Integer.parseInt(data[0]);
				t2[count] = Double.parseDouble(data[1]);
				xsts[count] = Double.parseDouble(data[2]);
				ysts[count] = Double.parseDouble(data[3]);
				zsts[count] = Double.parseDouble(data[4]);
				count++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		for(int i = 0; i < Data.altitude.length; i++) {			
			//Determines if the satellite is above the surface of the earth
			double u = 0.0, l = 0.0;
			u = xv1*xsts[i]+xv2*ysts[i]+xv3*zsts[i];
			l = xv1*xv1+xv2*xv2+xv3*xv3;
			count = 0;
			if(u > l) {
				//good satellite
				j[count] = i;
				count++;
			}
		}
		
		double x0 = 2.6567444497471053E7;
		double y0 = -210.25728967395753;
		double z0 = -300.27852915086174;
		double H, J, K, L;

			
		//Jacobian of x
		//Matrix A=F(x) Matrix B=Jacobian
		double[] a = new double[3];
		double[][] b = new double[3][3];
		double[] x = new double[3];
		Matrix A, B;
		
	for(int i = 0; i < 4; i++) {	
		H = (Math.sqrt( (xsts[j[0]]-x0)*(xsts[j[0]]-x0) + (ysts[j[0]]-y0)*(ysts[j[0]]-y0) + (zsts[j[0]]-z0)*(zsts[j[0]]-z0)));
		J = (Math.sqrt( (xsts[j[1]]-x0)*(xsts[j[1]]-x0) + (ysts[j[1]]-y0)*(ysts[j[1]]-y0) + (zsts[j[1]]-z0)*(zsts[j[1]]-z0)));
		K = (Math.sqrt( (xsts[j[2]]-x0)*(xsts[j[2]]-x0) + (ysts[j[2]]-y0)*(ysts[j[2]]-y0) + (zsts[j[2]]-z0)*(zsts[j[2]]-z0)));
		L = (Math.sqrt( (xsts[j[3]]-x0)*(xsts[j[3]]-x0) + (ysts[j[3]]-y0)*(ysts[j[3]]-y0) + (zsts[j[3]]-z0)*(zsts[j[3]]-z0)));

		a[0] =  (Math.sqrt( (xsts[j[1]]-x0)*(xsts[j[1]]-x0) + (ysts[j[1]]-y0)*(ysts[j[1]]-y0) + (zsts[j[1]]-z0)*(zsts[j[1]]-z0))-Math.sqrt( (xsts[j[0]]-x0)*(xsts[j[0]]-x0) + (ysts[j[0]]-y0)*(ysts[j[0]]-y0) + (zsts[j[0]]-z0)*(zsts[j[0]]-z0))-Data.c*(t2[j[0]]-t2[j[1]]));
		a[1] =  (Math.sqrt( (xsts[j[2]]-x0)*(xsts[j[2]]-x0) + (ysts[j[2]]-y0)*(ysts[j[2]]-y0) + (zsts[j[2]]-z0)*(zsts[j[2]]-z0))-Math.sqrt( (xsts[j[1]]-x0)*(xsts[j[1]]-x0) + (ysts[j[1]]-y0)*(ysts[j[1]]-y0) + (zsts[j[1]]-z0)*(zsts[j[1]]-z0))-Data.c*(t2[j[1]]-t2[j[2]]));
		a[2] =  (Math.sqrt( (xsts[j[3]]-x0)*(xsts[j[3]]-x0) + (ysts[j[3]]-y0)*(ysts[j[3]]-y0) + (zsts[j[3]]-z0)*(zsts[j[3]]-z0))-Math.sqrt( (xsts[j[2]]-x0)*(xsts[j[2]]-x0) + (ysts[j[2]]-y0)*(ysts[j[2]]-y0) + (zsts[j[2]]-z0)*(zsts[j[2]]-z0))-Data.c*(t2[j[2]]-t2[j[3]]));
		b[0][0] =	(xsts[j[0]]-x0)/H-(xsts[j[1]]-x0)/J;
		b[0][1] = 	(ysts[j[0]]-y0)/H-(ysts[j[1]]-y0)/J;
		b[0][2] = 	(zsts[j[0]]-z0)/H-(zsts[j[1]]-z0)/J;
		b[1][0] = 	(xsts[j[1]]-x0)/J-(xsts[j[2]]-x0)/K;
		b[1][1] = 	(ysts[j[1]]-y0)/J-(ysts[j[2]]-y0)/K;
		b[1][2] = 	(zsts[j[1]]-y0)/J-(zsts[j[2]]-y0)/K;
		b[2][0] = 	(xsts[j[2]]-x0)/K-(xsts[j[3]]-x0)/L;
		b[2][1] = 	(ysts[j[2]]-y0)/K-(ysts[j[3]]-y0)/L;
		b[2][2] = 	(zsts[j[2]]-y0)/K-(zsts[j[3]]-y0)/L;
		
//		A = new Matrix(a);
//		B = new Matrix(b);
		//X = B.solve(A); //3x1
		//LUDecomposition lu = new LUDecomposition(B);
		x = lsolve(b, a);

		x0 = x[0];
		y0 = x[1];
		z0 = x[2];
		
//		Print out vehicle position
		//System.out.println(ad + " " + am + " " + as + " " + ns + " " + bd + " " + bm + " " + bs + " " + ew + " " + h);
		
		//Log files to receiver.log
		try {
			PrintStream out = new PrintStream(new BufferedOutputStream(new FileOutputStream(log, true)));
			//Output to receiver.log: All inputs,
			out.println("Receiver log, Shaylee Pace and Mitch Norton data.dat:");
			out.println("pi = " + Data.pi);
			out.println("c = " + Data.c);
			out.println("r = " + Data.r);
			out.println("s = " + Data.s);
			out.println(" end data.dat\n");
			
			//out.println(ad + " " + am + " " + as + " " + ns + " " + bd + " " + bm + " " + bs + " " + ew + " " + h);
			
			out.close();
			
		} catch (FileNotFoundException event) {
			event.printStackTrace();
		}
	}
	
		
		
		
		
		//Call new x= (xk, yk, zk), lat, long
//		h = Math.sqrt(xk*xk+yk*yk+zk*zk);
//		Lat = Math.atan(zk/(Math.sqrt(xk*xk+yk*yk)));
//		Long = Math.atan(yk/xk); //if x>0 and y>0
//		Long = Data.pi + Math.atan(yk/xk); // if x<0
//		Long = 2*Data.pi + Math.atan(yk/xk); // if x>0 y<0
		

	}
	
	 private static final double EPSILON = 1e-10;

	    // Gaussian elimination with partial pivoting
	    public static double[] lsolve(double[][] A, double[] b) {
	        int N  = b.length;

	        for (int p = 0; p < N; p++) {

	            // find pivot row and swap
	            int max = p;
	            for (int i = p + 1; i < N; i++) {
	                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
	                    max = i;
	                }
	            }
	            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
	            double   t    = b[p]; b[p] = b[max]; b[max] = t;

	            // singular or nearly singular
	            if (Math.abs(A[p][p]) <= EPSILON) {
	                throw new RuntimeException("Matrix is singular or nearly singular");
	            }

	            // pivot within A and b
	            for (int i = p + 1; i < N; i++) {
	                double alpha = A[i][p] / A[p][p];
	                b[i] -= alpha * b[p];
	                for (int j = p; j < N; j++) {
	                    A[i][j] -= alpha * A[p][j];
	                }
	            }
	        }

	        // back substitution
	        double[] x = new double[N];
	        for (int i = N - 1; i >= 0; i--) {
	            double sum = 0.0;
	            for (int j = i + 1; j < N; j++) {
	                sum += A[i][j] * x[j];
	            }
	            x[i] = (b[i] - sum) / A[i][i];
	        }
	        return x;
	    }
}