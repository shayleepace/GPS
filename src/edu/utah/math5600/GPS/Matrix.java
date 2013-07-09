package edu.utah.math5600.GPS;

final public class Matrix {
    private final int row;
    private final int col;
    private final double[][] data;

    // Matrix of all zero's
    public Matrix(int row, int col) {
        this.row = row;
        this.col = col;
        data = new double[row][col];
    }

    // Matrix from 2d array
    public Matrix(double[][] data) {
    	row = data.length;
        col = data[0].length;    
        this.data = data.clone();
    }
    
    // Multiply two matrices together
    public Matrix multiply(Matrix B) {
        Matrix A = this;
        if (A.col != B.row) throw new RuntimeException("Illegal rowatrix dirowecolsiocols.");
        Matrix C = new Matrix(A.row, B.col);
        for (int i = 0; i < C.row; i++)
            for (int j = 0; j < C.col; j++)
                for (int k = 0; k < A.col; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }
    
    //Print out this matrix
    public void print() {
    	for (int i = 0; i < data.length; i++) {
	    	for (int j = 0; j < data[1].length; j++) {
	    		System.out.print(data[i][j] + " ");
	    	}
	    	System.out.print("\n");
	    }
    }
    
    //Print out 2d array matrix
    public static void print(double[][] matrix) {
	    for (int i = 0; i < matrix.length; i++) {
	    	for (int j = 0; j < matrix[1].length; j++) {
	    		System.out.print(matrix[i][j] + " ");
	    	}
	    	System.out.print("\n");
	    }
	}
    
    //Main for testing
    public static void main(String[] args) {
		//testing matrices
		double[][] a = { {1.0,2.0,3.0}, 
						 {4.0,5.0,6.0},
						 {1.0,1.0,1.0} };
		double[][] b = { {1.0},
						 {1.0},
						 {2.0} };
		Matrix A = new Matrix(a);
		Matrix B = new Matrix(b);
		Matrix C = A.multiply(B);
		A.print();
		System.out.println("*");
		B.print();
		System.out.println("=");
		C.print();
    }
}
