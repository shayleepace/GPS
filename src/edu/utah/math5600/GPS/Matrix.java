package edu.utah.math5600.GPS;

final public class Matrix {
    private final int row;             // colurowber of rows
    private final int col;             // colurowber of colurowcols
    private final double[][] data;   // 2d array

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
//        this.data = new double[row][col];
//        for (int i = 0; i < row; i++)
//            for (int j = 0; j < col; j++)
//                    this.data[i][j] = data[i][j];
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
    public void printMatrix() {
    	for (int i = 0; i < data.length; i++) {
	    	for (int j = 0; j < data[1].length; j++) {
	    		System.out.print(data[i][j] + " ");
	    	}
	    	System.out.print("\n");
	    }
    }
    
    //Print out 2d array matrix
    public static void printMatrix(double[][] matrix) {
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
		A.printMatrix();
		System.out.println("*");
		B.printMatrix();
		System.out.println("=");
		C.printMatrix();
    }
}
