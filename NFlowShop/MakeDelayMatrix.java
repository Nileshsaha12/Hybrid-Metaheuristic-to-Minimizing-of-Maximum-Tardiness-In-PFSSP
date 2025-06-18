package  NFlowShop;
public class MakeDelayMatrix {
    protected int[][] matrix;
    protected int[][] delayedMatrix;

    MakeDelayMatrix(int[][] matrix) {
        this.matrix = matrix;
        int Size = matrix.length;
        int mcs = matrix[0].length;
        delayedMatrix = new int[Size][Size];

        for (int i = 0; i < Size; i++) {
            for (int j = 0; j < Size; j++) {
                if (i == j) {
                    delayedMatrix[i][j] = 0;
                } else {
                   // System.out.println("j"+i+"j"+j);
                    int delay = 0; // FIX: Initialize with p_{i,1}
                    int temp1 = 0;
                    for (int x = 1; x < mcs; x++) {
                        int a = 0, b = 0;
                        for (int y = 1; y <= x; y++) {
                            a += matrix[i][y]; // Sum p_i from machine 1 to x
                        }
                        for (int z = 0; z <= x - 1; z++) {
                            b += matrix[j][z]; // Sum p_j from machine 0 to x-1
                        }
                        int temp = a - b;
                        temp1 = Math.max(temp, temp1);
                    }
                    delay = Math.max(temp1, 0); // Add max difference
                    delayedMatrix[i][j] = delay;
                }
            }
        }
    }

    public int[][] returnDelayMatrix() {
        return delayedMatrix;
    }
}


