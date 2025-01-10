package NFlowShop;


public class MakeDelayMatrix{
    protected int[][] matrix;
    protected int[][] delayedMatrix;
    MakeDelayMatrix(int[][] matrix){
        this.matrix=matrix;
        int Size=matrix.length;
        int mcs=matrix[0].length;
        delayedMatrix=new int[Size][Size];

        for(int i=0;i<Size;i++){
            int delay=0;
            for(int j=0;j<Size;j++){
                if(i==j){
                    delayedMatrix[i][j]=0;
                }
                else{
                    int temp=0;
                    int temp1=0;
                    delay=matrix[i][0];
                    for(int x=1;x<mcs;x++){
                        int a=0;
                        int b=0;

                        for(int y=1;y<=x;y++){
                            a+=matrix[i][y];
                        }
                        for(int z=0;z<=x-1;z++){
                            b+=matrix[j][z];
                        }
                        temp=a-b;
                        temp1=Math.max(temp,temp1);
                    }
                    delay+=Math.max(temp1,0);
                    delayedMatrix[i][j]=delay;
                }
            }
        }

    }
    public int[][] returnDelayMatrix(){
        return delayedMatrix;
    }
//    public  static int[][] getMakeSpan(int[][] Matrix,int mcs,int x,int y) {
//        int[][] SpanMatrix = new int[2][mcs];
//        for (int i = 0; i < 2; i++) {
//            if (i == 0) {
//                int a = 0;
//                for (int j = 0; j < mcs; j++) {
//                    SpanMatrix[i][j] = a + Matrix[x][j];
//                    a = SpanMatrix[i][j];
//                }
//            } else {
//                for (int j = 0; j < mcs; j++) {
//                    if (j == 0) {
//                        SpanMatrix[i][j] = SpanMatrix[i - 1][j] + Matrix[y][j];
//                    } else {
//                        if (SpanMatrix[i][j - 1] > SpanMatrix[i - 1][j]) {
//                            SpanMatrix[i][j] = SpanMatrix[i][j - 1] + Matrix[y][j];
//                        } else {
//                            SpanMatrix[i][j] = SpanMatrix[i - 1][j] + Matrix[y][j];
//                        }
//                    }
//                }
//            }
//
//
//        }
//        return SpanMatrix;
//    }
}


