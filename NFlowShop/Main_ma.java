package NFlowShop;
import NFlowShop.MakeDelayMatrix;
import NFlowShop.Scheduling;
public class Main_ma {
    public static void main(String[] args) {
        int jbs=5;
        int mcs=5;
        int[][] jbmatrix={{2,3,1,4,7},
                 {4,5,3,6,9},{4,7,2,3,2},{6,4,3,2,7},{6,5,3,6,4}
                  };
        int[] dueDates ={4,3,5,6,4};
        MakeDelayMatrix md=new MakeDelayMatrix(jbmatrix);
        int[][] delay = md.returnDelayMatrix();
        Scheduling sc=new Scheduling(jbs,mcs,jbmatrix,delay,dueDates );

}
}

