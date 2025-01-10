package NFlowShop;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

public class Demo {
    static int jobs, mcs;
    static final int MAXCOLS = 510;
    static long timeseed, ts;

    // Random number generator based on the Linear Congruential Generator (LCG)
    static int unif(long[] seed, int low, int high) {
        long m = 2147483647L, a = 16807L, b = 127773L, c = 2836L;
        int k;
        double val01;
        k = (int) (seed[0] / b);
        seed[0] = a * (seed[0] % b) - k * c;
        if (seed[0] < 0) {
            seed[0] = seed[0] + m;
        }
        val01 = (double) seed[0] / (double) m;
        return (low + (int) (val01 * (high - low + 1)));
    }

    public static void FileRead() throws IOException {
        int[][] a = new int[jobs][mcs];
        BufferedReader br = new BufferedReader(new FileReader(System.getProperty("user.dir") + "\\NFlowShop\\tlrdata.txt"));

        String line;
        String instance = "";

        while ((line = br.readLine()) != null) {
            StringTokenizer stk = new StringTokenizer(line);
            instance = stk.nextToken();
            jobs = Integer.parseInt(stk.nextToken());
            mcs = Integer.parseInt(stk.nextToken());
            timeseed = Long.parseLong(stk.nextToken());

            // Read job processing times
            for (int i = 0; i < jobs; i++) {
                for (int j = 0; j < mcs; j++) {
                    if (stk.hasMoreTokens()) {
                        a[i][j] = Integer.parseInt(stk.nextToken());
                    }
                }
            }
        }
        br.close();

        // Print the read data for verification
        System.out.println("Instance: " + instance);
        System.out.println("Jobs: " + jobs);
        System.out.println("Machines: " + mcs);
        System.out.println("Time Seed: " + timeseed);

        for (int i = 0; i < jobs; i++) {
            for (int j = 0; j < mcs; j++) {
                System.out.print(a[i][j] + " ");
            }
            System.out.println();
        }
    }

    public static void main(String[] args) {
        try {
            FileRead();
            // Additional logic can be added here
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
