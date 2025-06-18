package NFlowShop;
import javax.swing.*;
import java.util.*;
import java.util.stream.Collectors;
public class Scheduling {

    // Problem parameters
    int numJobs;
    int numMachines;
    int[][] processingTimes;
    int[] dueDates;
    int[][] delayMatrix;

    // Algorithm parameters
    int populationSize;
    int maxGenerationsGA;
    int maxIterationsJAYA;
    int maxIterationsSA;
    double crossoverRate;
    double mutationRate;
    double initialTemperature;
    double coolingRate;
    List<Integer> bestOverallSchedule;
    double bestOverallTardiness = Double.MAX_VALUE;
    int GApopulation;

    public Scheduling(int jobs, int mcs, int[][] jbmatrix, int[][] delayMatrix, int[] dueDates,
                      int InitialpopulationSize, int GApopulation, int maxGenerationsGA, int maxIterationsJAYA, int maxIterationsSA,
                      double crossoverRate, double mutationRate, double initialTemperature, double coolingRate) {
        this.numJobs = jobs;
        this.numMachines = mcs;
        this.processingTimes = jbmatrix;
        this.delayMatrix = delayMatrix;
        this.dueDates = dueDates;
        this.populationSize = InitialpopulationSize;
        this.maxGenerationsGA = maxGenerationsGA;
        this.maxIterationsJAYA = maxIterationsJAYA;
        this.maxIterationsSA = maxIterationsSA;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
        this.initialTemperature = initialTemperature;
        this.coolingRate = coolingRate;
        this.GApopulation= GApopulation;

    }




    public List<Integer> getBestSchedule() {
        return new ArrayList<>(bestOverallSchedule);
    }
    public double getBestTardiness() {
        return bestOverallTardiness;
    }
    private void updateBest(List<Integer> schedule) {
        double tardiness = calculateMaxTardiness(schedule, processingTimes, dueDates, delayMatrix);
        if (tardiness < bestOverallTardiness) {
            bestOverallSchedule = new ArrayList<>(schedule);
            bestOverallTardiness = tardiness;
        }
    }

    public static int[][] calculateCompletionTimes(List<Integer> schedule, int[][] processingTimes, int[][] delayMatrix) {
        // Validate the schedule
        for (int job : schedule) {
            if (job < 0 || job >= processingTimes.length) {
                throw new IllegalArgumentException("Invalid job index in schedule: " + job);
            }
        }
        int jobs = schedule.size();
        int mcs = processingTimes[0].length;
        int[][] SpanMatrix = new int[jobs][mcs];

        for (int i = 0; i < jobs; i++) {
            int n = schedule.get(i);
            if (i == 0) {
                int a = 0;
                for (int j = 0; j < mcs; j++) {
                    SpanMatrix[i][j] = a + processingTimes[n][j];
                    a = SpanMatrix[i][j];
                }
            } else {
                int prevJob = schedule.get(i - 1);
                int delay = delayMatrix[prevJob][n]; // Add delay from previous job to current job
                for (int j = 0; j < mcs; j++) {
                    if (j == 0) {
                        SpanMatrix[i][j] = SpanMatrix[i - 1][j] + processingTimes[n][j] + delay;
                    } else {
                        SpanMatrix[i][j] = SpanMatrix[i][j - 1] + processingTimes[n][j];
                    }
                }
            }
        }

        return SpanMatrix;
    }




    public static double calculateMaxTardiness(List<Integer> schedule, int[][] processingTimes, int[] dueDates, int[][] delayMatrix) {
        // Calculate completion times
        int[][] completionTimes = calculateCompletionTimes(schedule, processingTimes, delayMatrix);
        int numJobs = schedule.size();
        double maxTardiness = 0;

        for (int j = 0; j < numJobs; j++) {
            int job = schedule.get(j); // Get the job ID from the schedule
            // Access the completion time of the last machine for this job
            int lastMachineIndex = processingTimes[0].length - 1;
            double tardiness = Math.max(0, completionTimes[j][lastMachineIndex] - dueDates[job]);
            if (tardiness > maxTardiness) {
                maxTardiness = tardiness;
            }
        }

        return maxTardiness;
    }
    private List<List<Integer>> selectBest(List<List<Integer>> population, int count) {
        // Sort by tardiness (lower is better)
        population.sort(Comparator.comparingDouble(indiv ->
                calculateMaxTardiness(indiv, processingTimes, dueDates, delayMatrix)
        ));
        // Remove duplicate sequences
        List<List<Integer>> unique = population.stream()
                .distinct()
                .collect(Collectors.toList());
        // Return the top 'count' unique sequences
        return new ArrayList<>(unique.subList(0, Math.min(count, unique.size())));
    }



    public void optimize() {
        bestOverallTardiness = Double.MAX_VALUE;
        List<Integer> jayaBest = runJAYA();
        updateBest(jayaBest);
        System.out.println(jayaBest.toString());
        System.out.println(calculateMaxTardiness(jayaBest,processingTimes,dueDates,delayMatrix));
        List<Integer> gaBest = runGA(jayaBest);
        updateBest(gaBest);
        System.out.println(gaBest.toString());
        System.out.println(calculateMaxTardiness(gaBest,processingTimes,dueDates,delayMatrix));
        List<Integer> saBest = runSA(gaBest);
        updateBest(saBest);
        System.out.println(saBest.toString());
        System.out.println(calculateMaxTardiness(saBest,processingTimes,dueDates,delayMatrix));
    }




    public List<Integer> runJAYA() {
        // 1) Initialize population
        List<List<Integer>> population =  generateInitialPopulation(populationSize, numJobs,dueDates);
        // Sort population and remove duplicates
        population = selectBest(population, population.size());
//        for(List<Integer> pop:population){
//            System.out.print(pop +" :");
//            System.out.println(calculateMaxTardiness(pop,processingTimes,dueDates,delayMatrix));
//        }

        // Identify best and worst in the sorted population
        List<Integer> globalBest = new ArrayList<>(population.get(0));
        List<Integer> globalWorst = new ArrayList<>(population.get(population.size() - 1));
        double bestTardiness = calculateMaxTardiness(globalBest, processingTimes, dueDates, delayMatrix);
        double worstTardiness = calculateMaxTardiness(globalWorst, processingTimes, dueDates, delayMatrix);


        // 2) JAYA iterations
        for (int iteration = 0; iteration < maxIterationsJAYA; iteration++) {

            // Modify each solution based on best and worst
            for (int p = 0; p < population.size(); p++) {

                List<Integer> individual = population.get(p);
                List<Integer> newSolution = new ArrayList<>(individual);
                Random rand = new Random();

                // Apply JAYA update formula
                for (int i = 0; i < numJobs; i++) {
                    double r1 = rand.nextDouble();
                    double r2 = rand.nextDouble();
                    int current = individual.get(i);
                    int bestVal = globalBest.get(i);
                    int worstVal = globalWorst.get(i);
                    // JAYA update: move toward bestVal, away from worstVal
                    int newValue = (int) (current
                            + r1 * (bestVal - current)
                            - r2 * (worstVal - current));
                    // Ensure valid job index
                    newSolution.set(i, Math.abs(newValue) % numJobs);
                }

                // Repair the new solution
                newSolution = repairSolution(newSolution);

                // Check if the new solution is better
                double currentTardiness = calculateMaxTardiness(individual, processingTimes, dueDates, delayMatrix);




                double newTardiness = calculateMaxTardiness(newSolution, processingTimes, dueDates, delayMatrix);

                if (newTardiness < currentTardiness) {
                    population.set(p, new ArrayList<>(newSolution));
                }
            }

            // After updating the entire population, re-sort and pick best/worst again
            population = selectBest(population, population.size());
            globalBest = population.get(0);
            globalWorst = population.get(population.size() - 1);

            bestTardiness = calculateMaxTardiness(globalBest, processingTimes, dueDates, delayMatrix);
            worstTardiness = calculateMaxTardiness(globalWorst, processingTimes, dueDates, delayMatrix);
        }
        return globalBest;
    }

    public static List<List<Integer>> generateInitialPopulation(int populationSize, int numJobs, int[] dueDates) {
        List<List<Integer>> population = new ArrayList<>();

        // Create a list of job indices [0, 1, 2, ..., numJobs-1]
        List<Integer> sortedJobs = new ArrayList<>();
        for (int i = 0; i < numJobs; i++) {
            sortedJobs.add(i);
        }

        // Sort the jobs based on dueDates (lower due date comes first)
        sortedJobs.sort(Comparator.comparingInt(i -> dueDates[i]));

        // Use a Set to ensure unique sequences
        Set<List<Integer>> uniqueIndividuals = new HashSet<>();
        Random rand = new Random();

        // Keep generating individuals until we have the desired population size
        while (uniqueIndividuals.size() < populationSize ) {
            // Start with the sorted baseline
            List<Integer> individual = new ArrayList<>(sortedJobs);

            // Perform a random number of swaps (at least one) to randomize the order
            int numSwaps = rand.nextInt(numJobs) + 1;
            for (int s = 0; s < numSwaps; s++) {
                int idx1 = rand.nextInt(numJobs);
                int idx2 = rand.nextInt(numJobs);
                Collections.swap(individual, idx1, idx2);
            }






            // Add to the set if not already present
            uniqueIndividuals.add(individual);
            // System.out.println(individual);
        }

        // Convert the set to a list
        population.addAll(uniqueIndividuals);
        return population;
    }

    // Repair the solution to ensure it represents a valid permutation (no duplicates)
    private List<Integer> repairSolution(List<Integer> solution) {
        List<Integer> repaired = new ArrayList<>(solution);
        Set<Integer> seen = new HashSet<>();
        List<Integer> duplicateIndices = new ArrayList<>();

        // Identify duplicates
        for (int i = 0; i < repaired.size(); i++) {
            int job = repaired.get(i);
            if (!seen.add(job)) {
                duplicateIndices.add(i);
            }
        }

        // Find missing jobs
        List<Integer> missing = new ArrayList<>();
        for (int job = 0; job < numJobs; job++) {
            if (!seen.contains(job)) {
                missing.add(job);
            }
        }

        // Replace duplicates with missing jobs
        for (int i = 0; i < duplicateIndices.size() && i < missing.size(); i++) {
            repaired.set(duplicateIndices.get(i), missing.get(i));
        }
        return repaired;
    }

    private List<Integer> runGA(List<Integer> initialSolution) {

        // Generate initial population.
        List<List<Integer>> population =generateGAPopulation(initialSolution, GApopulation,dueDates);
        //System.out.println(population.toString());

        for (int generation = 0; generation < maxGenerationsGA; generation++) {
            // Phase 1: Select the best 'bestCount' individuals.
            List<List<Integer>> bestIndividuals = selectBest(population, GApopulation);





            // Phase 2: Apply Order Crossover.
            List<List<Integer>> offspringOrder = new ArrayList<>();
            for (int i = 0; i < bestIndividuals.size() / 2; i++) {
                // Pairing: using indices 2*i and 2*i+1
                List<Integer> parent1 = bestIndividuals.get(2 * i);
                List<Integer> parent2 = bestIndividuals.get(2*i +1);
                List<Integer> child1=orderCrossover(parent1, parent2);
                List<Integer> child2=orderCrossover(parent2, parent1);

                offspringOrder.add(parent1);
                offspringOrder.add(parent2);
                offspringOrder.add(child1);
                offspringOrder.add(child2);
            }
            // If there's an odd individual, pair it with a random best partner.
            if (bestIndividuals.size() % 2 == 1) {
                List<Integer> parent1 = bestIndividuals.get(bestIndividuals.size() - 1);
                List<Integer> parent2 = bestIndividuals.get(new Random().nextInt(bestIndividuals.size() - 1));
                offspringOrder.add(orderCrossover(parent1, parent2));
            }

            bestIndividuals = selectBest(offspringOrder, GApopulation);


            // Phase 3: Apply PMX Crossover.
            List<List<Integer>> offspringPMX = new ArrayList<>();
            for (int i = 0; i < bestIndividuals.size() / 2; i++) {
                List<Integer> parent1 = bestIndividuals.get( 2* i);
                List<Integer> parent2 = bestIndividuals.get(2*i + 1);
                List<Integer> child1 = pmxCrossover(parent1, parent2);
                List<Integer> child2 = pmxCrossover(parent2, parent1);

                offspringPMX.add(parent1);
                offspringPMX.add(parent2);
                offspringPMX.add(child1);
                offspringPMX.add(child2);
            }
            if (bestIndividuals.size() % 2 == 1) {
                List<Integer> parent1 = bestIndividuals.get(bestIndividuals.size() - 1);
                List<Integer> parent2 = bestIndividuals.get(new Random().nextInt(bestIndividuals.size() - 1));
                offspringPMX.add(pmxCrossover(parent1, parent2));
            }

            bestIndividuals = selectBest(offspringPMX, GApopulation);

            // Phase 4: Apply Inverse Mutation.
            List<List<Integer>> offspringInverse = new ArrayList<>();
            for (List<Integer> indiv : bestIndividuals) {
                List<Integer> child =inverseMutation(new ArrayList<>(indiv));




                offspringInverse.add(indiv);
                offspringInverse.add(child);
            }

            bestIndividuals = selectBest(offspringInverse, GApopulation);

            // Phase 5: Apply Pairwise Mutation.
            List<List<Integer>> offspringPairwise = new ArrayList<>();
            for (List<Integer> indiv : bestIndividuals) {
                List<Integer> child = pairwiseMutation(new ArrayList<>(indiv));

                offspringPairwise.add(indiv);
                offspringPairwise.add(child);
            }

            bestIndividuals = selectBest(offspringPairwise, GApopulation);

            // Update population for the next generation.
            population = bestIndividuals;
        }
        // Return the best schedule after all generations.
        return population.get(0);
    }

    private List<List<Integer>> generateGAPopulation(List<Integer> initialSolution, int populationSize, int[] dueDates) {

        Set<List<Integer>> uniqueIndividuals = new HashSet<>();
        Random rand = new Random();
        int numJobs = initialSolution.size();

        // Always include the original
        uniqueIndividuals.add(new ArrayList<>(initialSolution));

        // Keep going until we have enough distinct individuals
        while (uniqueIndividuals.size() < populationSize) {
            List<Integer> individual;

            do {
                // Start from the baseline
                individual = new ArrayList<>(initialSolution);

                // Pick distinct random indices
                int picks = Math.min(5, numJobs);
                Set<Integer> idxSet = new HashSet<>();
                while (idxSet.size() < picks) {
                    idxSet.add(rand.nextInt(numJobs));
                }
                List<Integer> idxList = new ArrayList<>(idxSet);




                // Extract jobs at these indices
                List<Integer> jobsToSort = new ArrayList<>();
                for (int idx : idxList) {
                    jobsToSort.add(individual.get(idx));
                }

                // Check if the extracted jobs are already sorted by due dates
                boolean isSorted = true;
                for (int i = 0; i < jobsToSort.size() - 1; i++) {
                    if (dueDates[jobsToSort.get(i)] > dueDates[jobsToSort.get(i + 1)]) {
                        isSorted = false;
                        break;
                    }
                }

                if (isSorted) {
                    // Shuffle to ensure a change if already sorted
                    Collections.shuffle(jobsToSort);
                } else {
                    // Sort by due dates otherwise
                    jobsToSort.sort(Comparator.comparingInt(job -> dueDates[job]));
                }

                // Replace the jobs back into the individual
                Collections.sort(idxList); // Ensure indices are in order
                for (int i = 0; i < idxList.size(); i++) {
                    individual.set(idxList.get(i), jobsToSort.get(i));
                }

            } while (!uniqueIndividuals.add(individual));
        }

        return new ArrayList<>(uniqueIndividuals);
    }

// ----- Genetic Operators -----

    // Order Crossover (OX)
    public static List<Integer> orderCrossover(List<Integer> parent1, List<Integer> parent2) {
        int length = parent1.size();
        Random rand = new Random();
        int start = rand.nextInt(length);
        int end = start + rand.nextInt(length - start);

        List<Integer> child = new ArrayList<>(Collections.nCopies(length, -1));
        // Copy a segment from parent1.
        for (int i = start; i < end; i++) {
            child.set(i, parent1.get(i));
        }




        // Fill remaining positions with genes from parent2 in order.
        int index = end;
        for (int gene : parent2) {
            if (!child.contains(gene)) {
                if (index >= length) {
                    index = 0;
                }
                child.set(index, gene);
                index++;
            }
        }
        return child;
    }

    // Partially Mapped Crossover (PMX)
    public static List<Integer> pmxCrossover(List<Integer> parent1, List<Integer> parent2) {
        int length = parent1.size();
        Random rand = new Random();
        int start = rand.nextInt(length);
        int end = start + rand.nextInt(length - start);

        List<Integer> child = new ArrayList<>(Collections.nCopies(length, -1));

        // Copy mapping segment from parent1.
        for (int i = start; i < end; i++) {
            child.set(i, parent1.get(i));
        }
        // Map remaining genes from parent2.
        for (int i = start; i < end; i++) {
            int gene = parent2.get(i);
            if (!child.contains(gene)) {
                int pos = i;
                while (true) {
                    int mappedGene = parent1.get(pos);
                    pos = parent2.indexOf(mappedGene);
                    if (child.get(pos) == -1) {
                        child.set(pos, gene);
                        break;
                    }
                }
            }
        }
        // Fill in the remaining positions.
        for (int i = 0; i < length; i++) {
            if (child.get(i) == -1) {
                child.set(i, parent2.get(i));
            }
        }
        return child;
    }


    // Inverse Mutation: Reverse a random segment.
    public static List<Integer> inverseMutation(List<Integer> individual) {
        int length = individual.size();
        Random rand = new Random();
        int i = rand.nextInt(length);
        int j = i + rand.nextInt(length - i);
        Collections.reverse(individual.subList(i, j + 1));
        return individual;
    }
    // Pairwise Mutation: Swap two random elements.
    public static List<Integer> pairwiseMutation(List<Integer> individual) {
        int length = individual.size();
        Random rand = new Random();
        int i = rand.nextInt(length);
        int j = rand.nextInt(length);
        while (i == j) {
            j = rand.nextInt(length);
        }
        Collections.swap(individual, i, j);
        return individual;
    }
    private List<Integer> runSA(List<Integer> initialSolution) {
        List<Integer> currentSolution = new ArrayList<>(initialSolution);
        List<Integer> bestSolution = new ArrayList<>(currentSolution);
        double currentTardiness = calculateMaxTardiness(currentSolution, processingTimes, dueDates, delayMatrix);
        double bestTardiness = currentTardiness;
        double temperature = initialTemperature;
        Random rand = new Random();

        for (int iteration = 0; iteration < maxIterationsSA; iteration++) {
            // Generate neighbors without causing infinite loops
            List<Integer> neighborForward = new ArrayList<>(currentSolution);
            List<Integer> neighborBackward = new ArrayList<>(currentSolution);

            // Forward shift: move job to a later position
            if (numJobs > 1) {
                int idx1 = rand.nextInt(numJobs - 1); // Ensures idx1 can have a valid idx2 > idx1
                int idx2 = idx1 + 1 + rand.nextInt(numJobs - idx1 - 1);
                Integer job = neighborForward.remove(idx1);
                neighborForward.add(idx2, job);
            }

            // Backward shift: move job to an earlier position
            if (numJobs > 1) {
                int idx3 = 1 + rand.nextInt(numJobs - 1); // Ensures idx3 > 0
                int idx4 = rand.nextInt(idx3);
                Integer job2 = neighborBackward.remove(idx3);
                neighborBackward.add(idx4, job2);
            }




            // Evaluate neighbors
            double tardinessForward = calculateMaxTardiness(neighborForward, processingTimes, dueDates, delayMatrix);
            double tardinessBackward = calculateMaxTardiness(neighborBackward, processingTimes, dueDates, delayMatrix);

            // Determine the best neighbor
            List<Integer> bestNeighbor;
            double bestNeighborTardiness;
            if (tardinessForward < tardinessBackward) {
                bestNeighbor = neighborForward;
                bestNeighborTardiness = tardinessForward;
            } else {
                bestNeighbor = neighborBackward;
                bestNeighborTardiness = tardinessBackward;
            }

            // Metropolis acceptance criterion
            double delta = bestNeighborTardiness - currentTardiness;
            if (delta < 0 || Math.exp(-delta / temperature) > rand.nextDouble()) {
                currentSolution = bestNeighbor;
                currentTardiness = bestNeighborTardiness;

                if (currentTardiness < bestTardiness) {
                    bestSolution = new ArrayList<>(currentSolution);
                    bestTardiness = currentTardiness;
                }
            }

            // Cool down temperature
            temperature *= coolingRate;
        }

        return bestSolution;
    }
}