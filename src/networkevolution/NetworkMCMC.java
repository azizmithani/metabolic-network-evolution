/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package networkevolution;

import cern.jet.random.Beta;
import cern.jet.random.Gamma;
import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import network.DirectedReaction;
import network.MetabolicNetwork;
import systemobject.PhyloTree;
import rahnumautilities.Utilities;
import systemobject.PhyloNode;

/**
 *
 * @author mithani
 */
public class NetworkMCMC {

    // Constants
    static public final int NO_OF_ALLOWED_MOVES = 4;
    static public final byte MOVE_ADD_EVENTS = 0;
    static public final byte MOVE_SWAP_EVENTS = 1;
    static public final byte MOVE_DEL_EVENTS = 2;
    static public final byte MOVE_PERMUTE_EVENTS = 3;
    static private final int PRECISION = 128;
    // MCMC Modes
    static public final int MCMC_MODE_PATH = (int) Math.pow(2, 0);
    static public final int MCMC_MODE_RATES = (int) Math.pow(2, 1);
    static public final int MCMC_MODE_BOTH = MCMC_MODE_PATH + MCMC_MODE_RATES;
    // priavate variables to hold data values that require time consuming calculations
    private BigDecimal natural_e;
    private BigDecimal epsilon;
    private BigDecimal[] factorial;
    // variables to hold results
//    private ArrayList pathList;
//    private BigDecimal[] pathProbList;
//    private ArrayList networkList;
//    private Double[] insRates;
//    private Double[] delRates;

    // variables to hold intermediary data
    private HashMap distinctPathList;
    private HashMap distinctNetworks;
    private HashMap distinctNetworkPairs;
    private HashMap distinctTrees;
    private HashMap mapTransProb;
    // variable to decide if output should be printed
    private boolean verbose = true;
    private boolean debug = true;
    private NetworkEvolver networkEvolver;

    // variable to generate random numbers from gamma distribution
    private Gamma gammaDist;
    // variable to generate random numbers from beta distribution
    private Beta betaDist;

    // for proposal sampling
    int decimalPlaces = 4;

    public NetworkMCMC() {
        natural_e = Utilities.naturalE(PRECISION);
        epsilon = Utilities.epsilon(PRECISION);
        factorial = Utilities.factorial(PRECISION);

        networkEvolver = new NetworkEvolver();

        betaDist = new Beta(1, 1, RandomEngine.makeDefault());
        gammaDist = new Gamma(1, 1, RandomEngine.makeDefault());
    }

    public NetworkEvolver getNetworkEvolver() {

        if (this.networkEvolver == null) {
            this.networkEvolver = new NetworkEvolver();
        }
        return this.networkEvolver;
    }

    public HashMap getDistinctNetworkPairs() {
        return distinctNetworkPairs;
    }

    public void setDistinctNetworkPairs(HashMap distinctNetworkPairs) {
        this.distinctNetworkPairs = distinctNetworkPairs;
    }

    public HashMap getDistinctTrees() {
        return distinctTrees;
    }

    public void setDistinctTrees(HashMap distinctTrees) {
        this.distinctTrees = distinctTrees;
    }

//    public Double[] getDelRates() {
//        return delRates;
//    }
//
//    public Double[] getInsRates() {
//        return insRates;
//    }
//
//    public ArrayList getNetworkList() {
//        return networkList;
//    }
//
//    public ArrayList getPathList() {
//        return pathList;
//    }
//
//    public BigDecimal[] getPathProbList() {
//        return pathProbList;
//    }
    public MCMCOutput networkMCMC(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, double insRate, double delRate, double evolTime, int nIter,
            int nBurning, boolean enumStates, int updateInterval, int evolutionMode, boolean sampleRates,
            double dependenceProbability, boolean sampleDependenceProbability) {

        String strPath;
        BigDecimal newPathProb;

        // keep track of number of moves accepted during MCMC
        int acceptanceCountPath = 0;
        int acceptanceCountParameter = 0;
        int acceptanceCountDP = 0;

        // initialise the variable
        MCMCOutput mcmcOutput = new MCMCOutput(nIter);

        ArrayList lstEdges = startNetwork.getDirectedReactions();
        // total number of edges
        int edgeCount = lstEdges.size();

        // find the number of differences between the two networks
        int nDifferences = MetabolicNetwork.findDifferences(startNetwork, endNetwork).size();

        // generate a path
        int[] path = generatePath(startNetwork, endNetwork, coreEdges, prohibEdges, evolutionMode, insRate, delRate, dependenceProbability);
        // Calculate the path probability
        BigDecimal pathProb = pathProbability(path, startNetwork, coreEdges, prohibEdges,
                insRate, delRate, evolTime, evolutionMode, dependenceProbability);

        // add the starting state to the list
        mcmcOutput.getPaths().add(0, path);
        mcmcOutput.getPathProbabilities()[0] = pathProb;
        mcmcOutput.getInsRates()[0] = insRate;
        mcmcOutput.getDelRates()[0] = delRate;

        // keep track of distinct paths visited so that we dont recalculate the probabilities
        HashMap lstDistinctPaths = new HashMap((int) (nIter * 0.5)); // Assumption: Only half of the time you see new path
        // keep track of paths that we have already processed
        HashSet processedPathList = new HashSet((int) (nIter * 0.4)); // Assumption: out of the times you see new path, not all of them are accepted
        // also keep track of networks visited 
        HashSet lstDistinctNetworks = new HashSet(20000); // initial capacity ?!?

        double[] insParameters = new double[]{1, edgeCount};
        double[] delParameters = new double[]{1, edgeCount};

        double dpAlpha = 1.0;
        double dpBeta = 1.0;
        if (sampleDependenceProbability) {
            double startNC = (double) Utilities.sum(startNetwork.getNeighboursCount()) / (double) edgeCount;
            double endNC = (double) Utilities.sum(endNetwork.getNeighboursCount()) / (double) edgeCount;

            Byte[] refNW = new Byte[edgeCount];
            for (int i = 0; i < refNW.length; i++) {
                refNW[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            }
            double refNC = (double) Utilities.sum(startNetwork.getNetwork(refNW).getNeighboursCount()) / (double) edgeCount;
            dpAlpha = (startNC + endNC) / 2 + 1;
            dpBeta = refNC + 1;
            betaDist = new Beta(dpAlpha, dpBeta, RandomEngine.makeDefault());
        }

        // add the starting state to the list
        lstDistinctPaths.put(getPathRatesString(path, insRate, delRate, dependenceProbability), pathProb);
        for (int iter = 1; iter <= nIter; iter++) {

            if (updateInterval > 0 && iter % updateInterval == 0) {
                System.out.print("Iteration: " + Integer.toString(iter));
                System.out.print("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate) + "\t" + Double.toString(dependenceProbability));
                System.out.println("\t" + Utilities.toString(path) + "\t" + pathProb.round(MathContext.DECIMAL128).toString());
            }

            if (sampleRates) {
                // data structures to hold insertion and deletion proposal parameters
                double[] newInsParameters = new double[2];
                double[] newDelParameters = new double[2];
                // sample new inserstion / deletion rates
                double[] newRates = proposeNewRates(path, startNetwork, coreEdges, prohibEdges,
                        newInsParameters, newDelParameters, evolutionMode);
                double newInsRate = newRates[0];
                double newDelRate = newRates[1];

                //System.out.println("Iteration: " + Integer.toString(iter) + "\tProposal:\t" + Double.toString(newInsRate) + "\t" + Double.toString(newDelRate));

                double newDependenceProbability;
                if (sampleDependenceProbability) {
                    newDependenceProbability = Utilities.round(betaDist.nextDouble(), decimalPlaces);
                } else {
                    newDependenceProbability = dependenceProbability;
                }

                strPath = getPathRatesString(path, newInsRate, newDelRate, newDependenceProbability);
                newPathProb = (BigDecimal) lstDistinctPaths.get(strPath);
                if (newPathProb == null) {
                    // Calculate the path probability
                    newPathProb = pathProbability(path, startNetwork, coreEdges, prohibEdges,
                            newInsRate, newDelRate, evolTime, evolutionMode, newDependenceProbability);
                    // add it to the list of distinct paths
                    lstDistinctPaths.put(strPath, newPathProb);
                }

                // calculate the acceptance probability
                BigDecimal alphaParameters = calculateParameterAcceptanceProbability(insParameters, delParameters,
                        newInsParameters, newDelParameters, insRate, delRate, newInsRate, newDelRate,
                        pathProb, newPathProb, dpAlpha, dpBeta, dependenceProbability, newDependenceProbability);

                // accept the new state with probability accept
                double acceptParameters = Math.min(1.0, alphaParameters.doubleValue());
                double uAcceptParameters = Math.random();

                if (uAcceptParameters <= acceptParameters) {
                    // set the new state
                    insRate = newInsRate;
                    delRate = newDelRate;
                    insParameters = newInsParameters;
                    delParameters = newDelParameters;
                    dependenceProbability = newDependenceProbability;

                    // also copy the new path probability
                    pathProb = newPathProb;

                    if (iter > nBurning) {
                        // increment the acceptance count
                        acceptanceCountDP++;
                        acceptanceCountParameter++;
                    }
                }
            }// end if sample rates
            // save the rates for the current iteration
            mcmcOutput.getInsRates()[iter] = insRate;
            mcmcOutput.getDelRates()[iter] = delRate;
            mcmcOutput.getDependenceProbabilities()[iter] = dependenceProbability;

            /*            if (sampleDependenceProbability) {
            double newDependenceProbability = Utilities.round(betaDist.nextDouble(), decimalPlaces);

            strPath = getPathRatesString(path, insRate, delRate, newDependenceProbability);
            newPathProb = (BigDecimal) lstDistinctPaths.get(strPath);
            if (newPathProb == null) {
            // Calculate the path probability
            newPathProb = pathProbability(path, startNetwork, coreEdges, prohibEdges,
            insRate, delRate, evolTime, evolutionMode, newDependenceProbability);
            // add it to the list of distinct paths
            lstDistinctPaths.put(strPath, newPathProb);
            }

            // calculate the acceptance probability
            BigDecimal alphaDP = calculateDPAcceptanceProbability(dpAlpha, dpBeta,
            dependenceProbability, newDependenceProbability, pathProb, newPathProb);

            // accept the new state with probability accept
            double acceptDP = Math.min(1.0, alphaDP.doubleValue());
            double uAcceptDP = Math.random();

            if (uAcceptDP <= acceptDP) {
            if (iter > nBurning) {// && dependenceProbability != newDependenceProbability) {
            // increment the acceptance count
            acceptanceCountDP++;
            }
            // set the new state
            dependenceProbability = newDependenceProbability;

            // also copy the new path probability
            pathProb = newPathProb;

            }
            }// end if sample rates
            // save the dp for the current iteration
            mcmcOutput.getDependenceProbabilities()[iter] = dependenceProbability;
             */

            //mcmcDiagnostics(mcmcOutput, iter, nBurning);


            // sample a new path
            int[] newPath = proposeNewPath(path, startNetwork, coreEdges, prohibEdges,
                    nDifferences);

            strPath = getPathRatesString(newPath, insRate, delRate, dependenceProbability);
            newPathProb = (BigDecimal) lstDistinctPaths.get(strPath);
            if (newPathProb == null) {
                // Calculate the path probability
                newPathProb = pathProbability(newPath, startNetwork, coreEdges, prohibEdges,
                        insRate, delRate, evolTime, evolutionMode, dependenceProbability);
                // add it to the list of distinct paths
                lstDistinctPaths.put(strPath, newPathProb);
            }
            // calculate the acceptance probability
            BigDecimal alpha = calculatePathAcceptanceProbability(path, newPath, edgeCount,
                    nDifferences, pathProb, newPathProb);

            // accept the new state with probability accept
            double accept = Math.min(1.0, alpha.doubleValue());
            double uAccept = Math.random();

            /*            System.out.print(path.length);
            System.out.print("\t" + pathProb.toString() + "\t");
            System.out.print(newPath.length);
            System.out.print("\t" + newPathProb.toString());
            System.out.println("\t" + accept.toString());
            System.out.print("\t" + Utilities.toString(path, ", "));
            System.out.println("\t" + Utilities.toString(newPath, ", "));
             */
            if (uAccept <= accept) {

                path = newPath;
                pathProb = newPathProb;

                // Throw away the first N entries specified by nBurning
                if (iter > nBurning) {
                    // incrememnt the acceptance count
                    acceptanceCountPath++;

                    // If we seen not this path before then process it
                    if (enumStates && processedPathList.add(strPath)) {
                        MetabolicNetwork theNetwork = startNetwork.clone();
                        for (int i = 0; i < path.length - 1; i++) { // dont need to process the last event. it will always lead to destination network
                            // get the edge corresponding to ith event
                            DirectedReaction reaction = (DirectedReaction) lstEdges.get(path[i] - 1);

                            // Update the network
                            theNetwork = theNetwork.updateNetwork(reaction);

                            // If we have not seen this network before then add it
                            // to the list
                            String strNetwork = Utilities.toString(theNetwork.getReactionSequence());
                            if (lstDistinctNetworks.add(strNetwork)) // if successfully added then it's a new network
                            {
                                mcmcOutput.getVisitedNetworks().add(theNetwork.getReactionSequence());
                            }
                        } // end for each event in the path
                    } // if enumerate states                
                } // if currentIter >burning
            } // if accepted

            // add the new state to the list
            mcmcOutput.getPaths().add(iter, path);
            mcmcOutput.getPathProbabilities()[iter] = pathProb;

        } // end for currentIter

        mcmcOutput.getVisitedNetworks().trimToSize();
        mcmcOutput.setAcceptCountPath(acceptanceCountPath);
        mcmcOutput.setAcceptCountParameter(acceptanceCountParameter);
        mcmcOutput.setAcceptCountDependenceProbability(acceptanceCountDP);
        return mcmcOutput;

    }

    public void networkMCMC(PhyloTree phylotree, MetabolicNetwork refNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, double insRate, double delRate, int nIter, int updateInterval,
            int evolutionModel, int nBurning, boolean sampleRates, int subNetworkSize,
            double dependenceProbability, boolean sampleDependenceProbability) {

        debug = false;
        // Initialise the intermediary variables
//        distinctPathList = new HashMap((int) (nIter * 0.5)); // Assumption: Only half of the time you see new path
        distinctNetworkPairs = new HashMap((int) (nIter * 0.2)); // Assumption: Number of networks visited << number of iterations
        mapTransProb = new HashMap((int) (nIter * 0.2));

        // acceptance count
        int acceptanceCountParameter = 0;
        int acceptanceCountDP = 0;

        // Initialise the MCMC output variables at the internal nodes of the phylogeny
        initialiseMCMCOutput(phylotree.getRoot(), nIter - 1);

        // set initial states
        MCMCStep_Initial(phylotree.getRoot(), refNetwork, insRate, delRate);

        // get indices of core and prohibited edges
        Iterator itCore = coreEdges.iterator();
        Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iCoreEdges.length; i++) {
            iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itCore.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itCore.next();
            iCoreEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }
        Iterator itProhib = prohibEdges.iterator();
        Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iProhibEdges.length; i++) {
            iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itProhib.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itProhib.next();
            iProhibEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }

        // get a local reference to the rate parameters
        MCMCOutput rootMCMCOutput = (MCMCOutput) phylotree.getRoot().getData();
        Double[] insRates = rootMCMCOutput.getInsRates();
        Double[] delRates = rootMCMCOutput.getDelRates();
        Double[] dependenceProbabilities = rootMCMCOutput.getDependenceProbabilities();

        double[] insParameters = new double[]{1, iCoreEdges.length};
        double[] delParameters = new double[]{1, iCoreEdges.length};

        double dpAlpha = 1.0;
        double dpBeta = 1.0;
        if (sampleDependenceProbability) {
            double phyloNC = getAverageNeighboursCount(phylotree, iCoreEdges.length);

            double refNC = (double) Utilities.sum(refNetwork.getNeighboursCount()) / (double) iCoreEdges.length;
            dpAlpha = phyloNC + 1;
            dpBeta = refNC + 1;
            betaDist = new Beta(dpAlpha, dpBeta, RandomEngine.makeDefault());
        }

        // Calculate the tree probability
        BigDecimal treeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(), refNetwork,
                iCoreEdges, iProhibEdges, insRate, delRate, evolutionModel, subNetworkSize, dependenceProbability)));

        // save the rates for the current iteration
        insRates[0] = insRate;
        delRates[0] = delRate;

        if (updateInterval > 0) {
            System.out.print("Iteration: 0");
            System.out.println("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate) + "\t" + Double.toString(dependenceProbability));
        }

        double newInsRate, newDelRate, newDependenceProbability;
        int nInternalNodes = phylotree.getLeaves().size() - 1;
        int alterableEdges = refNetwork.getReactionSequence().length - Utilities.sum(iCoreEdges) - Utilities.sum(iProhibEdges);
        // set the parameters
        double[] newInsParameters = new double[2];
        double[] newDelParameters = new double[2];

        for (int iter = 1; iter <= nIter; iter++) {

            if (sampleRates || sampleDependenceProbability) {
                if (debug) {
                    System.out.println("> Sampling rates ...");
                }

                if (sampleRates) {
                    // set the parameters
                    newInsParameters = new double[2];
                    newDelParameters = new double[2];
                    // sample new inserstion / deletion rates
//                double[] newRates = sampleRates(phylotree.getRoot(), refNetwork, iCoreEdges, iProhibEdges,
//                        evolutionModel, insParameters, delParameters);
                    double[] newRates = sampleRates(phylotree.getRoot(), refNetwork, iCoreEdges, iProhibEdges,
                            evolutionModel, newInsParameters, newDelParameters, nInternalNodes, alterableEdges);
                    newInsRate = newRates[0];
                    newDelRate = newRates[1];
                } else {
                    newInsRate = insRate;
                    newDelRate = delRate;
                }
//                System.out.println(Double.toString(newInsRate) + "\t" + Double.toString(newDelRate));
                if (sampleDependenceProbability) {
                    newDependenceProbability = Utilities.round(betaDist.nextDouble(), 5);
                } else {
                    newDependenceProbability = dependenceProbability;
                }

                String strTree;
                BigDecimal newTreeProb = (BigDecimal) distinctTrees.get(strTree = getTreeRatesString(phylotree, newInsRate, newDelRate, newDependenceProbability));
                if (newTreeProb == null) {
                    newTreeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
                            refNetwork, iCoreEdges, iProhibEdges, newInsRate, newDelRate,
                            evolutionModel, subNetworkSize, newDependenceProbability)));

                    distinctTrees.put(strTree, newTreeProb);

                }

                //System.out.println ("\tProposal: " + "\t" + newInsRate + "\t" + newDelRate + "\t" + newTreeProb.round(new MathContext(8)).toString());

                double accept, uAccept;
                if (newInsRate == 0 || newDelRate == 0) {
                    uAccept = 1.0;
                    accept = 0.0;
                } else {

                    //System.out.println(">\t" + treeProb.toString() + "\t" + newTreeProb.toString());
                    // calculate the acceptance probability
                    BigDecimal alpha = calculateParameterAcceptanceProbabilityPhylo(insParameters,
                            delParameters, newInsParameters, newDelParameters, insRate, delRate, newInsRate,
                            newDelRate, treeProb, newTreeProb, dpAlpha, dpBeta, dependenceProbability,
                            newDependenceProbability);

                    //System.out.println(">\t" + alphaRates.toString());

                    // accept the new state with probability accept
                    accept = Math.min(1.0, alpha.doubleValue());
                    uAccept = Math.random();
                }


                if (uAccept <= accept) {
                    if (iter > nBurning) {// && dependenceProbability != newDependenceProbability) {
                        // increment the acceptance count
                        acceptanceCountDP++;
                        acceptanceCountParameter++;
                    }
                    // save the new state
                    insRate = newInsRate;
                    delRate = newDelRate;
                    insParameters = newInsParameters;
                    delParameters = newDelParameters;
                    dependenceProbability = newDependenceProbability;

                    // also copy the new tree probability
                    treeProb = newTreeProb;

                }

            }// end if sample rates or sample DP

            // save the rates for the current iteration
            insRates[iter] = insRate;
            delRates[iter] = delRate;
            // save the dp for the current iteration
            dependenceProbabilities[iter] = dependenceProbability;

            boolean isDebug = debug;
            if (updateInterval > 0 && iter % updateInterval == 0) {
                System.out.print("Iteration: " + Integer.toString(iter));
                System.out.println("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate) + "\t" + Double.toString(dependenceProbability));// + "\t" + treeProb.round(new MathContext(8)).toString());
                debug = true;
            }

            // sample networks at internal nodes
            if (false) {
                System.out.println("> Sampling networks at internal nodes of the tree...");
            }

            sampleInternalNetworks(phylotree.getRoot(), refNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, dependenceProbability);

            // calculate the probability for the new tree
            if (sampleRates || sampleDependenceProbability) {
                String strTree;
                treeProb = (BigDecimal) distinctTrees.get(strTree = getTreeRatesString(phylotree, insRate, delRate, dependenceProbability));
                if (treeProb == null) {
                    treeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
                            refNetwork, iCoreEdges, iProhibEdges, insRate, delRate,
                            evolutionModel, subNetworkSize, dependenceProbability)));

                    distinctTrees.put(strTree, treeProb);

                }
            }

            debug = isDebug;

            if (debug) {
                System.out.println();
            }
        } // end for iter = 1 to nIter

        // set the parameter accepatance count
        ((MCMCOutput) phylotree.getRoot().getData()).setAcceptCountParameter(acceptanceCountParameter);
        ((MCMCOutput) phylotree.getRoot().getData()).setAcceptCountDependenceProbability(acceptanceCountDP);
        ((MCMCOutput) phylotree.getRoot().getData()).setDependenceProbabilities(dependenceProbabilities);

    }


    /*
    public void networkMCMC(PhyloTree phylotree, MetabolicNetwork refNetwork, ArrayList coreEdges,
    ArrayList prohibEdges, double insRate, double delRate, int nIter, int updateInterval,
    int evolutionModel, int nBurning, boolean sampleRates, int subNetworkSize,
    double dependenceProbability, boolean sampleDependenceProbability) {

    debug = false;
    // Initialise the intermediary variables
    //        distinctPathList = new HashMap((int) (nIter * 0.5)); // Assumption: Only half of the time you see new path
    distinctNetworkPairs = new HashMap((int) (nIter * 0.2)); // Assumption: Number of networks visited << number of iterations
    mapTransProb = new HashMap((int) (nIter * 0.2));

    // acceptance count
    int acceptanceCountParameter = 0;
    int acceptanceCountDP = 0;

    // Initialise the MCMC output variables at the internal nodes of the phylogeny
    initialiseMCMCOutput(phylotree.getRoot(), nIter - 1);

    // set initial states
    MCMCStep_Initial(phylotree.getRoot(), refNetwork, insRate, delRate);

    // get indices of core and prohibited edges
    Iterator itCore = coreEdges.iterator();
    Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
    for (int i = 0; i < iCoreEdges.length; i++) {
    iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
    }
    while (itCore.hasNext()) {
    DirectedReaction reaction = (DirectedReaction) itCore.next();
    iCoreEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
    }
    Iterator itProhib = prohibEdges.iterator();
    Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
    for (int i = 0; i < iProhibEdges.length; i++) {
    iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
    }
    while (itProhib.hasNext()) {
    DirectedReaction reaction = (DirectedReaction) itProhib.next();
    iProhibEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
    }

    // get a local reference to the rate parameters
    MCMCOutput rootMCMCOutput = (MCMCOutput) phylotree.getRoot().getData();
    Double[] insRates = rootMCMCOutput.getInsRates();
    Double[] delRates = rootMCMCOutput.getDelRates();
    Double[] dependenceProbabilities = rootMCMCOutput.getDependenceProbabilities();

    double[] insParameters = new double[]{1, iCoreEdges.length};
    double[] delParameters = new double[]{1, iCoreEdges.length};

    double dpAlpha = 1.0;
    double dpBeta = 1.0;
    if (sampleDependenceProbability) {
    double phyloNC = getAverageNeighboursCount(phylotree, iCoreEdges.length);

    double refNC = (double) Utilities.sum(refNetwork.getNeighboursCount()) / (double) iCoreEdges.length;
    dpAlpha = phyloNC + 1;
    dpBeta = refNC + 1;
    betaDist = new Beta(dpAlpha, dpBeta, RandomEngine.makeDefault());
    }

    // Calculate the tree probability
    BigDecimal treeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(), refNetwork,
    iCoreEdges, iProhibEdges, insRate, delRate, evolutionModel, subNetworkSize, dependenceProbability)));

    // save the rates for the current iteration
    insRates[0] = insRate;
    delRates[0] = delRate;

    if (updateInterval > 0) {
    System.out.print("Iteration: 0");
    System.out.println("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate) + "\t" + Double.toString(dependenceProbability));
    }

    int nInternalNodes = phylotree.getLeaves().size() - 1;
    for (int iter = 1; iter <= nIter; iter++) {

    if (sampleRates) {
    if (debug) {
    System.out.println("> Sampling rates ...");
    }

    // set the parameters
    double[] newInsParameters = new double[2];
    double[] newDelParameters = new double[2];
    // sample new inserstion / deletion rates
    //                double[] newRates = sampleRates(phylotree.getRoot(), refNetwork, iCoreEdges, iProhibEdges,
    //                        evolutionModel, insParameters, delParameters);
    double[] newRates = sampleRates(phylotree.getRoot(), refNetwork, iCoreEdges, iProhibEdges,
    evolutionModel, newInsParameters, newDelParameters, nInternalNodes, insRate, delRate);
    double newInsRate = newRates[0];
    double newDelRate = newRates[1];

    //                System.out.println(Double.toString(newInsRate) + "\t" + Double.toString(newDelRate));

    String strTree;
    BigDecimal newTreeProb = (BigDecimal) distinctTrees.get(strTree = getTreeRatesString(phylotree, newInsRate, newDelRate, dependenceProbability));
    if (newTreeProb == null) {
    newTreeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
    refNetwork, iCoreEdges, iProhibEdges, newInsRate, newDelRate,
    evolutionModel, subNetworkSize, dependenceProbability)));

    distinctTrees.put(strTree, newTreeProb);

    }

    //System.out.println ("\tProposal: " + "\t" + newInsRate + "\t" + newDelRate + "\t" + newTreeProb.round(new MathContext(8)).toString());

    double acceptRates, uAcceptRates;
    if (newInsRate == 0 || newDelRate == 0) {
    uAcceptRates = 1.0;
    acceptRates = 0.0;
    } else {

    //System.out.println(">\t" + treeProb.toString() + "\t" + newTreeProb.toString());
    // calculate the acceptance probability
    BigDecimal alphaRates = calculateRateAcceptanceProbabilityPhylo(insParameters, delParameters,
    newInsParameters, newDelParameters,
    insRate, delRate, newInsRate, newDelRate, treeProb, newTreeProb);

    //System.out.println(">\t" + alphaRates.toString());

    // accept the new state with probability accept
    acceptRates = Math.min(1.0, alphaRates.doubleValue());
    uAcceptRates = Math.random();
    }

    if (uAcceptRates < acceptRates) {
    // save the new state
    insRate = newInsRate;
    delRate = newDelRate;
    insParameters = newInsParameters;
    delParameters = newDelParameters;

    // also copy the new tree probability
    treeProb = newTreeProb;

    if (iter > nBurning) {
    // increment the acceptance count
    acceptanceCountParameter++;
    }

    }

    // save the rates for the current iteration
    insRates[iter] = insRate;
    delRates[iter] = delRate;

    }

    if (sampleDependenceProbability) {
    double newDependenceProbability = Utilities.round(betaDist.nextDouble(), 5);

    String strTree;
    BigDecimal newTreeProb = (BigDecimal) distinctTrees.get(strTree = getTreeRatesString(phylotree, insRate, delRate, newDependenceProbability));
    if (newTreeProb == null) {
    // Calculate the tree probability
    newTreeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
    refNetwork, iCoreEdges, iProhibEdges, insRate, delRate,
    evolutionModel, subNetworkSize, newDependenceProbability)));

    // add it to the list of distinct tree
    distinctTrees.put(strTree, newTreeProb);
    }

    // calculate the acceptance probability
    BigDecimal alphaDP = calculateDPAcceptanceProbability(dpAlpha, dpBeta,
    dependenceProbability, newDependenceProbability, treeProb, newTreeProb);

    // accept the new state with probability accept
    double acceptDP = Math.min(1.0, alphaDP.doubleValue());
    double uAcceptDP = Math.random();

    if (uAcceptDP <= acceptDP) {
    if (iter > nBurning) {// && dependenceProbability != newDependenceProbability) {
    // increment the acceptance count
    acceptanceCountDP++;
    }
    // set the new state
    dependenceProbability = newDependenceProbability;

    // also copy the new path probability
    treeProb = newTreeProb;

    }
    }// end if sample rates
    // save the dp for the current iteration
    dependenceProbabilities[iter] = dependenceProbability;

    boolean isDebug = debug;
    if (updateInterval > 0 && iter % updateInterval == 0) {
    System.out.print("Iteration: " + Integer.toString(iter));
    System.out.println("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate) + "\t" + Double.toString(dependenceProbability));// + "\t" + treeProb.round(new MathContext(8)).toString());
    debug = true;
    }

    // sample networks at internal nodes
    if (false) {
    System.out.println("> Sampling networks at internal nodes of the tree...");
    }

    sampleInternalNetworks(phylotree.getRoot(), refNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, dependenceProbability);

    // calculate the probability for the new tree
    if (sampleRates || sampleDependenceProbability) {
    String strTree;
    treeProb = (BigDecimal) distinctTrees.get(strTree = getTreeRatesString(phylotree, insRate, delRate, dependenceProbability));
    if (treeProb == null) {
    treeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
    refNetwork, iCoreEdges, iProhibEdges, insRate, delRate,
    evolutionModel, subNetworkSize, dependenceProbability)));

    distinctTrees.put(strTree, treeProb);

    }
    }

    debug = isDebug;

    if (debug) {
    System.out.println();
    }
    } // end for iter = 1 to nIter

    // set the parameter accepatance count
    ((MCMCOutput) phylotree.getRoot().getData()).setAcceptCountParameter(acceptanceCountParameter);
    ((MCMCOutput) phylotree.getRoot().getData()).setAcceptCountDependenceProbability(acceptanceCountDP);
    ((MCMCOutput) phylotree.getRoot().getData()).setDependenceProbabilities(dependenceProbabilities);

    }

     */
    private void initialiseMCMCOutput(PhyloNode phylonode, int nIter) {

        if (!phylonode.isLeaf()) {

            // initilise left and right sub-trees
            initialiseMCMCOutput(phylonode.getLeftSon(), nIter);
            initialiseMCMCOutput(phylonode.getRightSon(), nIter);

            // initilise the output object at the current node
            phylonode.setData(new MCMCOutput(nIter + 1));
        }
    }

    private void MCMCStep_Initial(PhyloNode phylonode, MetabolicNetwork refNetwork, double insRate, double delRate) {

        if (!phylonode.getLeftSon().isLeaf()) // run MCMC Step on left sub tree
        {
            MCMCStep_Initial(phylonode.getLeftSon(), refNetwork, insRate, delRate);
        }
        if (!phylonode.getRightSon().isLeaf()) // run MCMC Step on right sub tree
        {
            MCMCStep_Initial(phylonode.getRightSon(), refNetwork, insRate, delRate);
        }

        // get left & right sons' networks
        MetabolicNetwork leftNetwork, rightNetwork;
        if (phylonode.getLeftSon().isLeaf()) {
            leftNetwork = phylonode.getLeftSon().getMetabolicNetwork();
        } else {
            MCMCOutput leftMCMCOutput = (MCMCOutput) phylonode.getLeftSon().getData();
            leftNetwork = leftMCMCOutput.getCurrentNetwork();
        }
        if (phylonode.getRightSon().isLeaf()) {
            rightNetwork = phylonode.getRightSon().getMetabolicNetwork();
        } else {
            MCMCOutput rightMCMCOutput = (MCMCOutput) phylonode.getRightSon().getData();
            rightNetwork = rightMCMCOutput.getCurrentNetwork();
        }

        // Assign initial value to the ancestral network as follows:
        // 1. If the edge is present (or abesnt) in both then mark it as present (or absent)
        // 2. For different edge, randomly select 0 or 1
        Byte[] ancestralNetworkSeq = new Byte[leftNetwork.getReactionSequence().length];
        Byte[] leftNetworkSeq = leftNetwork.getReactionSequence();
        Byte[] rightNetworkSeq = rightNetwork.getReactionSequence();
        for (int i = 0; i < ancestralNetworkSeq.length; i++) {
            if (leftNetworkSeq[i].byteValue() == rightNetworkSeq[i].byteValue()) // edge has same status
            {
                ancestralNetworkSeq[i] = leftNetworkSeq[i];
            } else {
                int edgeStatus = Utilities.randsample(new double[]{0.5, 0.5});
                ancestralNetworkSeq[i] = (edgeStatus == 0 ? MetabolicNetwork.SEQ_ENTRY_ABSENT : MetabolicNetwork.SEQ_ENTRY_PRESENT);
            }
        }

        MetabolicNetwork ancestralNetwork = refNetwork.getNetwork(ancestralNetworkSeq);
        ancestralNetwork.setID("(" + leftNetwork.getID() + "," + rightNetwork.getID() + ")");
        ancestralNetwork.setName("(" + leftNetwork.getName() + ", " + rightNetwork.getName() + ")");
        // set the initial state
        MCMCOutput mcmcOutput = (MCMCOutput) phylonode.getData();
        // set the current network
        mcmcOutput.setCurrentNetwork(ancestralNetwork);

        // set the state for this iteration
        mcmcOutput.getNetworks().add(ancestralNetworkSeq);
        if (phylonode.getParent() == null) { // save the rates at root only.
            mcmcOutput.setCurrentRates(new double[]{insRate, delRate});
//            mcmcOutput.getInsRates()[0] = insRate;
//            mcmcOutput.getDelRates()[0] = delRate;
        }

    }

    private void sampleInternalNetworks(PhyloNode phylonode, MetabolicNetwork refNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int evolutionModel, double insRate, double delRate, double dependenceProbability) {

        if (!phylonode.getLeftSon().isLeaf()) // run MCMC Step on left sub tree
        {
            sampleInternalNetworks(phylonode.getLeftSon(), refNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, dependenceProbability);
        }
        if (!phylonode.getRightSon().isLeaf()) // run MCMC Step on right sub tree
        {
            sampleInternalNetworks(phylonode.getRightSon(), refNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, dependenceProbability);
        }

        // get a local copy of the network evolver
        NetworkEvolver nwEvolver = getNetworkEvolver();
        // get left & right sons' networks
        Byte[] leftNetworkSeq, rightNetworkSeq;
        if (phylonode.getLeftSon().isLeaf()) {
            leftNetworkSeq = phylonode.getLeftSon().getMetabolicNetwork().getReactionSequence();
        } else {
            MCMCOutput leftMCMCOutput = (MCMCOutput) phylonode.getLeftSon().getData();
            leftNetworkSeq = leftMCMCOutput.getCurrentNetwork().getReactionSequence();
        }
        if (phylonode.getRightSon().isLeaf()) {
            rightNetworkSeq = phylonode.getRightSon().getMetabolicNetwork().getReactionSequence();
        } else {
            MCMCOutput rightMCMCOutput = (MCMCOutput) phylonode.getRightSon().getData();
            rightNetworkSeq = rightMCMCOutput.getCurrentNetwork().getReactionSequence();
        }
        // get parent's network
        Byte[] parentNetworkSeq;
        if (phylonode.getParent() != null) {
            parentNetworkSeq = ((MCMCOutput) phylonode.getParent().getData()).getCurrentNetwork().getReactionSequence();
        } else {
            parentNetworkSeq = null;
        }

        // get evol times (branch lengths)
        double leftEvolTime, rightEvolTime, parentEvolTime;
        leftEvolTime = phylonode.getLeftBranchLength();
        rightEvolTime = phylonode.getRightBranchLength();
        if (phylonode.getParent() == null) {
            parentEvolTime = -1.0;
        } else if (phylonode.isLeftSon()) {
            parentEvolTime = phylonode.getParent().getLeftBranchLength();
        } else {
            parentEvolTime = phylonode.getParent().getRightBranchLength();
        }

        // add the networks and times in arraylists for iteration purpose 
        ArrayList neighbourNetworks = new ArrayList(3);
        ArrayList neighbourEvolTimes = new ArrayList(3);
        // left
        neighbourNetworks.add(leftNetworkSeq);
        neighbourEvolTimes.add(leftEvolTime);
        // right
        neighbourNetworks.add(rightNetworkSeq);
        neighbourEvolTimes.add(rightEvolTime);
        if (parentNetworkSeq != null) {
            neighbourNetworks.add(parentNetworkSeq);
            neighbourEvolTimes.add(parentEvolTime);
        }

        // sample the new network
        MCMCOutput mcmcOutput = (MCMCOutput) phylonode.getData();
        MetabolicNetwork currentNetwork = mcmcOutput.getCurrentNetwork();
        Iterator it = currentNetwork.getDirectedReactions().iterator();
        Byte[] newNetworkSeq = new Byte[currentNetwork.getReactionSequence().length];
        int idx = -1;
        while (it.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) it.next();
            idx++;

            if (coreEdges.contains(reaction)) { // core edge .. has to be present
                newNetworkSeq[idx] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            } else if (prohibEdges.contains(reaction)) { // prohibited edge .. has to be absent
                newNetworkSeq[idx] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            } else {
                // sameple the new state w={0,1} from the distribution
                //      P(w) prop_to pi(w) prod_k=1^3 P_t_k(Y_k(i)|w)
                double[][] edgeRateMatrix = nwEvolver.calculateEdgeRateMatrix(currentNetwork, refNetwork, reaction, idx, evolutionModel, insRate, delRate, dependenceProbability);

                Iterator itNetwork = neighbourNetworks.iterator();
                Iterator itEvolTime = neighbourEvolTimes.iterator();
                double[] logProb = new double[]{0.0, 0.0};
                while (itNetwork.hasNext()) {
                    Byte[] theNetworkSeq = (Byte[]) itNetwork.next();
                    double evolTime = (Double) itEvolTime.next();

                    double[][] transitionProb = nwEvolver.calculateTransitionProbabilityMatrix(edgeRateMatrix, evolTime);
                    int edgeStatus = theNetworkSeq[idx];
                    logProb[0] += transitionProb[0][edgeStatus];
                    logProb[1] += transitionProb[1][edgeStatus];
                }
                // calculate the eq. prob
                double[] eqProb = nwEvolver.calculateEquilibriumProbability(edgeRateMatrix);
                // calculate sampling probabilities : pi(w) prod_k=1^3 P_t_k(Y_k(i)|w)
                double[] samplingProb = new double[2];
                for (int i = 0; i < samplingProb.length; i++) {
                    samplingProb[i] = eqProb[i] * Math.exp(logProb[i]);
                }
                // normalise them
                samplingProb = Utilities.normalise(samplingProb);
                // if (verbose)
                //     System.out.println(Double.toString(samplingProb[0]) + " " + Double.toString(samplingProb[1]));

                // sample new state for the edge
                newNetworkSeq[idx] = (Utilities.randsample(samplingProb) == 0 ? MetabolicNetwork.SEQ_ENTRY_ABSENT : MetabolicNetwork.SEQ_ENTRY_PRESENT);
            }

        } // end while it.hasNext()


        MetabolicNetwork newNetwork = refNetwork.getNetwork(newNetworkSeq);
        newNetwork.setID(mcmcOutput.getCurrentNetwork().getID());
        newNetwork.setName(mcmcOutput.getCurrentNetwork().getName());
        // set the initial state
        // set the current network
        mcmcOutput.setCurrentNetwork(newNetwork);

        // set the state for this iteration
        mcmcOutput.getNetworks().add(newNetworkSeq);
        if (phylonode.getParent() == null) { // save the rates at root only.
            mcmcOutput.setCurrentRates(new double[]{insRate, delRate});
//            mcmcOutput.getInsRates()[0] = insRate;
//            mcmcOutput.getDelRates()[0] = delRate;
        }

        if (debug) {
            System.out.println(newNetwork.getID() + ": " + Utilities.toString(newNetworkSeq, " "));
        }
    }

//    private double[] sampleRates(PhyloNode phylonode, MetabolicNetwork refNetwork, Byte[] iCoreEdges,
//            Byte[] iProhibEdges, int evolutionModel, double[] insParameters, double[] delParameters) {
//
//        int[][] changeSummary = new int[2][2]; // {[insertion-actual, insertion-allowed];[deletion-actual, deletion-allowed]}       
//
//        int networkCount = getEdgeChangeSummary(phylonode, refNetwork, iCoreEdges, iProhibEdges, evolutionModel, changeSummary);
//        double[][] avgChangeSummary = new double[changeSummary.length][changeSummary[0].length];
//        for (int i = 0; i < avgChangeSummary.length; i++) {
//            for (int j = 0; j < avgChangeSummary[0].length; j++) {
//                avgChangeSummary[i][j] = (double) changeSummary[i][j] / (double) networkCount;
//            }
//        }
//
//        // sample new rates
//        double[] newRates = new double[2];
//        Gamma gamma = new Gamma(1, 1);
//        int decimalPlaces = 6;
//        // insertion rate
//        newRates[0] = gamma.nextDouble(avgChangeSummary[0][0] + 1, avgChangeSummary[0][1]);
//        // deletion rate
//        newRates[1] = gamma.nextDouble(avgChangeSummary[1][0] + 1, avgChangeSummary[1][1]);
//
//        newRates[0] = Utilities.round((avgChangeSummary[0][1] == 0 ? 0.0 : newRates[0]), decimalPlaces);
//        newRates[1] = Utilities.round((avgChangeSummary[1][1] == 0 ? 0.0 : newRates[1]), decimalPlaces);
//
//        insParameters[0] = avgChangeSummary[0][0] + 1;
//        insParameters[1] = avgChangeSummary[0][1];
//        delParameters[0] = avgChangeSummary[1][0] + 1;
//        delParameters[1] = avgChangeSummary[1][1];
//
//        return newRates;
//    }
    private double[] sampleRates(PhyloNode phylonode, MetabolicNetwork refNetwork, Byte[] iCoreEdges,
            Byte[] iProhibEdges, int evolutionModel, double[] insParameters, double[] delParameters,
            int nInternalNodes, int alterableEdges) {

        int[][] changeSummary = new int[2][2]; // {[insertion-actual, insertion-allowed];[deletion-actual, deletion-allowed]}       
        getEdgeChangeSummary(phylonode, refNetwork, iCoreEdges, iProhibEdges, evolutionModel, changeSummary);

        // get the proportion of edges changed
        double[][] propChangeSummary = new double[changeSummary.length][changeSummary[0].length];
        for (int i = 0; i < propChangeSummary.length; i++) {
            for (int j = 0; j < propChangeSummary[0].length; j++) {
                propChangeSummary[i][j] = (double) changeSummary[i][j] / (double) (alterableEdges);
            }
        }


        // get the avg change summary by dividing by the number of network pairs (2*internal nodes)
        double[][] avgChangeSummary = new double[changeSummary.length][changeSummary[0].length];
        for (int i = 0; i < avgChangeSummary.length; i++) {
            for (int j = 0; j < avgChangeSummary[0].length; j++) {
                avgChangeSummary[i][j] = (double) propChangeSummary[i][j] / (double) (nInternalNodes * 2);
            }
        }

        // sample new rates
        double[] newRates = new double[2];
//        Gamma gamma = new Gamma(1, 1, RandomEngine.makeDefault());
        int decimalPlaces = 5;

        insParameters[0] = propChangeSummary[0][0] + 1;
        insParameters[1] = propChangeSummary[0][1];
        delParameters[0] = propChangeSummary[1][0] + 1;
        delParameters[1] = propChangeSummary[1][1];
        // insertion rate
        newRates[0] = gammaDist.nextDouble(insParameters[0], insParameters[1]);
        // deletion rate
        newRates[1] = gammaDist.nextDouble(delParameters[0], delParameters[1]);
        newRates[0] = Utilities.round((propChangeSummary[0][1] <= 0 ? 0.0 : newRates[0]), decimalPlaces);
        newRates[1] = Utilities.round((propChangeSummary[1][1] <= 0 ? 0.0 : newRates[1]), decimalPlaces);

//        Normal normal = new Normal(0, 1, RandomEngine.makeDefault());
//        newRates[0] = insRate + normal.nextDouble();
//        newRates[1] = delRate + normal.nextDouble();
//
//        newRates[0] = Utilities.round((newRates[0] < 0 ? 0.0 : newRates[0]), decimalPlaces);
//        newRates[1] = Utilities.round((newRates[1] < 0 ? 0.0 : newRates[1]), decimalPlaces);


        return newRates;
    }

    private int getEdgeChangeSummary(PhyloNode phylonode, MetabolicNetwork refNetwork, Byte[] iCoreEdges,
            Byte[] iProhibEdges, int evolutionModel, int[][] changeSummary) {

        int leftCount = 0, rightCount = 0;
        if (!phylonode.getLeftSon().isLeaf()) // run MCMC Step on left sub tree
        {
            leftCount = getEdgeChangeSummary(phylonode.getLeftSon(), refNetwork, iCoreEdges, iProhibEdges, evolutionModel, changeSummary);
        }
        if (!phylonode.getRightSon().isLeaf()) // run MCMC Step on right sub tree
        {
            rightCount = getEdgeChangeSummary(phylonode.getRightSon(), refNetwork, iCoreEdges, iProhibEdges, evolutionModel, changeSummary);
        }

        // get current, left & right sons' networks
        Byte[] currentNetworkSeq, leftNetworkSeq, rightNetworkSeq;
        currentNetworkSeq = ((MCMCOutput) phylonode.getData()).getCurrentNetwork().getReactionSequence();
        if (phylonode.getLeftSon().isLeaf()) {
            leftNetworkSeq = phylonode.getLeftSon().getMetabolicNetwork().getReactionSequence();
        } else {
            MCMCOutput leftMCMCOutput = (MCMCOutput) phylonode.getLeftSon().getData();
            leftNetworkSeq = leftMCMCOutput.getCurrentNetwork().getReactionSequence();
        }
        if (phylonode.getRightSon().isLeaf()) {
            rightNetworkSeq = phylonode.getRightSon().getMetabolicNetwork().getReactionSequence();
        } else {
            MCMCOutput rightMCMCOutput = (MCMCOutput) phylonode.getRightSon().getData();
            rightNetworkSeq = rightMCMCOutput.getCurrentNetwork().getReactionSequence();
        }

        int[][] changeSummaryLeft = getEdgeChangeSummary(currentNetworkSeq, leftNetworkSeq, iCoreEdges, iProhibEdges);
        int[][] changeSummaryRight = getEdgeChangeSummary(currentNetworkSeq, rightNetworkSeq, iCoreEdges, iProhibEdges);

        for (int i = 0; i < changeSummary.length; i++) {
            for (int j = 0; j < changeSummary[i].length; j++) {
                changeSummary[i][j] += changeSummaryLeft[i][j] + changeSummaryRight[i][j];
            }
        }

        return leftCount + rightCount + 2;
    }

    private int[][] getEdgeChangeSummary(Byte[] parentNetwork, Byte[] childNetwork,
            Byte[] iCoreEdges, Byte[] iProhibEdges) {

        int[][] changeSummary = new int[2][2]; // {[insertion-actual, insertion-allowed];[deletion-actual, deletion-allowed]}       
        for (int i = 0; i < parentNetwork.length; i++) {
            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT) || iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                continue;
            } else if (parentNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT)) { // case: insertion
                changeSummary[0][1]++;
                if (childNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                    changeSummary[0][0]++;
                }
            } else { // case: deletion
                changeSummary[1][1]++;
                if (childNetwork[i].equals(MetabolicNetwork.SEQ_ENTRY_ABSENT)) {
                    changeSummary[1][0]++;
                }
            }
        }

        return changeSummary;
    }

    private double treeProbability(PhyloNode phylonode, MetabolicNetwork refNetwork,
            Byte[] iCoreEdges, Byte[] iProhibEdges, double insRate, double delRate,
            int evolutionModel, int subNetworkSize, double dependenceProbability) {

        double leftTreeTP = 0.0, rightTreeTP = 0.0;
        if (!phylonode.getLeftSon().isLeaf()) // run MCMC Step on left sub tree
        {
            leftTreeTP = treeProbability(phylonode.getLeftSon(), refNetwork, iCoreEdges, iProhibEdges, insRate, delRate, evolutionModel, subNetworkSize, dependenceProbability);
        }
        if (!phylonode.getRightSon().isLeaf()) // run MCMC Step on right sub tree
        {
            rightTreeTP = treeProbability(phylonode.getRightSon(), refNetwork, iCoreEdges, iProhibEdges, insRate, delRate, evolutionModel, subNetworkSize, dependenceProbability);
        }

        // get a local copy of the network evolver
        NetworkEvolver nwEvolver = getNetworkEvolver();

        // get current, left & right sons' networks
        MetabolicNetwork currentNetwork, leftNetwork, rightNetwork;
        currentNetwork = ((MCMCOutput) phylonode.getData()).getCurrentNetwork();
        if (phylonode.getLeftSon().isLeaf()) {
            leftNetwork = phylonode.getLeftSon().getMetabolicNetwork();
        } else {
            MCMCOutput leftMCMCOutput = (MCMCOutput) phylonode.getLeftSon().getData();
            leftNetwork = leftMCMCOutput.getCurrentNetwork();
        }
        if (phylonode.getRightSon().isLeaf()) {
            rightNetwork = phylonode.getRightSon().getMetabolicNetwork();
        } else {
            MCMCOutput rightMCMCOutput = (MCMCOutput) phylonode.getRightSon().getData();
            rightNetwork = rightMCMCOutput.getCurrentNetwork();
        }
        // get evol times (branch lengths)
        double leftEvolTime, rightEvolTime;
        leftEvolTime = phylonode.getLeftBranchLength();
        rightEvolTime = phylonode.getRightBranchLength();

        // check if tp for these pairs are already calculated? if not recalculate them
        // left side
        String strNetworksRates = getNetworksRatesString(currentNetwork, leftNetwork, insRate, delRate, dependenceProbability);
        Double leftTP = (Double) distinctNetworkPairs.get(strNetworksRates);
        if (leftTP == null) {
            // Calculate the probability
            leftTP = nwEvolver.approximateTransitionProbability(currentNetwork, leftNetwork,
                    refNetwork, iCoreEdges, iProhibEdges, evolutionModel, insRate, delRate,
                    leftEvolTime, subNetworkSize, dependenceProbability);
            // add it to the list of distinct network pairs
            distinctNetworkPairs.put(strNetworksRates, leftTP);
        }
        // right side
        strNetworksRates = getNetworksRatesString(currentNetwork, rightNetwork, insRate, delRate, dependenceProbability);
        Double rightTP = (Double) distinctNetworkPairs.get(strNetworksRates);
        if (rightTP == null) {
            // Calculate the probability
            rightTP = nwEvolver.approximateTransitionProbability(currentNetwork, rightNetwork,
                    refNetwork, iCoreEdges, iProhibEdges, evolutionModel, insRate, delRate,
                    rightEvolTime, subNetworkSize, dependenceProbability);
            // add it to the list of distinct network pairs
            distinctNetworkPairs.put(strNetworksRates, rightTP);
        }

        return leftTreeTP + rightTreeTP + Math.log(leftTP) + Math.log(rightTP);
    }

//    public Double[] estimateRates(PhyloTree phylotree, MetabolicNetwork refNetwork, ArrayList coreEdges,
//            ArrayList prohibEdges, int nIter, int evolutionModel, int nBurning, int subNetworkSize) {
//        // estimate rates using the fixed networks at internal nodes. The function assumes that the 
//        // internal networks to be used to estimate rates are set at each internal node.
//
//        // get indices of core and prohibited edges
//        Iterator itCore = coreEdges.iterator();
//        Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
//        for (int i = 0; i < iCoreEdges.length; i++)
//            iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//        while (itCore.hasNext()) {
//            DirectedReaction reaction = (DirectedReaction) itCore.next();
//            iCoreEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
//        }
//        Iterator itProhib = prohibEdges.iterator();
//        Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
//        for (int i = 0; i < iProhibEdges.length; i++)
//            iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
//        while (itProhib.hasNext()) {
//            DirectedReaction reaction = (DirectedReaction) itProhib.next();
//            iProhibEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
//        }
//
//        Double[] mlEstimate = new Double[3];
//        for (int i = 0; i < mlEstimate.length; i++)
//            mlEstimate[i] = new Double(0.0);
//
//        double insRate = Math.random();
//        double delRate = Math.random();
//
//        // Calculate the tree probability
//        BigDecimal treeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(), refNetwork,
//                iCoreEdges, iProhibEdges, insRate, delRate, evolutionModel, subNetworkSize)));
//
//        int nInternalNodes = phylotree.getLeaves().size() - 1;
//        for (int iter = 0; iter < nIter; iter++) {
//
//            // set the parameters
//            double[] insParameters = new double[2];
//            double[] delParameters = new double[2];
//            // sample new inserstion / deletion rates
//            double[] newRates = sampleRates(phylotree.getRoot(), refNetwork, iCoreEdges, iProhibEdges,
//                    evolutionModel, insParameters, delParameters, nInternalNodes);
//            double newInsRate = newRates[0];
//            double newDelRate = newRates[1];
//
//            BigDecimal newTreeProb;
//            String strTree = getTreeRatesString(phylotree, newInsRate, newDelRate);
//            newTreeProb = (BigDecimal) distinctTrees.get(strTree);
//            if (newTreeProb == null) {
//                newTreeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
//                        refNetwork, iCoreEdges, iProhibEdges, newInsRate, newDelRate,
//                        evolutionModel, subNetworkSize)));
//
//                distinctTrees.put(strTree, newTreeProb);
//                
//            }
////            BigDecimal newTreeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
////                    refNetwork, iCoreEdges, iProhibEdges, newInsRate, newDelRate,
////                    evolutionModel, subNetworkSize)));
//
//
//            // calculate the acceptance probability
//            BigDecimal alphaRates = calculateRateAcceptanceProbability(insParameters, delParameters,
//                    insRate, delRate, newInsRate, newDelRate, treeProb, newTreeProb);
//
//            // accept the new state with probability accept
//            double acceptRates = Math.min(1.0, alphaRates.doubleValue());
//            double uAcceptRates = Math.random();
//
//            if (newInsRate == 0 || newDelRate == 0) {
//            // do nothing
//            } else if (uAcceptRates <= acceptRates) {
//                // return the new state
//                insRate = newInsRate;
//                delRate = newDelRate;
//
//                // also copy the new path probability
//                treeProb = newTreeProb;
//            }
//
//            if (iter > nBurning && mlEstimate[2] < treeProb.doubleValue()) {
//                mlEstimate[0] = insRate;
//                mlEstimate[1] = delRate;
//                mlEstimate[2] = treeProb.doubleValue();
//            }
//
//        } // end for iter
//
//        return mlEstimate;
//    }
    public void estimateRates(PhyloTree phylotree, MetabolicNetwork refNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int nIter, int updateInterval, int evolutionModel, int nBurning,
            int subNetworkSize, double dependenceProbability) {

        // get indices of core and prohibited edges
        Iterator itCore = coreEdges.iterator();
        Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iCoreEdges.length; i++) {
            iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itCore.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itCore.next();
            iCoreEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }
        Iterator itProhib = prohibEdges.iterator();
        Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iProhibEdges.length; i++) {
            iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itProhib.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itProhib.next();
            iProhibEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }

        Double[][] rates = new Double[nIter - nBurning + 1][2];

        double insRate = Math.random();
        double delRate = Math.random();

        // Calculate the path probability
        BigDecimal treeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(), refNetwork,
                iCoreEdges, iProhibEdges, insRate, delRate, evolutionModel, subNetworkSize, dependenceProbability)));

        if (updateInterval > 0 && 0 % updateInterval == 0) {
            System.out.print("Iteration: 0");
            System.out.println("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate));
        }

        setNetworksAtInternalNodes(phylotree.getRoot(), refNetwork, 0);
        int alterableEdges = refNetwork.getReactionSequence().length - Utilities.sum(iCoreEdges) - Utilities.sum(iProhibEdges);
        int nInternalNodes = phylotree.getLeaves().size() - 1;
        for (int iter = 1; iter <= nIter; iter++) {

            // set the parameters
            double[] insParameters = new double[2];
            double[] delParameters = new double[2];
            // sample new inserstion / deletion rates
            double[] newRates = sampleRates(phylotree.getRoot(), refNetwork, iCoreEdges, iProhibEdges,
                    evolutionModel, insParameters, delParameters, nInternalNodes, alterableEdges);
            double newInsRate = newRates[0];
            double newDelRate = newRates[1];

            BigDecimal newTreeProb;
            String strTree = getTreeRatesString(phylotree, newInsRate, newDelRate, dependenceProbability);
            newTreeProb = (BigDecimal) distinctTrees.get(strTree);
            if (newTreeProb == null) {
                newTreeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
                        refNetwork, iCoreEdges, iProhibEdges, newInsRate, newDelRate,
                        evolutionModel, subNetworkSize, dependenceProbability)));

                distinctTrees.put(strTree, newTreeProb);

            }

//            BigDecimal newTreeProb = new BigDecimal(Math.exp(treeProbability(phylotree.getRoot(),
//                    refNetwork, iCoreEdges, iProhibEdges, newInsRate, newDelRate,
//                    evolutionModel, subNetworkSize)));

            // calculate the acceptance probability
            BigDecimal alphaRates = calculateRateAcceptanceProbability(insParameters, delParameters,
                    insParameters, delParameters, insRate, delRate, newInsRate, newDelRate, treeProb, newTreeProb);

            // accept the new state with probability accept
            double acceptRates = Math.min(1.0, alphaRates.doubleValue());
            double uAcceptRates = Math.random();

            if (newInsRate == 0 || newDelRate == 0) {
                // do nothing
            } else if (uAcceptRates <= acceptRates) {
                // return the new state
                insRate = newInsRate;
                delRate = newDelRate;

                // also copy the new path probability
                treeProb = newTreeProb;
            }

            if (iter > nBurning) {
                // save the rates for the current iteration
                rates[iter - nBurning - 1][0] = insRate;
                rates[iter - nBurning - 1][1] = delRate;
            }

            if (updateInterval > 0 && iter % updateInterval == 0) {
                System.out.print("Iteration: " + Integer.toString(iter));
                System.out.println("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate));
            }

            setNetworksAtInternalNodes(phylotree.getRoot(), refNetwork, iter);

        } // end for iter = 1 to nIter
    }

    public void setNetworksAtInternalNodes(PhyloNode phylonode, MetabolicNetwork refNetwork, int index) {

        if (!phylonode.getLeftSon().isLeaf()) {
            setNetworksAtInternalNodes(phylonode.getLeftSon(), refNetwork, index);
        }
        if (!phylonode.getRightSon().isLeaf()) {
            setNetworksAtInternalNodes(phylonode.getRightSon(), refNetwork, index);
        }

        MCMCOutput mcmcOutput = (MCMCOutput) phylonode.getData();
        mcmcOutput.setCurrentNetwork(refNetwork.getNetwork((Byte[]) mcmcOutput.getNetworks().get(index)));
    }

    public int[] proposeNewPath(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int nDifferences) {

        // Get the probabilities of selecting the next move
        BigDecimal[] moveProb = pathLengthProposalProbability(path.length, nDifferences);
        // select a move
        int selectedMove = Utilities.randsample(moveProb);

        int[] newPath;
        switch (selectedMove) {
            case MOVE_ADD_EVENTS:
                newPath = addEvents(path, startNetwork, coreEdges, prohibEdges, null);
                break;
            case MOVE_SWAP_EVENTS:
                newPath = swapEvents(path, startNetwork, coreEdges, prohibEdges, null);
                break;
            case MOVE_DEL_EVENTS:
                newPath = deleteEvents(path, startNetwork, coreEdges, prohibEdges, null);
                break;
            case MOVE_PERMUTE_EVENTS:
                newPath = permuteEvents(path, startNetwork, coreEdges, prohibEdges);
                break;
            default:
                newPath = new int[]{};
        }

        return newPath;
    }

    private int[] addEvents(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int[] parameters) {

        return addEvents(path, startNetwork, coreEdges, prohibEdges, parameters, -1, -1, -1);
    }

    public int[] addEvents(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int[] parameters, int position1, int position2, int edge) {

        ArrayList lstEdges = startNetwork.getDirectedReactions();
        // total number of edges
        int edgeCount = lstEdges.size();

        // initialise the variable to store the new path
        int[] newPath = new int[path.length + 2];

        // initialise a variable to keep track to edges already selected
        int sIdx = 0;
        int[] selectedEdges = new int[edgeCount];

        if (edge <= 0) { // no edge was specified. we need to select an edge.
            Iterator it = coreEdges.iterator();
            while (it.hasNext()) {
                // get the next edge
                DirectedReaction theEdge = (DirectedReaction) it.next();
                // mark this edge as already selected
                selectedEdges[lstEdges.indexOf(theEdge)] = 1;
                sIdx++;
            } // while it.hasNext()
            it = prohibEdges.iterator();
            while (it.hasNext()) {
                // get the next edge
                DirectedReaction theEdge = (DirectedReaction) it.next();
                // mark this edge as already selected
                selectedEdges[lstEdges.indexOf(theEdge)] = 1;
                sIdx++;
            } // while it.hasNext()
        }

        // Loop until a valid path is found or all possible edges have been 
        // check for the possible possitions. 
        while (true) {

            int eventEdge;
            if (edge <= 0) {
                // get the edges that can be inserted
                int idx = 0;
                int[] insertableEdges = new int[edgeCount - sIdx];
                for (int i = 0; i < selectedEdges.length; i++) {
                    if (selectedEdges[i] != 1) {
                        insertableEdges[idx++] = i + 1;
                    }
                }

                // select an edge
                eventEdge = insertableEdges[(int) Math.floor(Math.random() * insertableEdges.length)];
                // mark the edge as selected
                selectedEdges[eventEdge - 1] = 1;
                sIdx++;

            } else {
                eventEdge = edge;
            }
            // Initialize the variable to keep track of the visited positions.
            // Events are appended after the selected position
            // Edges can be inserted on the either end points as well so we need
            // keep track of those positinos as well.
            int[] visitedPos = new int[path.length + 1];

            // Loop until a valid path is found or all possible positions 
            // for the edge has been checked. 
            while (true) {

                int selectedPosition1, selectedPosition2;
                if (position1 < 0) {
                    // select the positions to add the events
                    selectedPosition1 = (int) Math.floor(Math.random() * (path.length + 1)) - 1;
                    selectedPosition2 = (int) Math.floor(Math.random() * (path.length + 1)) - 1;

                    // mark the positions just visited
                    visitedPos[selectedPosition1 + 1] = 1;
                    visitedPos[selectedPosition2 + 1] = 1;
                } else {
                    selectedPosition1 = position1;
                    selectedPosition2 = position2;
                }

                // Append the first event at position specified by selectedPosition1
                if (selectedPosition1 > selectedPosition2) {
                    int temp = selectedPosition1;
                    selectedPosition1 = selectedPosition2;
                    selectedPosition2 = temp;
                } else if (selectedPosition2 > selectedPosition1) // We add one to position 2 as the positions have been shifted
                // by 1 due to insertion of the event at position 1
                {
                    selectedPosition2 = selectedPosition2 + 1;
                }

                // create the new path
                int pIdx = 0;
                for (int i = 0; i < selectedPosition1; i++) {
                    newPath[pIdx++] = path[i];
                }
                newPath[pIdx++] = eventEdge;
                selectedPosition1 = Math.max(0, selectedPosition1);
                selectedPosition2 = Math.min(path.length, selectedPosition2);
                for (int i = selectedPosition1; i < selectedPosition2; i++) {
                    newPath[pIdx++] = path[i];
                }
                newPath[pIdx++] = eventEdge;
                selectedPosition2 = Math.max(0, selectedPosition2);
                for (int i = selectedPosition2; i < path.length; i++) {
                    newPath[pIdx++] = path[i];
                }

                if (isValidPath(newPath, startNetwork, coreEdges, prohibEdges)) {
                    if (parameters != null) {
                        // set the parameters
                        parameters[0] = selectedPosition1;
                        parameters[1] = selectedPosition2;
                        parameters[2] = eventEdge;
                    }
                    return newPath;

                } else {
                    int prod = 0;
                    if (position1 < 0) // positions were not specified. check if we have seen all positions
                    {
                        for (int i = 0; i < visitedPos.length; i++) {
                            prod *= visitedPos[i];
                        }
                    } else {
                        prod = 1;
                    }
                    if (prod == 1) {
                        break;
                    }
                }
            } // end while true - positions

            int prod = 0;
            if (edge <= 0) // no edge was specified. check if we have seen all edges
            {
                for (int i = 0; i < selectedEdges.length; i++) {
                    prod *= selectedEdges[i];
                }
            } else {
                prod = 1;
            }
            if (prod == 1) {
                break;
            }
        } // while true - edge

        // we come here if no valid path was found. 
        return null;
    }

    public int[] swapEvents(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int[] parameters) {

        return swapEvents(path, startNetwork, coreEdges, prohibEdges, parameters, -1, -1);
    }

    public int[] swapEvents(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int[] parameters, int position1, int position2) {

        // initialise the variable to store the new path
        int[] newPath = new int[path.length];

        // Initialize the variable to keep track of the visited positions.
        int[] visitedPos = new int[path.length];

        while (true) {
            int pos1, pos2, eventEdge, eventEdge2;

            if (position1 < 0) {
                // select a random position
                pos1 = (int) Math.floor(Math.random() * path.length);
                // get the event edge at this position
                eventEdge = path[pos1];
                // mark the position
                visitedPos[pos1] = 1;

                // select second position. should have a different edge 
                int[] secondPos = new int[path.length];
                int sIdx = 0;
                for (int i = 0; i < path.length; i++) {
                    if (path[i] != eventEdge) {
                        secondPos[sIdx++] = i;
                    }
                }

                if (sIdx == 0) {
                    break;
                }

                pos2 = secondPos[(int) Math.floor(Math.random() * sIdx)];
                // get the event edge at this position
                eventEdge2 = path[pos2];
                // mark the position
                visitedPos[pos2] = 1;
            } else {
                pos1 = position1;
                pos2 = position2;
                // get the event edge at these position
                eventEdge = path[pos1];
                eventEdge2 = path[pos2];

                if (path[pos1] == path[pos2]) {
                    return null;
                }
            }
            if (pos1 > pos2) {
                int temp = pos1;
                pos1 = pos2;
                pos2 = temp;

                // also swap the events
                temp = eventEdge;
                eventEdge = eventEdge2;
                eventEdge2 = temp;
            }

            // create the new path
            int pIdx = 0;
            for (int i = 0; i < pos1; i++) {
                newPath[pIdx++] = path[i];
            }
            newPath[pIdx++] = eventEdge2;
            for (int i = pos1 + 1; i < pos2; i++) {
                newPath[pIdx++] = path[i];
            }
            newPath[pIdx++] = eventEdge;
            for (int i = pos2 + 1; i < path.length; i++) {
                newPath[pIdx++] = path[i];
            }

            if (isValidPath(newPath, startNetwork, coreEdges, prohibEdges)) {
                if (parameters != null) {
                    // set the parameters
                    parameters[0] = pos1;
                    parameters[1] = pos2;
                    parameters[2] = eventEdge;
                }
                return newPath;
            } else {
                int prod = 0;
                for (int i = 0; i < visitedPos.length; i++) {
                    prod *= visitedPos[i];
                }

                if (prod == 1) {
                    break;
                }
            }

        } // end while true

        // we come here if no valid path was found. 
        return null;
    }

    private int[] deleteEvents(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int[] parameters) {

        return deleteEvents(path, startNetwork, coreEdges, prohibEdges, parameters, -1, -1, -1);
    }

    public int[] deleteEvents(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int[] parameters, int position1, int position2, int edge) {

        ArrayList lstEdges = startNetwork.getDirectedReactions();
        // total number of edges
        int edgeCount = lstEdges.size();

        // initialise the variable to store the new path
        int[] newPath = new int[path.length - 2];

        // get the edge frequencies
        int[] edgeFreq = new int[edgeCount];
        for (int i = 0; i < path.length; i++) {
            edgeFreq[path[i] - 1]++;
        }

        // NOTE: We dont need to remove core edges as they dont occur in the path

        // Initialize the variable to keep track of the visited positions. We
        // add all the positions which correnspond to the events occuring once
        // only as visited
        int[] visitedPos = new int[path.length];
        int[] events = new int[path.length]; // variable to track of the positions with egdes occuring more than once
        int eIdx = 0;
        // initialise a variable to count the number of edges that are present more than once
        for (int i = 0; i < path.length; i++) {
            if (edgeFreq[path[i] - 1] <= 1) {
                visitedPos[i] = 1;
            } else if (edge > 0 && path[i] != edge) {
                visitedPos[i] = 1;
            } else {
                events[eIdx++] = i;
            }
        }

        // get the deletable edges
        int[] deletableEvents = new int[eIdx];
        for (int i = 0; i < eIdx; i++) {
            deletableEvents[i] = events[i];
        }

        if (deletableEvents.length == 0) {
            System.out.println(Utilities.toString(startNetwork.getReactionSequence()));
            System.out.println(Utilities.toString(path));
            return new int[]{};
        }

        while (true) {

            int pos1, pos2, eventEdge;
            if (position1 < 0) {
                // Get the positions corresponding to edges which occur
                // more than onces
                pos1 = deletableEvents[(int) Math.floor(Math.random() * deletableEvents.length)];
                // get the edge for the event
                eventEdge = path[pos1];
                // select the second position
                int[] secondPos = new int[deletableEvents.length];
                int sIdx = 0;
                for (int i = 0; i < path.length; i++) {
                    if (i != pos1 && path[i] == eventEdge) {
                        secondPos[sIdx++] = i;
                    }
                }
                pos2 = secondPos[(int) Math.floor(Math.random() * sIdx)];
            } else {
                if (edge > 0 && path[position1] != edge) {
                    return null;
                }
                pos1 = position1;
                pos2 = position2;
                // get the edge for the event
                eventEdge = path[pos1];

            }

            if (pos1 > pos2) {
                int temp = pos1;
                pos1 = pos2;
                pos2 = temp;
            }

            // create the new path
            int pIdx = 0;
            for (int i = 0; i < pos1; i++) {
                newPath[pIdx++] = path[i];
            }
            for (int i = pos1 + 1; i < pos2; i++) {
                newPath[pIdx++] = path[i];
            }
            for (int i = pos2 + 1; i < path.length; i++) {
                newPath[pIdx++] = path[i];
            }

            if (isValidPath(newPath, startNetwork, coreEdges, prohibEdges)) {
                if (parameters != null) {
                    // set the parameters
                    parameters[0] = pos1;
                    parameters[1] = pos2;
                    parameters[2] = eventEdge;
                }
                return newPath;
            } else {
                int prod = 0;
                if (position1 < 0) // positions were not specified. check if we have seen all positions
                {
                    for (int i = 0; i < visitedPos.length; i++) {
                        prod *= visitedPos[i];
                    }
                } else {
                    prod = 1;
                }

                if (prod == 1) {// no 0 present i.e. all positions have been visited
                    break;
                }
            }
        } // end while true

        // we come here if no valid path was found. 
        return null;
    }

    public int[] permuteEvents(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges) {

        // check if all events are the same in the path
        int firstEvent = path[0];
        boolean allSame = true;
        for (int i = 1; i < path.length; i++) {
            if (path[i] != firstEvent) {
                allSame = false;
                break;
            }
        }
        if (allSame) // permutation not possible
        {
            return path;
        }

        // initialise the variable to store the new path
        int[] newPath = new int[path.length];

        // Initialize the variable to keep track of the visited permutations.
        HashSet visitedPermutations = new HashSet();
        // add the given path to it
        visitedPermutations.add(Utilities.toString(path));

        while (true) {
            newPath = Utilities.permute(path);

            if (visitedPermutations.add(Utilities.toString(newPath))) {
                if (isValidPath(newPath, startNetwork, coreEdges, prohibEdges)) {
                    return newPath;
                } else {
                    if (visitedPermutations.size() == Integer.MAX_VALUE) {
                        break;
                    }
                }
            }
        } // end while true

        // we come here if no valid path was found. 
        return null;
    }

    private boolean isValidPath(int[] newPath, MetabolicNetwork startNetwork, ArrayList coreEdges, ArrayList prohibEdges) {
        return true;
    }

    private BigDecimal[] pathLengthProposalProbability(int pathLen, int nDifferences) {
        /** 
         *  Probabilities of proposing a different path length
         */

        // Get the 128 Decimal Math Context
        MathContext mc = MathContext.DECIMAL128;

        BigDecimal[] prob = new BigDecimal[NO_OF_ALLOWED_MOVES];
        int nAllowedMoves = 0; // no of allowed moves for the current case

        if (pathLen <= 1) { // Only addition is allowed if there is only one element in the path 
            prob[MOVE_ADD_EVENTS] = new BigDecimal(1.0, mc);
            nAllowedMoves = 1;
        } else if (pathLen == nDifferences) {
            // Only addition and swapping is allowed if the path path length equals the number of differences (i.e. length of most parsimonious path)
            prob[MOVE_ADD_EVENTS] = new BigDecimal(1.0, mc);
//            prob[MOVE_SWAP_EVENTS] = new BigDecimal(1.0, mc);
            prob[MOVE_PERMUTE_EVENTS] = new BigDecimal(1.0, mc);
            nAllowedMoves = 2;
        } else if (pathLen >= Math.max(nDifferences + 10, nDifferences * 3)) { // dont explore too long paths
            // Only deletion, swapping and permutation is allowed 
//            prob[MOVE_SWAP_EVENTS] = new BigDecimal(1.0, mc);
            prob[MOVE_DEL_EVENTS] = new BigDecimal(1.0, mc);
            prob[MOVE_PERMUTE_EVENTS] = new BigDecimal(1.0, mc);
            nAllowedMoves = 2;
        } else if (pathLen > nDifferences) { // otherwise all three operations are allowed
            prob[MOVE_ADD_EVENTS] = new BigDecimal(1.0, mc);
//            prob[MOVE_SWAP_EVENTS] = new BigDecimal(1.0, mc);
            prob[MOVE_DEL_EVENTS] = new BigDecimal(1.0, mc);
            prob[MOVE_PERMUTE_EVENTS] = new BigDecimal(1.0, mc);
            nAllowedMoves = 3;
        }

        // Each operation has equal probability
        for (int i = 0; i < prob.length; i++) {
            if (prob[i] != null) {
                prob[i] = prob[i].divide(new BigDecimal(nAllowedMoves, mc), mc);
            } else {
                prob[i] = new BigDecimal(0.0);
            }
        }

        return prob;
    }

    public BigDecimal pathProbability(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, double insRate, double delRate, double evolTime, int evolutionMode,
            double dependenceProbability) {

        // Get the 128 Decimal Math Context
        MathContext mc = MathContext.DECIMAL128;

        // Get the vectors of rates and total rates 
        double[] rates = new double[path.length];
        double[] totRates = new double[path.length + 1];
        getChainParameters(path, startNetwork, coreEdges, prohibEdges, evolutionMode, insRate, delRate, rates, totRates, dependenceProbability);

        // initialize the log probability
        BigDecimal prob = new BigDecimal(1.0);
        if (evolTime == -1) // Calculate jump probability
        {
            for (int i = 0; i < rates.length; i++) {
                BigDecimal tp = new BigDecimal(rates[i], mc).divide(new BigDecimal(totRates[i], mc), mc);
                prob = prob.multiply(tp, mc);
            }
        } else {
            prob = trajectoryLikelihood(rates, totRates, evolTime);
            if (prob.compareTo(new BigDecimal(0)) == -1) {
                prob = new BigDecimal(0);
            } else if (prob.compareTo(new BigDecimal(1)) == 1) {
                prob = new BigDecimal(0);
            }

        }
        return prob;
    }

    public void getChainParameters(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int evolutionMode, double insRate, double delRate,
            double[] rates, double[] totRates, double dependenceProbability) {

        // copy the starting network
        MetabolicNetwork theNetwork = startNetwork.clone();

        for (int i = 0; i < path.length; i++) {
            // get the edge corresponding to ith event
            DirectedReaction reaction = (DirectedReaction) theNetwork.getDirectedReactions().get(path[i] - 1);

            // Get the event rate & total rate
            double[] rate = getJumpParameters(theNetwork, reaction, coreEdges, prohibEdges, evolutionMode, insRate, delRate, dependenceProbability);

            rates[i] = rate[0];
            totRates[i] = rate[1];

            // Update the network
            theNetwork = theNetwork.updateNetwork(reaction);
        } // end for each event in the path

        // Get the total rate from the final network
        double[] rate = getJumpParameters(theNetwork, null, coreEdges, prohibEdges, evolutionMode, insRate, delRate, dependenceProbability);
        totRates[totRates.length - 1] = rate[1];

    }

    private double[] getJumpParameters(MetabolicNetwork theNetwork, DirectedReaction eventEdge,
            ArrayList coreEdges, ArrayList prohibEdges, int evolutionMode, double insRate, double delRate,
            double dependenceProbability) {

        double[] parameters = new double[2]; // rate, totRate

        // calculate edge rates
        double[] edgeRates = getNetworkEvolver().calculateHyperedgeChangeRates(theNetwork, coreEdges, prohibEdges, evolutionMode, insRate, delRate, dependenceProbability);

        // extract the rate for the event edge
        double eventEdgeRate;
        if (eventEdge != null) {
            int edgeIdx = theNetwork.getDirectedReactions().indexOf(eventEdge);
            eventEdgeRate = edgeRates[edgeIdx];
        } else {
            eventEdgeRate = 0.0;
        }

        // set the parameter values to be returned
        parameters[0] = eventEdgeRate; // rate
        parameters[1] = Utilities.sum(edgeRates); // total Rate

        return parameters;
    }

    public void getRateParameters(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, int[] nEvents, int[] nTotalEvents) {


        int nIns = 0; // number of insertion events in the path
        int nDel = 0; // number of deletion events in the path
        int[] insCounts = new int[path.length]; // number of insertable edges at each state of the path
        int[] delCounts = new int[path.length]; // number of deletable edges at each state of the path
        // copy the starting network
        MetabolicNetwork theNetwork = startNetwork.clone();

        for (int i = 0; i < path.length; i++) {
            // Get the edges that can be inserted into the given network (inactive edges - prohibited edges)
            ArrayList insertableEdges = theNetwork.getInactiveDirectedReactions();
            insertableEdges.removeAll(prohibEdges);
            insCounts[i] = insertableEdges.size();

            // Get the edges that can be deleted from the given network (active edges - core edges)
            ArrayList deletableEdges = theNetwork.getActiveDirectedReactions();
            deletableEdges.removeAll(coreEdges);
            delCounts[i] = deletableEdges.size();

            // increment the appropriate event count
            if (theNetwork.getReactionSequence()[path[i] - 1].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                nDel++;
            } else {
                nIns++;
            }

            // get the edge corresponding to ith event
            DirectedReaction reaction = (DirectedReaction) theNetwork.getDirectedReactions().get(path[i] - 1);
            // Update the network
            theNetwork = theNetwork.updateNetwork(reaction);
        } // end for each event in the path

        // set the return values
        nEvents[0] = nIns;
        nEvents[1] = nDel;
        nTotalEvents[0] = Utilities.sum(insCounts);
        nTotalEvents[1] = Utilities.sum(delCounts);

    }

    public BigDecimal trajectoryLikelihood(double[] rates, double[] totExitRates, double evolTime) {

        // sorting the totalExitRates
        Arrays.sort(totExitRates);

        MathContext mc = new MathContext(PRECISION);

        BigDecimal[] zeta = new BigDecimal[totExitRates.length];
        for (int i = 0; i < totExitRates.length; i++) {
            zeta[i] = BigDecimal.valueOf(totExitRates[i]);
        }

        // Allocate space
        BigDecimal[][] coefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[][] tempCoefficients = new BigDecimal[zeta.length][zeta.length];
        BigDecimal[] ksi = new BigDecimal[zeta.length];
        int[] degeneracy = new int[zeta.length];

        int M = 0;
        ksi[0] = zeta[0];
        degeneracy[0] = 0;
        coefficients[0][0] = new BigDecimal(1.0);

        for (int i = 1; i < zeta.length; i++) {

            for (int n = 0; n < M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    tempCoefficients[n][k] = new BigDecimal(0.0);
                    BigDecimal factor = new BigDecimal(1.0).divide(ksi[n].subtract(zeta[i]), mc);
                    for (int j = k; j <= degeneracy[n]; j++) {
                        tempCoefficients[n][k] = tempCoefficients[n][k].add(coefficients[n][j].multiply(factor));
                        factor = factor.multiply(new BigDecimal(j + 1.0).divide(ksi[n].subtract(zeta[i]), mc));
                    }
                    tempCoefficients[n][k] = tempCoefficients[n][k].negate();
                }
            }

            if (zeta[i].compareTo(ksi[M]) != 0) {

                for (int k = 0; k <= degeneracy[M]; k++) {
                    tempCoefficients[M][k] = new BigDecimal(0.0);
                    BigDecimal factor = new BigDecimal(1.0).divide(ksi[M].subtract(zeta[i]), mc);
                    for (int j = k; j <= degeneracy[M]; j++) {
                        tempCoefficients[M][k] = tempCoefficients[M][k].add(coefficients[M][j].multiply(factor));
                        factor = factor.multiply(new BigDecimal(j + 1.0).divide(ksi[M].subtract(zeta[i]), mc));
                    }
                    tempCoefficients[M][k] = tempCoefficients[M][k].negate();
                }

                ksi[++M] = zeta[i];
                degeneracy[M] = 0;

                tempCoefficients[M][0] = new BigDecimal(0.0);
                for (int n = 0; n < M; n++) {
                    if (ksi[n].compareTo(zeta[i]) != 0) {
                        BigDecimal factor = new BigDecimal(1.0).divide(ksi[n].subtract(zeta[i]), mc);
                        for (int k = 0; k <= degeneracy[n]; k++) {
                            tempCoefficients[M][0] = tempCoefficients[M][0].add(coefficients[n][k].multiply(factor));
                            factor = factor.multiply(new BigDecimal(k + 1.0).divide(ksi[n].subtract(zeta[i]), mc));
                        }
                    }
                }
            } else { // not new
                degeneracy[M]++;
                for (int k = 1; k <= degeneracy[M]; k++) {
                    tempCoefficients[M][k] = coefficients[M][k - 1].divide(new BigDecimal(k), mc);
                }

                tempCoefficients[M][0] = new BigDecimal(0.0);
                for (int n = 0; n < M; n++) {
                    BigDecimal factor = new BigDecimal(1.0).divide(ksi[n].subtract(zeta[i]), mc);
                    for (int k = 0; k <= degeneracy[n]; k++) {
                        tempCoefficients[M][0] = tempCoefficients[M][0].add(coefficients[n][k].multiply(factor));
                        factor = factor.multiply(new BigDecimal(k + 1.0).divide(ksi[n].subtract(zeta[i]), mc));
                    }
                }
            }

            for (int n = 0; n <= M; n++) {
                for (int k = 0; k <= degeneracy[n]; k++) {
                    coefficients[n][k] = tempCoefficients[n][k];
                //report for debugging reasons:
//                    System.out.print(coefficients[n][k].toString() + ",\t");
                }
//                System.out.println();
            }
        } //for (int i = 1; i < zeta.length; i++)

        BigDecimal jump = new BigDecimal(1.0);
        for (int i = 0; i < rates.length; i++) {
            jump = jump.multiply(BigDecimal.valueOf(rates[i]));
        }

        BigDecimal temp = new BigDecimal(0.0);
        for (int n = 0; n <= M; n++) {
            BigDecimal cT = new BigDecimal(0.0);
            for (int k = 0; k <= degeneracy[n]; k++) {
                cT = cT.add(coefficients[n][k].multiply(new BigDecimal(evolTime).pow(k)));
            }
            BigDecimal factor = Utilities.exp(ksi[n].multiply(new BigDecimal(evolTime)).negate(), this.natural_e, this.epsilon, PRECISION, this.factorial);
            temp = temp.add(factor.multiply(cT));

        }
        return jump.multiply(temp, mc);

    }

    public BigDecimal calculatePathAcceptanceProbability(BigDecimal pathProb, BigDecimal newPathProb) {

        // Acceptance probability
        BigDecimal numerator = newPathProb;
        BigDecimal denominator = pathProb;

        BigDecimal accept;
        if (denominator.compareTo(BigDecimal.ZERO) == 0) {
            accept = new BigDecimal(0.0);
        } else {
            accept = numerator.divide(denominator, MathContext.DECIMAL128);
        }

        return accept;
    }

    public BigDecimal calculatePathAcceptanceProbability(BigDecimal pathProb, BigDecimal newPathProb, BigDecimal startNetworkProb, BigDecimal newStartNetworkProb) {

        // Acceptance probability
        BigDecimal numerator = newPathProb.multiply(newStartNetworkProb, MathContext.DECIMAL128);
        BigDecimal denominator = pathProb.multiply(startNetworkProb, MathContext.DECIMAL128);

        BigDecimal accept;
        if (denominator.compareTo(BigDecimal.ZERO) == 0) {
            accept = new BigDecimal(0.0);
        } else {
            accept = numerator.divide(denominator, MathContext.DECIMAL128);
        }

        return accept;
    }

    public BigDecimal calculatePathAcceptanceProbability(int[] path, int[] newPath, int edgeCount,
            int nDifferences, BigDecimal pathProb, BigDecimal newPathProb) {


        int pathLen = path.length;
        int newPathLen = newPath.length;

        BigDecimal accept;
        if (newPathLen == pathLen) // Acceptance probability is simply the ratio of likelihoods (path probabilities)
        {
            if (pathProb.compareTo(BigDecimal.ZERO) == 0) {
                accept = BigDecimal.ONE;
            } else {
                accept = newPathProb.divide(pathProb, MathContext.DECIMAL128);
            }
        } else {

            // Probability for selecting new path length
            BigDecimal pathLengthProb = getPathLengthProbability(pathLen, newPathLen, nDifferences);
            BigDecimal newPathLengthProb = getPathLengthProbability(newPathLen, pathLen, nDifferences); //reverse move

            // Density of proposed path
            BigDecimal propDensity = getPathProposalDensity(path, newPath, edgeCount);
            BigDecimal newPropDensity = getPathProposalDensity(newPath, path, edgeCount); //reverse move

            // Acceptance probability
            BigDecimal numerator = newPathProb.multiply(newPathLengthProb, MathContext.DECIMAL128).multiply(newPropDensity, MathContext.DECIMAL128);
            BigDecimal denominator = pathProb.multiply(pathLengthProb, MathContext.DECIMAL128).multiply(propDensity, MathContext.DECIMAL128);

            if (denominator.compareTo(BigDecimal.ZERO) == 0) {
                accept = new BigDecimal(0.0);
            } else {
                accept = numerator.divide(denominator, MathContext.DECIMAL128);
            }


        } // end if newPathLen = pathLen

        return accept;
    }

    public BigDecimal calculateRateAcceptanceProbability(double[] insParameters, double[] delParameters,
            double[] newInsParameters, double[] newDelParameters, double insRate, double delRate,
            double newInsRate, double newDelRate, BigDecimal pathProb, BigDecimal newPathProb) {

        double alpha, beta;
        // Probability of proposed rates
        // insertion rate
        alpha = insParameters[0];
        beta = insParameters[1];
        BigDecimal insPropDensity = getRateProposalDensity(insRate, alpha, beta);
        alpha = newInsParameters[0];
        beta = newInsParameters[1];
        BigDecimal revInsPropDensity = getRateProposalDensity(newInsRate, alpha, beta);//reverse move
        // deletion rate
        alpha = delParameters[0];
        beta = delParameters[1];
        BigDecimal delPropDensity = getRateProposalDensity(delRate, alpha, beta);
        alpha = newDelParameters[0];
        beta = newDelParameters[1];
        BigDecimal revDelPropDensity = getRateProposalDensity(newDelRate, alpha, beta);//reverse move


        // Acceptance probability
        BigDecimal accept;
        BigDecimal numerator = newPathProb.multiply(revInsPropDensity, MathContext.DECIMAL128).multiply(revDelPropDensity, MathContext.DECIMAL128);
        BigDecimal denominator = pathProb.multiply(insPropDensity, MathContext.DECIMAL128).multiply(delPropDensity, MathContext.DECIMAL128);


        // temporary .. gamma prior 1
//        BigDecimal prior = new BigDecimal (new Gamma(newInsRate*insParameters[0], 1.0/insParameters[0]).nextDouble() * new Gamma(newDelRate*delParameters[0], 1.0/delParameters[0]).nextDouble()/
//                (new Gamma(insRate*insParameters[0], 1.0/insParameters[0]).nextDouble() * new Gamma(delRate*delParameters[0], 1.0/delParameters[0]).nextDouble()));
//        numerator = numerator.multiply(prior);

//        // temporary .. gamma prior
//        BigDecimal prior = new BigDecimal (new Gamma(newInsRate*insParameters[1], 1.0/insParameters[1]).nextDouble() * new Gamma(newDelRate*delParameters[1], 1.0/delParameters[1]).nextDouble()/
//                (new Gamma(insRate*insParameters[1], 1.0/insParameters[1]).nextDouble() * new Gamma(delRate*delParameters[1], 1.0/delParameters[1]).nextDouble()));
//        numerator = numerator.multiply(prior);

//        // temporary .. gamma prior
//        BigDecimal prior = new BigDecimal (new Gamma(newInsRate, 1.0).nextDouble() * new Gamma(newDelRate, 1.0).nextDouble()/
//                (new Gamma(insRate, 1.0).nextDouble() * new Gamma(delRate, 1.0).nextDouble()));
//        numerator = numerator.multiply(prior);


        if (denominator.compareTo(BigDecimal.ZERO) == 0) {
            accept = new BigDecimal(0.0);
        } else {
            accept = numerator.divide(denominator, MathContext.DECIMAL128);
        }


        return accept;
    }

    public BigDecimal calculateParameterAcceptanceProbability(double[] insParameters, double[] delParameters,
            double[] newInsParameters, double[] newDelParameters, double insRate, double delRate,
            double newInsRate, double newDelRate, BigDecimal pathProb, BigDecimal newPathProb,
            double dpAlpha, double dpBeta, double dependenceProbability, double newDependenceProbability) {

        double alpha, beta;
        // Probability of proposed rates
        // insertion rate
        alpha = insParameters[0];
        beta = insParameters[1];
        BigDecimal insPropDensity = getRateProposalDensity(insRate, alpha, beta);
        alpha = newInsParameters[0];
        beta = newInsParameters[1];
        BigDecimal revInsPropDensity = getRateProposalDensity(newInsRate, alpha, beta);//reverse move
        // deletion rate
        alpha = delParameters[0];
        beta = delParameters[1];
        BigDecimal delPropDensity = getRateProposalDensity(delRate, alpha, beta);
        alpha = newDelParameters[0];
        beta = newDelParameters[1];
        BigDecimal revDelPropDensity = getRateProposalDensity(newDelRate, alpha, beta);//reverse move

        // Probability of proposed dependence probability
        BigDecimal DPPropDensity = getDPProposalDensity(dependenceProbability, dpAlpha, dpBeta);
        BigDecimal revDPPropDensity = getDPProposalDensity(newDependenceProbability, dpAlpha, dpBeta);//reverse move


        // Acceptance probability
        BigDecimal accept;
        BigDecimal numerator = newPathProb.multiply(revInsPropDensity, MathContext.DECIMAL128).multiply(revDelPropDensity, MathContext.DECIMAL128).multiply(revDPPropDensity, MathContext.DECIMAL128);
        BigDecimal denominator = pathProb.multiply(insPropDensity, MathContext.DECIMAL128).multiply(delPropDensity, MathContext.DECIMAL128).multiply(DPPropDensity, MathContext.DECIMAL128);

        if (denominator.compareTo(BigDecimal.ZERO) == 0) {
            accept = new BigDecimal(0.0);
        } else {
            accept = numerator.divide(denominator, MathContext.DECIMAL128);
        }


        return accept;
    }

    public BigDecimal calculateRateAcceptanceProbabilityPhylo(double[] insParameters, double[] delParameters,
            double[] newInsParameters, double[] newDelParameters, double insRate, double delRate,
            double newInsRate, double newDelRate, BigDecimal treeProb, BigDecimal newTreeProb) {

        double alpha, beta;
        // Probability of proposed rates
        // insertion rate
        alpha = insParameters[0];
        beta = insParameters[1];
        BigDecimal insPropDensity = getRateProposalDensityPhylo(insRate, alpha, beta);
        alpha = newInsParameters[0];
        beta = newInsParameters[1];
        BigDecimal revInsPropDensity = getRateProposalDensityPhylo(newInsRate, alpha, beta);//reverse move
        // deletion rate
        alpha = delParameters[0];
        beta = delParameters[1];
        BigDecimal delPropDensity = getRateProposalDensityPhylo(delRate, alpha, beta);
        alpha = newDelParameters[0];
        beta = newDelParameters[1];
        BigDecimal revDelPropDensity = getRateProposalDensityPhylo(newDelRate, alpha, beta);//reverse move


        // Acceptance probability
        BigDecimal accept;
        BigDecimal numerator = newTreeProb.multiply(revInsPropDensity, MathContext.DECIMAL128).multiply(revDelPropDensity, MathContext.DECIMAL128);
        BigDecimal denominator = treeProb.multiply(insPropDensity, MathContext.DECIMAL128).multiply(delPropDensity, MathContext.DECIMAL128);

        if (denominator.compareTo(BigDecimal.ZERO) == 0) {
            accept = new BigDecimal(0.0);
        } else {
            accept = numerator.divide(denominator, MathContext.DECIMAL128);
        }

//        System.out.println(">\t\t\t\t\t\t\t" + accept.toString());

//        if (newInsRate != 0.0 && newDelRate != 0.0) {

//            BigDecimal prior = new BigDecimal(new Gamma(newInsRate, insParameters[0]).nextDouble() * new Gamma(newDelRate, delParameters[0]).nextDouble() /
//                    (new Gamma(insRate, insParameters[0]).nextDouble() * new Gamma(delRate, delParameters[0]).nextDouble()));
        //accept = accept.multiply(prior);

//            System.out.println(">\t" + prior.toString());

        // inverse chi sq
        //double dPrior = (Math.pow(insRate,2) /Math.pow(newInsRate,2))*(Math.pow(delRate,2) /Math.pow(newDelRate,2)) * Math.exp(-0.5*(1/newInsRate - 1/insRate)-0.5*(1/newDelRate - 1/delRate));
        // chi sq
//        double dPrior = Math.exp(-0.5*(newInsRate - insRate)-0.5*(newDelRate - delRate));
//        accept = accept.multiply(new BigDecimal(Math.abs(dPrior)));
//        }

//        System.out.println(">\t" + accept.toString());

        return accept;
    }

    public BigDecimal calculateParameterAcceptanceProbabilityPhylo(double[] insParameters, double[] delParameters,
            double[] newInsParameters, double[] newDelParameters, double insRate, double delRate,
            double newInsRate, double newDelRate, BigDecimal treeProb, BigDecimal newTreeProb,
            double dpAlpha, double dpBeta, double dependenceProbability, double newDependenceProbability) {

        double alpha, beta;
        // Probability of proposed rates
        // insertion rate
        alpha = insParameters[0];
        beta = insParameters[1];
        BigDecimal insPropDensity = getRateProposalDensityPhylo(insRate, alpha, beta);
        alpha = newInsParameters[0];
        beta = newInsParameters[1];
        BigDecimal revInsPropDensity = getRateProposalDensityPhylo(newInsRate, alpha, beta);//reverse move
        // deletion rate
        alpha = delParameters[0];
        beta = delParameters[1];
        BigDecimal delPropDensity = getRateProposalDensityPhylo(delRate, alpha, beta);
        alpha = newDelParameters[0];
        beta = newDelParameters[1];
        BigDecimal revDelPropDensity = getRateProposalDensityPhylo(newDelRate, alpha, beta);//reverse move

        // Probability of proposed dependence probability
        BigDecimal DPPropDensity = getDPProposalDensity(dependenceProbability, dpAlpha, dpBeta);
        BigDecimal revDPPropDensity = getDPProposalDensity(newDependenceProbability, dpAlpha, dpBeta);//reverse move


        // Acceptance probability
        BigDecimal accept;
        BigDecimal numerator = newTreeProb.multiply(revInsPropDensity, MathContext.DECIMAL128).multiply(revDelPropDensity, MathContext.DECIMAL128).multiply(revDPPropDensity, MathContext.DECIMAL128);
        BigDecimal denominator = treeProb.multiply(insPropDensity, MathContext.DECIMAL128).multiply(delPropDensity, MathContext.DECIMAL128).multiply(DPPropDensity, MathContext.DECIMAL128);

        if (denominator.compareTo(BigDecimal.ZERO) == 0) {
            accept = new BigDecimal(0.0);
        } else {
            accept = numerator.divide(denominator, MathContext.DECIMAL128);
        }

//        System.out.println(">\t\t\t\t\t\t\t" + accept.toString());

//        if (newInsRate != 0.0 && newDelRate != 0.0) {

//            BigDecimal prior = new BigDecimal(new Gamma(newInsRate, insParameters[0]).nextDouble() * new Gamma(newDelRate, delParameters[0]).nextDouble() /
//                    (new Gamma(insRate, insParameters[0]).nextDouble() * new Gamma(delRate, delParameters[0]).nextDouble()));
        //accept = accept.multiply(prior);

//            System.out.println(">\t" + prior.toString());

        // inverse chi sq
        //double dPrior = (Math.pow(insRate,2) /Math.pow(newInsRate,2))*(Math.pow(delRate,2) /Math.pow(newDelRate,2)) * Math.exp(-0.5*(1/newInsRate - 1/insRate)-0.5*(1/newDelRate - 1/delRate));
        // chi sq
//        double dPrior = Math.exp(-0.5*(newInsRate - insRate)-0.5*(newDelRate - delRate));
//        accept = accept.multiply(new BigDecimal(Math.abs(dPrior)));
//        }

//        System.out.println(">\t" + accept.toString());

        return accept;
    }

    public BigDecimal calculateDPAcceptanceProbability(double dpAlpha, double dpBeta,
            double dependenceProbability, double newDependenceProbability,
            BigDecimal pathProb, BigDecimal newPathProb) {

        // Probability of proposed dependence probability
        BigDecimal propDensity = getDPProposalDensity(dependenceProbability, dpAlpha, dpBeta);
        BigDecimal revPropDensity = getDPProposalDensity(newDependenceProbability, dpAlpha, dpBeta);//reverse move

        // Acceptance probability
        BigDecimal accept;
        BigDecimal numerator = newPathProb.multiply(revPropDensity, MathContext.DECIMAL128);
        BigDecimal denominator = pathProb.multiply(propDensity, MathContext.DECIMAL128);


        if (denominator.compareTo(BigDecimal.ZERO) == 0) {
            accept = new BigDecimal(0.0);
        } else {
            accept = numerator.divide(denominator, MathContext.DECIMAL128);
        }


        return accept;
    }

    public BigDecimal getPathLengthProbability(int pathLen, int newPathLen, int minLen) {
        /**
         *  Probability of switching to another model using given operation
         */

        // Get the probabilities
        BigDecimal[] P = pathLengthProposalProbability(pathLen, minLen);

        // find the operation applied
        byte operation;
        if (pathLen < newPathLen) {
            operation = MOVE_ADD_EVENTS;
        } else if (pathLen > newPathLen) {
            operation = MOVE_DEL_EVENTS;
        } else {
            operation = MOVE_SWAP_EVENTS;
        }

        // return log probability    
        return P[operation];

    }

    public BigDecimal getPathProposalDensity(int[] path, int[] newPath, int edgeCount) {

        BigDecimal prob;
        int pathLen = path.length;
        int newPathLen = newPath.length;

        // Get the 128 Decimal Math Context
        MathContext mc = MathContext.DECIMAL128;

        // Get the pathProposal density depending on the operation performed
        if (newPathLen > pathLen) {// Events were added
            // get the edge frequencies
            int[] edgeFreq = new int[edgeCount];
            for (int i = 0; i < path.length; i++) {
                edgeFreq[path[i] - 1]++;
            }

            int denominator = 0;
            for (int i = 0; i < edgeCount; i++) {
                int e = i + 1;
                if (edgeFreq[e - 1] == 0) { // the edge is not present in the path
                    denominator += (path.length + 1) * (path.length + 2);
                } else if (edgeFreq[e - 1] == 1) { // the edge is present only once in the path
                    denominator += path.length * (path.length + 1);
                } else { // the edge is present more than once
                    // create a path with streches of edge 'e' removed
                    int[] temp = new int[path.length];
                    int idx = 0;
                    int stretchCount = 0;
                    boolean newStretch = true;
                    temp[idx] = path[idx++];
                    for (int j = 1; j < path.length; j++) {
                        if (path[j] == e && path[j] == path[j - 1]) {// check if the edge is same as the previous position
                            if (newStretch) {
                                stretchCount++;
                                newStretch = false;
                            }
                        } else { // different edge
                            temp[idx++] = path[j];
                            newStretch = true;
                        }
                    } // end while
                    // copy the path to a new matrix with exact size
                    int[] path1 = new int[idx];
                    for (int j = 0; j < idx; j++) {
                        path1[j] = temp[j];
                    }
                    // now count the positions in this new path where the edge 'e' does not occur
                    int count = 0;
                    for (int j = 0; j < path1.length; j++) {
                        if (path1[j] != e) {
                            count++;
                        }
                    }
                    denominator += (count + 1) * (count + 2);
                } // end if edge is absent, present once or occurs multiples times
            } // end for each edge
            denominator /= 2;

            if (denominator == 0) {
                prob = new BigDecimal(0.0);
            } else {
                prob = new BigDecimal(1.0, mc).divide(new BigDecimal(denominator), mc);
            }
        } else if (newPathLen == pathLen) { // Events were swapped
            // get the edge frequencies
            int[] edgeFreq = new int[edgeCount];
            for (int i = 0; i < path.length; i++) {
                edgeFreq[path[i] - 1]++;
            }

            int denominator = 0;
            for (int i = 0; i < edgeFreq.length; i++) {
                if (edgeFreq[i] == pathLen) { // does all events the path correspond to the same edge?
                    denominator = 0; // swapping not possible
                    break;
                }
                denominator += edgeFreq[i] * (pathLen - edgeFreq[i]);
            }
            denominator /= 2;

            if (denominator == 0) {
                prob = new BigDecimal(0.0);
            } else {
                prob = new BigDecimal(1.0, mc).divide(new BigDecimal(denominator), mc);
            }
        } else { // Events were deleted
            // create a path with streches of same edge removed
            int[] temp = new int[path.length];
            int idx = 0;
            int stretchCount = 0;
            boolean newStretch = true;
            temp[idx] = path[idx++];
            for (int i = 1; i < path.length; i++) {
                if (path[i] == path[i - 1]) {// check if the edge is same as the previous position
                    if (newStretch) {
                        stretchCount++;
                        newStretch = false;
                    }
                } else { // different edge
                    temp[idx++] = path[i];
                    newStretch = true;
                }
            } // end while
            // copy the path to a new matrix with exact size
            int[] path1 = new int[idx];
            for (int i = 0; i < idx; i++) {
                path1[i] = temp[i];
            }

            // get the edge frequencies
            int[] edgeFreq = new int[edgeCount];
            for (int i = 0; i < path1.length; i++) {
                edgeFreq[path1[i] - 1]++;
            }

            int denominator = 0;
            for (int i = 0; i < edgeFreq.length; i++) {
                denominator += edgeFreq[i] * Math.max(edgeFreq[i] - 1, 0);
            }
            denominator /= 2;
            denominator += stretchCount;

            if (denominator == 0) {
                prob = new BigDecimal(0.0);
            } else {
                prob = new BigDecimal(1.0, mc).divide(new BigDecimal(denominator), mc);
            }
        }

        return prob;
    }

    public BigDecimal getRateProposalDensity(double rate, double alpha, double beta) {
        /**
         * g(x,alpha, beta) is proportional to x^(alpha-1) exp(-beta*x)
         */
        double prob = Math.pow(rate, alpha - 1) * Math.exp(-beta * rate); // gamma
//        double prob = 1/beta * Math.exp(-Math.pow((rate - alpha)/beta,2)); // normal
        return new BigDecimal(prob);
    }

    public BigDecimal getRateProposalDensityPhylo(double rate, double alpha, double beta) {
        /**
         * g(x,alpha, beta) is proportional to x^(alpha-1) exp(-beta*x)
         */
        double prob = Math.pow(rate, alpha - 1) * Math.exp(-beta * rate); // gamma
        //double prob = (1/beta) *  Math.exp(- Math.pow(rate - alpha,2) /(2*beta*beta)); // normal
        return new BigDecimal(prob);
    }

    public BigDecimal getDPProposalDensity(double dp, double alpha, double beta) {
        /**
         * beta(x,alpha, beta) is proportional to x^(alpha-1) (1 - x)^(beta-1)
         */
        double prob = Math.pow(dp, alpha - 1) * Math.pow(1 - dp, beta - 1); // beta
        return new BigDecimal(prob);
    }

    private int[] generatePath(MetabolicNetwork frmNetwork, MetabolicNetwork toNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, double insRate, double delRate,
            double dependenceProbability) {

        return generatePath(frmNetwork, toNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, -1, dependenceProbability);
    }

    private int[] generatePath(MetabolicNetwork frmNetwork, MetabolicNetwork toNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, double insRate, double delRate,
            int extraEvents, double dependenceProbability) {


        int[] path;

        if (extraEvents < 0) {
            extraEvents = 5;
        }

        ArrayList lstEdges = frmNetwork.getDirectedReactions();
        int edgeCount = lstEdges.size();

        // number of alterable edges
        int nAlterableEdges = edgeCount - coreEdges.size() - prohibEdges.size();
        // find the edges that are different
        ArrayList lstDifferences = MetabolicNetwork.findDifferences(frmNetwork, toNetwork);

        // repeat until a valid path is generated
        while (true) {
            // generate random number of redundant events in the path (10 max. events)
            int nExtraEdges = (int) Math.floor(Math.random() * extraEvents);
            // initialise a new path;
            path = new int[nExtraEdges * 2 + lstDifferences.size()];

            // indices available for assigning events
            int[] availableIdx = new int[path.length];
            for (int i = 0; i < availableIdx.length; i++) {
                availableIdx[i] = i;
            }

            for (int i = 0; i < lstDifferences.size(); i++) {
                DirectedReaction differentEdge = (DirectedReaction) lstDifferences.get(i);
                // select a random position
                int idx = (int) Math.floor(Math.random() * availableIdx.length);
                // assign the current edge as an event
                path[availableIdx[idx]] = lstEdges.indexOf(differentEdge) + 1;

                // remove this index from the available indices
                availableIdx = Utilities.removeElement(availableIdx, idx);
            }

            // list of alterable edges
            int[] edges = new int[edgeCount];
            Iterator it = coreEdges.iterator();
            while (it.hasNext()) {
                // get the next edge
                DirectedReaction theEdge = (DirectedReaction) it.next();
                // mark this edge as already selected
                edges[lstEdges.indexOf(theEdge)] = 1;
            } // while it.hasNext()
            it = prohibEdges.iterator();
            while (it.hasNext()) {
                // get the next edge
                DirectedReaction theEdge = (DirectedReaction) it.next();
                // mark this edge as already selected
                edges[lstEdges.indexOf(theEdge)] = 1;
            } // while it.hasNext()

            int aeIdx = 0;
            int[] alterableEdges = new int[nAlterableEdges];
            for (int i = 0; i < edges.length; i++) {
                if (edges[i] == 0) {
                    alterableEdges[aeIdx++] = i + 1;
                }
            }

            // get the edge probabilities
            double[] edgeProbs = getNetworkEvolver().calculateHyperedgeProbabilities(frmNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, dependenceProbability);
            // extract the probabilities for alterable edges
            double[] alterableEdgeProbs = new double[alterableEdges.length];
            for (int i = 0; i < alterableEdges.length; i++) {
                alterableEdgeProbs[i] = edgeProbs[alterableEdges[i] - 1];
            }

            for (int i = 0; i < nExtraEdges; i++) {

                // select a random edge                
//                int eIdx = (int) Math.floor(Math.random() * nAlterableEdges);
                int eIdx = Utilities.randsample(alterableEdgeProbs);
                int eventEdge = alterableEdges[eIdx];

                // select a random position
                int idx = (int) Math.floor(Math.random() * availableIdx.length);
                // assign the current edge as an event
                path[availableIdx[idx]] = eventEdge;
                // remove this index from the available indices
                availableIdx = Utilities.removeElement(availableIdx, idx);

                // we need to cancel the effect of this effect. Add the same edge again
                // select a random position
                idx = (int) Math.floor(Math.random() * availableIdx.length);
                // assign the current edge as an event
                path[availableIdx[idx]] = eventEdge;
                // remove this index from the available indices
                availableIdx = Utilities.removeElement(availableIdx, idx);
            }

            if (isValidPath(path, frmNetwork, coreEdges, prohibEdges)) {
                break;
            }
        } // end while true

        return path;
    }

    private double[] proposeNewRates(int[] path, MetabolicNetwork startNetwork, ArrayList coreEdges,
            ArrayList prohibEdges, double[] insParameters, double[] delParameters, int evolutionModel) {

        int[] nEvents = new int[2];
        int[] nTotalEvents = new int[2];
        double[] newRates = new double[2];

        // get the rate parameters
        getRateParameters(path, startNetwork, coreEdges, prohibEdges, nEvents, nTotalEvents);
        double nIns = (double) nEvents[0];
        double nDel = (double) nEvents[1];
        double totInsCount = (double) nTotalEvents[0];
        double totDelCount = (double) nTotalEvents[1];

//        Gamma gamma = new Gamma(1, 1, RandomEngine.makeDefault());

//        newRates[0] = Utilities.round(insRate + new Random().nextGaussian()*nIns/totInsCount + insRate, decimalPlaces); // insertion rate
//        newRates[1] = Utilities.round(delRate + new Random().nextGaussian()*nDel/totDelCount + delRate, decimalPlaces); // deletion rate
//        // set the parameters
//        insParameters[0] = insRate;
//        insParameters[1] = nIns/totInsCount;
//        delParameters[0] = delRate;
//        delParameters[1] = nDel/totDelCount;

//        newRates[0] = Utilities.round(gamma.nextDouble(nIns+1, totInsCount), decimalPlaces); // insertion rate
//        newRates[1] = Utilities.round(gamma.nextDouble(nDel+1, totDelCount), decimalPlaces); // deletion rate
//        // set the parameters
//        insParameters[0] = nIns+1;
//        insParameters[1] = totInsCount;
//        delParameters[0] = nDel+1;
//        delParameters[1] = totDelCount;

        newRates[0] = Utilities.round(gammaDist.nextDouble(nIns + 1, totInsCount / path.length), decimalPlaces); // insertion rate
        newRates[1] = Utilities.round(gammaDist.nextDouble(nDel + 1, totDelCount / path.length), decimalPlaces); // deletion rate
        // set the parameters
        insParameters[0] = nIns + 1;
        insParameters[1] = totInsCount / path.length;
        delParameters[0] = nDel + 1;
        delParameters[1] = totDelCount / path.length;

//// working .. kind of      
//        newRates[0] = Utilities.round(gamma.nextDouble(nIns+path.length+1, totInsCount/path.length), decimalPlaces); // insertion rate
//        newRates[1] = Utilities.round(gamma.nextDouble(nDel+path.length+1, totDelCount/path.length), decimalPlaces); // deletion rate
//        // set the parameters
//        insParameters[0] = nIns+path.length+1;
//        insParameters[1] = totInsCount/path.length;
//        delParameters[0] = nDel+path.length+1;
//        delParameters[1] = totDelCount/path.length;

//        // set the parameters
//        // insertion
//        insParameters[0] = nIns+path.length+1;
//        //insParameters[0] = nIns + k;
//        insParameters[1] = totInsCount/path.length;
//        // deletion
//        delParameters[0] = nDel+path.length+1;
//        //delParameters[0] = nDel + k;
//        delParameters[1] = totDelCount/path.length;
        return newRates;
    }

    private String getRatesString(double insRate, double delRate) {
        return Double.toString(insRate) + "_" + Double.toString(delRate);
    }

    private String getRatesString(double insRate, double delRate, double dependenceProbability) {
        return Double.toString(insRate) + "_" + Double.toString(delRate) + "_" + Double.toString(dependenceProbability);
    }

    private String getPathRatesString(int[] path, double insRate, double delRate) {
        return Utilities.toString(path) + "_" + Double.toString(insRate) + "_" + Double.toString(delRate);
    }

    private String getPathRatesString(int[] path, double insRate, double delRate, double dependenceProbability) {
        return Utilities.toString(path) + "_" + Double.toString(insRate) + "_" + Double.toString(delRate) + "_" + Double.toString(dependenceProbability);
    }

    private String getNetworkRatesString(MetabolicNetwork theNetwork, double insRate, double delRate) {
        return Utilities.toString(theNetwork.getReactionSequence(), "") + "_" + Double.toString(insRate) + "_" + Double.toString(delRate);
    }

    private String getNetworkRatesString(Byte[] theNetwork, double insRate, double delRate) {
        return Utilities.toString(theNetwork, "") + "_" + Double.toString(insRate) + "_" + Double.toString(delRate);
    }

    private String getNetworksRatesString(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, double insRate, double delRate, double dependenceProbability) {
        return Utilities.toString(startNetwork.getReactionSequence(), "") + "_" + Utilities.toString(endNetwork.getReactionSequence(), "") + "_" + Double.toString(insRate) + "_" + Double.toString(delRate) + "_" + Double.toString(dependenceProbability);
    }

    private String getNetworksRatesString(Byte[] startNetwork, Byte[] endNetwork, double insRate, double delRate, double dependenceProbability) {
        return Utilities.toString(startNetwork, "") + "_" + Utilities.toString(endNetwork, "") + "_" + Double.toString(insRate) + "_" + Double.toString(delRate) + "_" + Double.toString(dependenceProbability);
    }

    private String getNetworksPathRatesString(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, int[] path, double insRate, double delRate) {
        return Utilities.toString(startNetwork.getReactionSequence(), "") + "_" + Utilities.toString(endNetwork.getReactionSequence(), "") + "_" + Utilities.toString(path) + "_" + Double.toString(insRate) + "_" + Double.toString(delRate);
    }

    private String getNetworksPathRatesString(MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, int[] path, double insRate, double delRate, double dependenceProbability) {
        return Utilities.toString(startNetwork.getReactionSequence(), "") + "_" + Utilities.toString(endNetwork.getReactionSequence(), "") + "_" + Utilities.toString(path) + "_" + Double.toString(insRate) + "_" + Double.toString(delRate) + "_" + Double.toString(dependenceProbability);
    }

    private String getNetworksString(Byte[] startNetwork, Byte[] endNetwork) {
        return Utilities.toString(startNetwork, "") + "_" + Utilities.toString(endNetwork, "");
    }

    private String getTreeRatesString(PhyloTree tree, double insRate, double delRate, double dependenceProbability) {
        return getTreeString(tree.getRoot()) + "_" + Double.toString(insRate) + "_" + Double.toString(delRate) + "_" + Double.toString(dependenceProbability);
    }

    private String getTreeString(PhyloNode phylonode) {
        String leftString = "", rightString = "";
        if (phylonode.getLeftSon() != null) // get string on left sub tree
        {
            leftString = getTreeString(phylonode.getLeftSon());
        }
        if (phylonode.getRightSon() != null) // get string on right sub tree
        {
            rightString = getTreeString(phylonode.getRightSon());
        }

        // get current network
        Byte[] currentNetworkSeq;
        if (phylonode.isLeaf()) {
            currentNetworkSeq = phylonode.getMetabolicNetwork().getReactionSequence();
            return Utilities.toString(currentNetworkSeq, "");
        } else {
            currentNetworkSeq = ((MCMCOutput) phylonode.getData()).getCurrentNetwork().getReactionSequence();
            return "(" + leftString + "," + Utilities.toString(currentNetworkSeq, "") + "," + rightString + ")";
        }
    }

    public HashMap calculateLikelihood(ArrayList pathList, BigDecimal[] pathProbList, Double[] insRates, Double[] delRates, int nBurning, boolean printOutput) {
        return calculateLikelihood(null, null, null, -1, pathList, null, null, pathProbList, insRates, delRates, nBurning, printOutput);
    }

    private HashMap calculateLikelihood(MetabolicNetwork refNetwork, ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, ArrayList pathList, ArrayList startNetworksList, ArrayList endNetworksList, BigDecimal[] pathProbList, Double[] insRates, Double[] delRates, int nBurning, boolean printOutput) {

        HashSet processedRatesPaths = new HashSet((int) (pathList.size() * 0.4));
        HashSet processedRates = new HashSet((int) (pathList.size() * 0.4));
        HashMap likelihood = new HashMap((int) (pathList.size() * 0.4));
        int itIdx = -1;
        // process all paths 
        Iterator itPath = pathList.iterator();
        Iterator itStartNetwork = null, itEndNetwork = null;
        if (refNetwork != null) {
            itStartNetwork = startNetworksList.iterator();
            itEndNetwork = endNetworksList.iterator();
        }
        while (itPath.hasNext()) {
            int[] path = (int[]) itPath.next();
            if (++itIdx <= nBurning) {
                continue;
            }
            MetabolicNetwork startNetwork = null, endNetwork = null;
            if (refNetwork != null) {
                startNetwork = refNetwork.getNetwork((Byte[]) itStartNetwork.next());
                endNetwork = refNetwork.getNetwork((Byte[]) itEndNetwork.next());
            }

            Double insRate = insRates[itIdx];
            Double delRate = delRates[itIdx];
            String strRatesPath, strRates;
            if (refNetwork != null) {
                strRates = getNetworksRatesString(startNetwork, endNetwork, insRate, delRate, 1.0);
                strRatesPath = getNetworksPathRatesString(startNetwork, endNetwork, path, insRate, delRate);
            } else {
                strRates = getRatesString(insRate, delRate);
                strRatesPath = getPathRatesString(path, insRates[itIdx], delRates[itIdx]);
            }

            if (processedRatesPaths.add(strRatesPath)) {
                // we haven't seen this rates-path(-networks) combination before .. process it

                BigDecimal lh = pathProbList[itIdx];
//                if (this.lstDistinctNetworks != null && startNetwork != null) {
//                    String strNetworkRates = getNetworkRatesString(startNetwork, insRate, delRate);
//                    Double eqLogProb = (Double) this.lstDistinctNetworks.get(strNetworkRates);
//                    double eqProb;
//                    if (eqLogProb != null) {
//                        eqProb = Math.exp(eqLogProb.doubleValue());
//                    } else {
//                        eqProb = Math.exp(new NetworkEvolver().approximateEquilibriumProbability(startNetwork, refNetwork, coreEdges, prohibEdges, evolutionModel, insRate.doubleValue(), delRate.doubleValue()));
//                    }
//                    lh = lh.multiply(new BigDecimal(eqProb));
//                }

                if (!processedRates.add(strRates)) {
                    // we have seen this rates(-networks) combination before. Add to the likelihood
                    BigDecimal currLH = (BigDecimal) likelihood.get(strRates);
                    lh = lh.add(currLH);
                }

                // add the lh to the likelihood table .. it will replace the old value (if one exists)
                likelihood.put(strRates, lh);
            }
        } // while it.hasNext()


        HashMap likelihoodTable = new HashMap(likelihood.size());
        // process all likelihood entries
        Iterator itKey = likelihood.keySet().iterator();
        while (itKey.hasNext()) {
            String strRate = (String) itKey.next();
            BigDecimal lh = (BigDecimal) likelihood.get(strRate);

            //Double[] rates = new Double[]{Double.valueOf(strRates[strRates.length - 2]), Double.valueOf(strRates[strRates.length - 1])};
            likelihoodTable.put(strRate, lh);

            String[] strRates = strRate.split("_");
            if (printOutput) {
                if (refNetwork != null) {
                    System.out.print(strRates[0] + "\t" + strRates[1] + "\t");
                    System.out.print(strRates[2].toString() + "\t" + strRates[3].toString());
                } else {
                    System.out.print(strRates[0].toString() + "\t" + strRates[1].toString());
                }
                System.out.println("\t" + lh.round(MathContext.DECIMAL128).toString());
            }
        } // while itKey.hasNext()

        return likelihoodTable;
    }

    public HashMap calculateLikelihoodDP(ArrayList pathList, BigDecimal[] pathProbList, Double[] insRates, Double[] delRates, Double[] dependenceProbabilities, int nBurning, boolean printOutput) {
        return calculateLikelihoodDP(null, null, null, -1, pathList, null, null, pathProbList, insRates, delRates, dependenceProbabilities, nBurning, printOutput);
    }

    private HashMap calculateLikelihoodDP(MetabolicNetwork refNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, ArrayList pathList,
            ArrayList startNetworksList, ArrayList endNetworksList, BigDecimal[] pathProbList,
            Double[] insRates, Double[] delRates, Double[] dependenceProbabilities,
            int nBurning, boolean printOutput) {

        HashSet processedRatesPaths = new HashSet((int) (pathList.size() * 0.4));
        HashSet processedRates = new HashSet((int) (pathList.size() * 0.4));
        HashMap likelihood = new HashMap((int) (pathList.size() * 0.4));
        int itIdx = -1;
        // process all paths
        Iterator itPath = pathList.iterator();
        Iterator itStartNetwork = null, itEndNetwork = null;
        if (refNetwork != null) {
            itStartNetwork = startNetworksList.iterator();
            itEndNetwork = endNetworksList.iterator();
        }
        while (itPath.hasNext()) {
            int[] path = (int[]) itPath.next();
            if (++itIdx <= nBurning) {
                continue;
            }
            MetabolicNetwork startNetwork = null, endNetwork = null;
            if (refNetwork != null) {
                startNetwork = refNetwork.getNetwork((Byte[]) itStartNetwork.next());
                endNetwork = refNetwork.getNetwork((Byte[]) itEndNetwork.next());
            }

            Double insRate = insRates[itIdx];
            Double delRate = delRates[itIdx];
            Double dependenceProbability = dependenceProbabilities[itIdx];
            String strRatesPath, strRates;
            if (refNetwork != null) {
                strRates = getNetworksRatesString(startNetwork, endNetwork, insRate, delRate, dependenceProbability);
                strRatesPath = getNetworksPathRatesString(startNetwork, endNetwork, path, insRate, delRate, dependenceProbability);
            } else {
                strRates = getRatesString(insRate, delRate, dependenceProbability);
                strRatesPath = getPathRatesString(path, insRates[itIdx], delRates[itIdx], dependenceProbabilities[itIdx]);
            }

            if (processedRatesPaths.add(strRatesPath)) {
                // we haven't seen this rates-path(-networks) combination before .. process it

                BigDecimal lh = pathProbList[itIdx];
//                if (this.lstDistinctNetworks != null && startNetwork != null) {
//                    String strNetworkRates = getNetworkRatesString(startNetwork, insRate, delRate);
//                    Double eqLogProb = (Double) this.lstDistinctNetworks.get(strNetworkRates);
//                    double eqProb;
//                    if (eqLogProb != null) {
//                        eqProb = Math.exp(eqLogProb.doubleValue());
//                    } else {
//                        eqProb = Math.exp(new NetworkEvolver().approximateEquilibriumProbability(startNetwork, refNetwork, coreEdges, prohibEdges, evolutionModel, insRate.doubleValue(), delRate.doubleValue()));
//                    }
//                    lh = lh.multiply(new BigDecimal(eqProb));
//                }

                if (!processedRates.add(strRates)) {
                    // we have seen this rates(-networks) combination before. Add to the likelihood
                    BigDecimal currLH = (BigDecimal) likelihood.get(strRates);
                    lh = lh.add(currLH);
                }

                // add the lh to the likelihood table .. it will replace the old value (if one exists)
                likelihood.put(strRates, lh);
            }
        } // while it.hasNext()


        HashMap likelihoodTable = new HashMap(likelihood.size());
        // process all likelihood entries
        Iterator itKey = likelihood.keySet().iterator();
        while (itKey.hasNext()) {
            String strRate = (String) itKey.next();
            BigDecimal lh = (BigDecimal) likelihood.get(strRate);

            //Double[] rates = new Double[]{Double.valueOf(strRates[strRates.length - 2]), Double.valueOf(strRates[strRates.length - 1])};
            likelihoodTable.put(strRate, lh);

            String[] strRates = strRate.split("_");
            if (printOutput) {
                if (refNetwork != null) {
                    System.out.print(strRates[0] + "\t" + strRates[1] + "\t");
                    System.out.print(strRates[2].toString() + "\t" + strRates[3].toString() + "\t" + strRates[4].toString());
                } else {
                    System.out.print(strRates[0].toString() + "\t" + strRates[1].toString() + "\t" + strRates[2].toString());
                }
                System.out.println("\t" + lh.round(MathContext.DECIMAL128).toString());
            }
        } // while itKey.hasNext()

        return likelihoodTable;
    }

    public HashMap calculateLikelihood(PhyloTree phylotree, MetabolicNetwork refNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, int nBurning,
            int subNetworkSize, boolean printOutput, double dependenceProbability) {

        // get indices of core and prohibited edges
        Iterator itCore = coreEdges.iterator();
        Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iCoreEdges.length; i++) {
            iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itCore.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itCore.next();
            iCoreEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }
        Iterator itProhib = prohibEdges.iterator();
        Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iProhibEdges.length; i++) {
            iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itProhib.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itProhib.next();
            iProhibEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }

        // extract insertion and deletion rates from the root
        MCMCOutput mcmcOutput = (MCMCOutput) phylotree.getRoot().getData();
        Double[] insRates = mcmcOutput.getInsRates();
        Double[] delRates = mcmcOutput.getDelRates();

        return calculateLikelihood(phylotree, refNetwork, iCoreEdges, iProhibEdges, evolutionModel,
                insRates, delRates, nBurning, subNetworkSize, printOutput, dependenceProbability);

    }

    private HashMap calculateLikelihood(PhyloTree phylotree, MetabolicNetwork refNetwork,
            Byte[] iCoreEdges, Byte[] iProhibEdges, int evolutionModel, Double[] insRates,
            Double[] delRates, int nBurning, int subNetworkSize, boolean printOutput,
            double dependenceProbability) {

        distinctNetworks = new HashMap((int) (insRates.length * 0.2)); // Assumption: Number of networks visited << number of iterations
        HashSet processedRates = new HashSet((int) (insRates.length * 0.2));
        HashMap likelihood = new HashMap((int) (insRates.length * 0.2));

        double insRate, delRate;
        for (int iter = nBurning + 1; iter < insRates.length - 1; iter++) {

            // get current insertion and deletion rates
            insRate = insRates[iter];
            delRate = delRates[iter];

            String strRates = getRatesString(insRate, delRate);
            if (processedRates.add(strRates)) { // we havent seen this rate combination before
//                System.out.println(strRates);

                // calculate the likelihood for this rate combination
                Double lh = calculateLikelihood(phylotree.getRoot(), refNetwork, null, iCoreEdges, iProhibEdges,
                        insRate, delRate, evolutionModel, nBurning, subNetworkSize, dependenceProbability);

                // add the lh to the likelihood table 
                likelihood.put(strRates, lh);

                if (printOutput) {
                    System.out.print(Double.toString(insRate) + "\t" + Double.toString(delRate));
                    System.out.println("\t" + lh.toString());
                }

            } // if unseen rate combination
        } // end for

        return likelihood;
    }

    private double calculateLikelihood(PhyloNode phylonode, MetabolicNetwork refNetwork,
            Byte[] parentNetwork, Byte[] iCoreEdges, Byte[] iProhibEdges,
            double insRate, double delRate, int evolutionModel, int nBurning, int subNetworkSize,
            double dependenceProbability) {

        HashSet processedNetworks = new HashSet();
        HashSet processedNetworkPairs = new HashSet();

        HashSet lstNetworks;
        if (phylonode.isLeaf()) {
            lstNetworks = new HashSet();
            lstNetworks.add(phylonode.getMetabolicNetwork().getReactionSequence());
        } else {
            MCMCOutput mcmcOutput = (MCMCOutput) phylonode.getData();
            lstNetworks = new HashSet(mcmcOutput.getNetworks().subList(nBurning + 1, mcmcOutput.getNetworks().size()));
        }

        // get evol time (branch length)
        double parentEvolTime;
        if (parentNetwork == null) {
            parentEvolTime = -1.0;
        } else if (phylonode.isLeftSon()) {
            parentEvolTime = phylonode.getParent().getLeftBranchLength();
        } else {
            parentEvolTime = phylonode.getParent().getRightBranchLength();
        }

        // get a local copy
        NetworkEvolver nwEvolver = getNetworkEvolver();
        Iterator it = lstNetworks.iterator();
        double treeProb = 0.0;
        while (it.hasNext()) {
            // get the current network
            Byte[] currentNetwork = (Byte[]) it.next();

            if (!processedNetworks.add(Utilities.toString(currentNetwork, ""))) {
                continue;
            }

            Double parentTP;
            if (parentNetwork == null) { // this is the root .. calculate eq. probability, if asked
                String strNetworkRates;
                // calculate transition probability from parent to current network
                strNetworkRates = getNetworkRatesString(currentNetwork, insRate, delRate);
                parentTP = (Double) distinctNetworks.get(strNetworkRates);
                if (parentTP == null) {
                    parentTP = nwEvolver.approximateEquilibriumProbability(currentNetwork, refNetwork,
                            evolutionModel, insRate, delRate, subNetworkSize, dependenceProbability);
                    // add it to the list of distinct network pairs
                    distinctNetworks.put(strNetworkRates, parentTP);
                }
            } else {
                // process this pair if we havent seen it before
                String strNetworks = getNetworksString(parentNetwork, currentNetwork);
                if (!processedNetworkPairs.add(strNetworks)) {
                    continue;
                }

                // calculate transition probability from parent to current network
                String strNetworksRates = getNetworksRatesString(parentNetwork, currentNetwork, insRate, delRate, dependenceProbability);
                parentTP = (Double) distinctNetworkPairs.get(strNetworksRates);
                if (parentTP == null) {
                    boolean useOnlyVisitedNetworks = false;
                    if (useOnlyVisitedNetworks) {
                        parentTP = 1.0;
                        continue;
                    } else {
                        // Calculate the probability
                        parentTP = nwEvolver.approximateTransitionProbability(parentNetwork, currentNetwork,
                                refNetwork, iCoreEdges, iProhibEdges, evolutionModel, insRate, delRate,
                                parentEvolTime, subNetworkSize, dependenceProbability);
                        // add it to the list of distinct network pairs
                        distinctNetworkPairs.put(strNetworksRates, parentTP);
                    }
                }

            } // end if parentNetwork == null

            Double leftTP = 1.0, rightTP = 1.0;
            if (phylonode.getLeftSon() != null) // calculate transition probability from current network to Left son
            {
                leftTP = calculateLikelihood(phylonode.getLeftSon(), refNetwork, currentNetwork, iCoreEdges,
                        iProhibEdges, insRate, delRate, evolutionModel, nBurning, subNetworkSize, dependenceProbability);
            }
            if (phylonode.getRightSon() != null) // calculate transition probability from current network to Right son
            {
                rightTP = calculateLikelihood(phylonode.getRightSon(), refNetwork, currentNetwork, iCoreEdges,
                        iProhibEdges, insRate, delRate, evolutionModel, nBurning, subNetworkSize, dependenceProbability);
            }

            // add the log-values
            treeProb += Math.exp(Math.log(parentTP) + Math.log(leftTP) + Math.log(rightTP));
        } // end while


        return treeProb;
    }

    public HashMap calculateLikelihood(PhyloTree phylotree, MetabolicNetwork refNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, int nBurning,
            boolean printOutput, double dependenceProbability) {
        /** this function calculates the likelihood by using full network. 
         *  Working: For given rates,  it first calculates the exponential of the matrix and 
         *           saves it for future use. It then recursively calculates the likelihood of 
         *           the tree by summing over all possible networks
         *  Caution: Not to be used on large networks
         */

        // extract insertion and deletion rates from the root
        MCMCOutput mcmcOutput = (MCMCOutput) phylotree.getRoot().getData();
        Double[] insRates = mcmcOutput.getInsRates();
        Double[] delRates = mcmcOutput.getDelRates();

        HashSet processedRates = new HashSet((int) (insRates.length * 0.2));
        HashMap likelihood = new HashMap((int) (insRates.length * 0.2));

        double insRate, delRate;
        for (int iter = nBurning + 1; iter < insRates.length - 1; iter++) {

            // get current insertion and deletion rates
            insRate = insRates[iter];
            delRate = delRates[iter];

            String strRates = getRatesString(insRate, delRate);
            if (processedRates.add(strRates)) { // we havent seen this rate combination before
//                System.out.println(strRates);

                // calculate the likelihood for this rate combination
                Double lh = calculateLikelihood(phylotree, refNetwork, coreEdges, prohibEdges,
                        insRate, delRate, evolutionModel, nBurning, dependenceProbability);

                // add the lh to the likelihood table 
                likelihood.put(strRates, lh);

                if (printOutput) {
                    System.out.print(Double.toString(insRate) + "\t" + Double.toString(delRate));
                    System.out.println("\t" + lh.toString());
                }

            } // if unseen rate combination
        } // end for

        return likelihood;

    }

    public double calculateLikelihood(PhyloTree phylotree, MetabolicNetwork refNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate,
            int evolutionModel, int nBurning, double dependenceProbability) {
        /** this function calculates the likelihood by using full network. 
         *  Working: For given rates,  it first calculates the exponential of the matrix and 
         *           saves it for future use. It then recursively calculates the likelihood of 
         *           the tree by summing over all possible networks
         *  Caution: Not to be used on large networks
         */

        // get indices of core and prohibited edges
        Iterator itCore = coreEdges.iterator();
        Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iCoreEdges.length; i++) {
            iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itCore.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itCore.next();
            iCoreEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }
        Iterator itProhib = prohibEdges.iterator();
        Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iProhibEdges.length; i++) {
            iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        }
        while (itProhib.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itProhib.next();
            iProhibEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }

        NetworkEvolver nwEvolver = getNetworkEvolver();

        int edgeCount = refNetwork.getDirectedReactions().size();
        // array to store unalterable edges
        Byte[] unalterableEdges = new Byte[edgeCount];
        // define core and prohibited edges as unalterable edges.
        // all the rest are alterable edges.
        for (int i = 0; i < edgeCount; i++) {
            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            } else if (iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            } else {
                unalterableEdges[i] = -1;
            }
        } // end for all edges

        // enumerate all possible networks
        ArrayList lstAllNetworks = nwEvolver.enumerateAllNetworks(edgeCount, unalterableEdges);
        // create the rate matrix
        double[][] rateMatrix = nwEvolver.generateRateMatrix(lstAllNetworks, refNetwork, evolutionModel, insRate, delRate, dependenceProbability);

        // calculate the transition probability matrix
        double[][] transProbs = nwEvolver.calculateTransitionProbabilityMatrix(rateMatrix, 1.0);
        // calculate the equilibrium probabilities
        double[] eqProbs = nwEvolver.calculateEquilibriumProbabilityVector(rateMatrix);

        return calculateLikelihood(phylotree.getRoot(), refNetwork, null, insRate, delRate,
                nBurning, transProbs, eqProbs);
    }

    private double calculateLikelihood(PhyloNode phylonode, MetabolicNetwork refNetwork,
            Byte[] parentNetwork, double insRate, double delRate, int nBurning,
            double[][] transProbs, double[] eqProbs) {
        /** this function calculates the likelihood by using full network. 
         *  Working: For given rates,  it first calculates the exponential of the matrix and 
         *           saves it for future use. It then recursively calculates the likelihood of 
         *           the tree by summing over all possible networks
         *  Caution: Not to be used on large networks
         */
        HashSet processedNetworks = new HashSet();
        HashSet processedNetworkPairs = new HashSet();

        HashSet lstNetworks;
        if (phylonode.isLeaf()) {
            lstNetworks = new HashSet();
            lstNetworks.add(phylonode.getMetabolicNetwork().getReactionSequence());
        } else {
            MCMCOutput mcmcOutput = (MCMCOutput) phylonode.getData();
            lstNetworks = new HashSet(mcmcOutput.getNetworks().subList(nBurning + 1, mcmcOutput.getNetworks().size()));
        }

        Iterator it = lstNetworks.iterator();
        double treeProb = 0.0;
        while (it.hasNext()) {
            // get the current network
            Byte[] currentNetwork = (Byte[]) it.next();

            if (!processedNetworks.add(Utilities.toString(currentNetwork, ""))) {
                continue;
            }

            Double parentTP;
            if (parentNetwork == null) { // this is the root .. calculate eq. probability
                parentTP = eqProbs[NetworkEvolver.getNetworkIndex(currentNetwork)];
            } else {
                // process this pair if we havent seen it before
                String strNetworks = getNetworksString(parentNetwork, currentNetwork);
                if (!processedNetworkPairs.add(strNetworks)) {
                    continue;
                }

                /* AMT                // see if this network pair was visited
                String strNetworksRates = getNetworksRatesString(parentNetwork, currentNetwork, insRate, delRate);
                parentTP = (Double) distinctNetworkPairs.get(strNetworksRates);
                if (parentTP == null) {
                boolean useOnlyVisitedNetworks = true;
                if (useOnlyVisitedNetworks) {
                parentTP = 1.0;
                continue;
                }
                }
                 */
                // transition probability from parent to current network
                parentTP = transProbs[NetworkEvolver.getNetworkIndex(parentNetwork)][NetworkEvolver.getNetworkIndex(currentNetwork)];

            } // end if parentNetwork == null

            Double leftTP = 1.0, rightTP = 1.0;
            if (phylonode.getLeftSon() != null) // calculate transition probability from current network to Left son
            {
                leftTP = calculateLikelihood(phylonode.getLeftSon(), refNetwork, currentNetwork,
                        insRate, delRate, nBurning, transProbs, eqProbs);
            }
            if (phylonode.getRightSon() != null) // calculate transition probability from current network to Right son
            {
                rightTP = calculateLikelihood(phylonode.getRightSon(), refNetwork, currentNetwork,
                        insRate, delRate, nBurning, transProbs, eqProbs);
            }

            // add the log-values
            treeProb += Math.exp(Math.log(parentTP) + Math.log(leftTP) + Math.log(rightTP));
        } // end while


        return treeProb;
    }

    private HashMap consolidateLikelihood(HashMap networkLikelihood, MetabolicNetwork refNetwork, ArrayList coreEdges, ArrayList prohibEdges, int evolutionModel, boolean useEquillibriumProb, double dependenceProbability) {
        Set keySet = networkLikelihood.keySet();
        HashSet processedRates = new HashSet((int) (keySet.size() * 0.4));
        HashMap likelihood = new HashMap((int) (keySet.size() * 0.8));
        Iterator itKey = keySet.iterator();
        while (itKey.hasNext()) {
            // get next networks-rates combination
            String strNetworkRates = (String) itKey.next();
            // Extract the different values
            String[] networkRates = strNetworkRates.split("_");
            Double insRate = new Double(networkRates[2]);
            Double delRate = new Double(networkRates[3]);
            BigDecimal lh = (BigDecimal) networkLikelihood.get(strNetworkRates);
            MetabolicNetwork startNetwork;
            if (refNetwork != null) {
                startNetwork = refNetwork.getNetwork(MetabolicNetwork.getReactionSequence(networkRates[0]));

                if (useEquillibriumProb && this.distinctNetworks != null) {
                    Double eqLogProb = (Double) this.distinctNetworks.get(strNetworkRates);
                    double eqProb;
                    if (eqLogProb != null) {
                        eqProb = Math.exp(eqLogProb.doubleValue());
                    } else {
                        eqProb = Math.exp(getNetworkEvolver().approximateEquilibriumProbability(startNetwork, refNetwork, coreEdges, prohibEdges, evolutionModel, insRate.doubleValue(), delRate.doubleValue(), dependenceProbability));
                    }
                    lh = lh.multiply(new BigDecimal(eqProb));
                }
            }

            String strRates = getRatesString(insRate, delRate);
            if (!processedRates.add(strRates)) {
                // we have seen this rates combination before. Add the current network combination to the likelihood
                BigDecimal currLH = (BigDecimal) likelihood.get(strRates);
                lh = lh.add(currLH);
            }
            // add the lh to the likelihood table .. it will replace the old value (if one exists)
            likelihood.put(strRates, lh);

        } // end while it.hasNext()

        return likelihood;
    }

    private void mcmcDiagnostics(MCMCOutput mcmcOutput, int iter, int nBurning) {

//        if (iter < nBurning+2)
//            return;

        int lag = 10;
        double[] auto_corr = calculateParameterAutocorrelation(mcmcOutput, iter, lag);
        System.out.println("Autocorrelation: " + Double.toString(auto_corr[0]) + "\t" + Double.toString(auto_corr[1]));
    }

    public double[][] calculateAutocorrelation(MCMCOutput mcmcOutput, boolean printOutput) {

        int startLag = 1;
        int endLag = 200;

        double[][] auto_corr = new double[endLag - startLag + 1][5];
        for (int lag = startLag; lag <= endLag; lag++) {
            double[] corr = calculateParameterAutocorrelation(mcmcOutput, mcmcOutput.getInsRates().length, lag);
            double[] corr_path = calculatePathAutocorrelation(mcmcOutput, mcmcOutput.getInsRates().length, lag);

            // save the values
            int idx = lag - startLag;
            auto_corr[idx][0] = lag;
            auto_corr[idx][1] = corr[0];
            auto_corr[idx][2] = corr[1];
            auto_corr[idx][3] = corr_path[0];
            auto_corr[idx][4] = corr_path[1];
        }

        if (printOutput) {
            System.out.println("Autocorrelation:");
            System.out.println("Lag\tInsertion Rate\tDeletion Rate");
            for (int i = 0; i < auto_corr.length; i++) {
                System.out.println(Long.toString(Math.round(auto_corr[i][0])) + "\t" + Double.toString(auto_corr[i][1]) + "\t" + Double.toString(auto_corr[i][2]) + "\t" + Double.toString(auto_corr[i][3]) + "\t" + Double.toString(auto_corr[i][4]));
            }
            System.out.println();
        }
        return auto_corr;
    }

    public double[][] calculateAutocorrelationPhylo(MCMCOutput mcmcOutput, boolean printOutput) {

        int startLag = 1;
        int endLag = 200;

        double[][] auto_corr = new double[endLag - startLag + 1][5];
        for (int lag = startLag; lag <= endLag; lag++) {
            double[] corr = calculateParameterAutocorrelation(mcmcOutput, mcmcOutput.getInsRates().length, lag);

            // save the values
            int idx = lag - startLag;
            auto_corr[idx][0] = lag;
            auto_corr[idx][1] = corr[0];
            auto_corr[idx][2] = corr[1];
        }

        if (printOutput) {
            System.out.println("Autocorrelation:");
            System.out.println("Lag\tInsertion Rate\tDeletion Rate");
            for (int i = 0; i < auto_corr.length; i++) {
                System.out.println(Long.toString(Math.round(auto_corr[i][0])) + "\t" + Double.toString(auto_corr[i][1]) + "\t" + Double.toString(auto_corr[i][2]));
            }
            System.out.println();
        }
        return auto_corr;
    }

    private double[] calculateParameterAutocorrelation(MCMCOutput mcmcOutput, int iter, int lag) {

        Double[] insRates = mcmcOutput.getInsRates();
        Double[] delRates = mcmcOutput.getDelRates();

        // averages
        double ir_avg_1 = 0, ir_avg_2 = 0;
        double dr_avg_1 = 0, dr_avg_2 = 0;
        for (int i = 0; i < iter - lag; i++) {
            ir_avg_1 += insRates[i];
            dr_avg_1 += delRates[i];

            ir_avg_2 += insRates[i + lag];
            dr_avg_2 += delRates[i + lag];
        }
        int nSample = iter - lag;
        ir_avg_1 /= (double) nSample;
        dr_avg_1 /= (double) nSample;
        ir_avg_2 /= (double) nSample;
        dr_avg_2 /= (double) nSample;

        double ir_A = 0, ir_B1 = 0, ir_B2 = 0;
        double dr_A = 0, dr_B1 = 0, dr_B2 = 0;
        for (int i = 0; i < iter - lag; i++) {
            ir_A += (insRates[i].doubleValue() - ir_avg_1) * (insRates[i + lag].doubleValue() - ir_avg_2);
            dr_A += (delRates[i].doubleValue() - dr_avg_1) * (delRates[i + lag].doubleValue() - dr_avg_2);

            ir_B1 += Math.pow(insRates[i].doubleValue() - ir_avg_1, 2);
            dr_B1 += Math.pow(delRates[i].doubleValue() - dr_avg_1, 2);
            ir_B2 += Math.pow(insRates[i + lag].doubleValue() - ir_avg_2, 2);
            dr_B2 += Math.pow(delRates[i + lag].doubleValue() - dr_avg_2, 2);
        }

        double corr_ir = ir_A / Math.sqrt(ir_B1 * ir_B2);
        double corr_dr = dr_A / Math.sqrt(dr_B1 * dr_B2);

        return new double[]{corr_ir, corr_dr};
    }

    private double[] calculatePathAutocorrelation(MCMCOutput mcmcOutput, int iter, int lag) {

        BigDecimal[] pathProbs = mcmcOutput.getPathProbabilities();
        ArrayList pathList = mcmcOutput.getPaths();

        // get path lengths
        double[] pathLengths = new double[pathList.size()];
        int idx = 0;
        Iterator it = pathList.iterator();
        while (it.hasNext()) {
            int[] path = (int[]) it.next();
            pathLengths[idx++] = (double) path.length;
        }

        // averages
        double path_avg_1 = 0, path_avg_2 = 0;
        double len_avg_1 = 0, len_avg_2 = 0;
        for (int i = 0; i < iter - lag; i++) {
            path_avg_1 += pathProbs[i].doubleValue();
            len_avg_1 += pathLengths[i];

            path_avg_2 += pathProbs[i + lag].doubleValue();
            len_avg_2 += pathLengths[i + lag];
        }
        int nSample = iter - lag;
        path_avg_1 /= (double) nSample;
        path_avg_2 /= (double) nSample;
        len_avg_1 /= (double) nSample;
        len_avg_2 /= (double) nSample;

        double path_A = 0, path_B1 = 0, path_B2 = 0;
        double len_A = 0, len_B1 = 0, len_B2 = 0;
        for (int i = 0; i < iter - lag; i++) {
            path_A += (pathProbs[i].doubleValue() - path_avg_1) * (pathProbs[i + lag].doubleValue() - path_avg_2);
            len_A += (pathLengths[i] - len_avg_1) * (pathLengths[i + lag] - len_avg_2);

            path_B1 += Math.pow(pathProbs[i].doubleValue() - path_avg_1, 2);
            len_B1 += Math.pow(pathLengths[i] - len_avg_1, 2);
            path_B2 += Math.pow(pathProbs[i + lag].doubleValue() - path_avg_2, 2);
            len_B2 += Math.pow(pathLengths[i + lag] - len_avg_2, 2);
        }

        double corr_path = path_A / Math.sqrt(path_B1 * path_B2);
        double corr_len = len_A / Math.sqrt(len_B1 * len_B2);

        return new double[]{corr_path, corr_len};
    }

    private double getAverageNeighboursCount(PhyloTree tree, double edgeCount) {

        return getAverageNeighboursCount(tree.getRoot(), edgeCount) / (double) tree.getLeaves().size();

    }

    private double getAverageNeighboursCount(PhyloNode phylonode, double edgeCount) {
        double avgNC = 0.0;
        if (!phylonode.isLeaf()) {
            // go to left sub tree
            avgNC += getAverageNeighboursCount(phylonode.getLeftSon(), edgeCount);
            // go to right sub tree
            avgNC += getAverageNeighboursCount(phylonode.getRightSon(), edgeCount);
        } else { // leaf
            MetabolicNetwork theNetwork = phylonode.getMetabolicNetwork();
            avgNC = (double) Utilities.sum(theNetwork.getNeighboursCount()) / edgeCount;
        }

        return avgNC;
    }
}
