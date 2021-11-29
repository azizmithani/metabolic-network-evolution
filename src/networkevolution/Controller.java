/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package networkevolution;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import network.DirectedReaction;
import network.MetabolicNetwork;
import network.Metabolite;
import network.NetworkObject;
import network.Reaction;
import rahnumadatabase.Database;
import rahnumautilities.Email;
import rahnumautilities.Utilities;
import systemobject.Organism;
import systemobject.PhyloNode;
import systemobject.PhyloTree;

/**
 *
 * @author mithani
 */
public class Controller {

    static public String DB_DRIVER = "com.mysql.jdbc.Driver";
    // test environment
    static public String DB_URL = "jdbc:mysql://localhost:3306/";
    static public String DB_NAME = "rahnumatest";
    static public String DB_USER = "rahnumatestdba";
    static public String DB_PASSWORD = "oveaQue5";
//    // production environment
//    static public String DB_URL = "jdbc:mysql://peng.stats.ox.ac.uk:3306/";
//    static public String DB_NAME = "rahnuma";
//    static public String DB_USER = "rahnumadba";
//    static public String DB_PASSWORD = "eetha5Ku";

    // 
    static private boolean saveData = false;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        boolean external = false;

        if (external) {
            runExternalCommand(args);
            Email email = new Email("markov.stats.ox.ac.uk", "25");
            email.postMail(new String[]{"mithani@stats.ox.ac.uk"}, "Job Completed", Utilities.toString(args, " "), "mithani@stats.ox.ac.uk", null);
            System.exit(0);
        }

//        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
//        database.readApplicationParameters();

        String[] orgIds1, orgIds2, pathwayCodes;
//        orgIds1 = new String[]{"pae","pau","ppu","ppf","pst","psb","psp","pfl","pfo","pen","pmy"};

        orgIds1 = new String[]{"pae", "pau", "ppu", "pst", "psb", "psp", "pfl", "pfo", "pen"};
//        createDistanceMatrix(orgIds1);

        int evolutionModel = NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT;
        orgIds1 = new String[]{"pae"};
        orgIds2 = new String[]{"pst"};
        String[] orgIdsForCoreAndProhib = new String[]{"pae", "pau", "ppu", "pst", "psb", "psp", "pfl", "pfo", "pen"};
        pathwayCodes = new String[]{"00030"};
        //runGibbsSampler(evolutionModel, orgIds1, orgIds2, pathwayCodes, orgIdsForCoreAndProhib);
        orgIds1 = new String[]{"pfl"};
        orgIds2 = new String[]{"pfo"};
        //runGibbsSampler(evolutionModel, orgIds1, orgIds2, pathwayCodes, orgIdsForCoreAndProhib);
        //networkProperties(orgIdsForCoreAndProhib, pathwayCodes);
        //getNeighboursCount(orgIdsForCoreAndProhib, pathwayCodes, true);

//        printMetabolicNetworks(new String[]{"psb", "psp", "pst"}, new String[]{"00340"}, orgIdsForCoreAndProhib);
//        printMetabolicNetworks(new String[]{"pae", "pfl", "pfo", "pst"}, new String[]{"00340"}, orgIdsForCoreAndProhib);
//        printMetabolicNetworks(new String[]{"pfl", "pfo", "psb", "psp", "pst"}, new String[]{"00010", "00020", "00030"}, orgIdsForCoreAndProhib);
//        printMetabolicNetworks(new String[]{"pae", "pau", "pen", "ppu", "pfl", "pfo", "psb", "psp", "pst"}, pathwayCodes, orgIdsForCoreAndProhib, true);

        pathwayCodes = new String[]{"00251", "00252", "00260", "00271", "00272", "00280", "00290", "00300", "00310", "00330", "00340", "00350", "00360", "00380", "00400", "00220", "00010", "00020", "00030", "00620", "00630", "00640", "00910", "00410", "00430", "00440", "00450", "00460", "00471", "00472", "00473", "00480"};
        //getNeighboursCount(orgIds2, pathwayCodes, true);
        //pathwayCodes = new String[]{"00040","00051","00052","00053","00500","00530","00520","00650","00660","00031","00680","00920","00061","00071","00100","00561","00564","00230","00240","00730","00740","00750","00760","00770","00780","00790","00670","00860","00130"};
        //pathwayCodes = new String[]{"00562"};
        //printMetabolicNetworks(new String[]{"pae", "pau", "pen", "ppu", "pfl", "pfo", "psb", "psp", "pst"}, pathwayCodes, orgIdsForCoreAndProhib, true);
        for (int p = 0; p < pathwayCodes.length; p++) {
            //printMetabolicNetworks(new String[]{"pae", "pau", "pen", "ppu", "pfl", "pfo", "psb", "psp", "pst"}, new String[]{pathwayCodes[p]}, orgIdsForCoreAndProhib, true);
            //printReferenceReactionsAndStatus(new String[]{pathwayCodes[p]}, orgIdsForCoreAndProhib, true);
            String str = "nohup java -server -jar NetworkEvolution.jar phylo 2 . ref_" + pathwayCodes[p] + ".txt pae_" + pathwayCodes[p] + ".txt,pau_" + pathwayCodes[p] + ".txt,pen_" + pathwayCodes[p] + ".txt,ppu_" + pathwayCodes[p] + ".txt,pfl_" + pathwayCodes[p] + ".txt,pfo_" + pathwayCodes[p] + ".txt,psb_" + pathwayCodes[p] + ".txt,psp_" + pathwayCodes[p] + ".txt,pst_" + pathwayCodes[p] + ".txt \"(((ppu:0.09054,pen:0.07823):0.06648,((pfo:0.09127,pfl:0.08092):0.01186,((psp:0.05003,psb:0.04650):0.01840,pst:0.07138):0.07489):0.04414):0.13239,(pau:0.00563,pae:0.00494):0.13239)\"" + " core_" + pathwayCodes[p] + ".txt prohibited_" + pathwayCodes[p] + ".txt -1 -1 3 60000 10000 10 -1 > PhyloGibbs_all_" + pathwayCodes[p] + "_br.out &";

             //System.out.println(str);
             //System.out.println();
        }

        orgIds1 = new String[]{"pfl"};
        orgIds2 = new String[]{"pfo"};
        pathwayCodes = new String[]{"00360"};
//        runGibbsSampler(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, orgIds1, orgIds2, pathwayCodes, orgIdsForCoreAndProhib);

/*        if (args != null && args.length > 0) {
            orgIds1 = new String[]{args[0]};
            orgIds2 = new String[]{args[1]};
            pathwayCodes = new String[]{args[2]};
            if (args.length > 3 && args[3].startsWith("T"))
                orgIdsForCoreAndProhib = new String[]{"pae", "pau", "ppu", "pst", "psb", "psp", "pfl", "pfo", "pen"};
            else
                orgIdsForCoreAndProhib = null;
        }
*/
 //        runGibbsSampler(evolutionModel, orgIds1, orgIds2, pathwayCodes, orgIdsForCoreAndProhib);

//        networkProperties(new String[]{"pst"}, pathwayCodes);
//        networkProperties(new String[]{"pfo"}, pathwayCodes);
//        networkProperties(new String[]{"pfl"}, pathwayCodes);
//        networkProperties(new String[]{"pst"}, pathwayCodes);

        //pathwayCodes = new String[]{"00340"}; // histidine
        //runGibbsSampler(NetworkEvolver.EVOL_MODEL_HYBRID, orgIds1, orgIds2, pathwayCodes, orgIdsForCoreAndProhib);
        //runGibbsSampler(NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE, orgIds1, orgIds2, pathwayCodes);

        //studyKOCascading(evolutionModel, orgIds1, pathwayCodes);

        //printPathProbabilities(pathwayCodes, "20085526_105502_0.5_0.3model_2.out");
        //printSummaries(orgIds1, pathwayCodes, "20080604_191210_0.5_0.3_model_2.out");

//        String[] orgIdsCore, orgIdsProhib;
//        orgIdsCore = new String[]{"pae","pau","ppu","ppf","pst","psb","psp","pfl","pfo","pen","pmy"};
//        orgIdsProhib = Utilities.extractColumn(database.getOrganisms("Gamma/others"), 0);
//        NetworkSimulator.simulateEvolution(orgIds1, null, orgIdsCore, orgIdsProhib);
        //NetworkSimulator.simulateEvolution(null, pathwayCodes, null, null);

        // Build the phylogenetic tree

        String dirName;
        try {
            dirName = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            dirName = null;
        }
        String refNWFileName = "toy_network_G.txt";
        String startNWFileName = "toy_network_H1.txt";
        String endNWFileName = "toy_network_H2.txt";
//        simulateNetworkEvolution(NetworkEvolver.EVOL_MODEL_HYBRID, dirName, refNWFileName, startNWFileName, null, null, 0.5,0.3, 1.0);

        refNWFileName = "toy_network_A.txt";
        String thirdNWFileName = "toy_network_A.txt";
        double insRate = 0.1;
        double delRate = 0.2;
        equilibriumProbability(evolutionModel, dirName,refNWFileName, thirdNWFileName, insRate, delRate, 2, 1.0);
        equilibriumProbability(evolutionModel, dirName,refNWFileName, thirdNWFileName, insRate, delRate, 4, 1.0);
        equilibriumProbability(evolutionModel, dirName,refNWFileName, thirdNWFileName, insRate, delRate, 6, 1.0);
        equilibriumProbability(evolutionModel, dirName,refNWFileName, thirdNWFileName, insRate, delRate, 8, 1.0);
        equilibriumProbability(evolutionModel, dirName,refNWFileName, thirdNWFileName, insRate, delRate, 10, 1.0);
//        testPathSamplerForPathLengths(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName);
//        testGibbsSampler(evolutionModel, dirName, refNWFileName, startNWFileName);
//        runGibbsSampler(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, null, null, -1, -1, -1);
//        generateRateMatrix_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, 0.514, 0.131, 0.0);
        double[] dp = new double[]{0.1,0.3,0.5,0.7,0.9};
        for (int d = 0; d < dp.length; d++) {
          //generateRateMatrices_LH_Surface_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, dp[d]);
        }
        for (int i = 0; i < 10; i++) {
        //transitionProbability(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, null, null, 0.1, 0.1, 1, i+1);
        }
//        transitionProbability(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, null, null, 0.5, 0.3, 1, 7);

        String nw3Filename = "toy_network_H3.txt";
        String[] nwFilenames = new String[]{startNWFileName, endNWFileName, nw3Filename};
        evolutionModel = NetworkEvolver.EVOL_MODEL_HYBRID;
//        runGibbsSamplerOnPhylogeny(evolutionModel, dirName, refNWFileName, nwFilenames, -1.0);
        //testGibbsSamplerOnPhylogeny(evolutionModel, dirName, refNWFileName, nwFilenames, 1.0);
        //runGibbsSamplerOnPhylogeny_LH_Surface(evolutionModel, dirName, refNWFileName, nwFilenames);
        //generateRateMatrices1_LH_Surface_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName);
        //calculate_LH_Surface_Phylo(evolutionModel, dirName, refNWFileName, nwFilenames);

        refNWFileName = "ref_00310.txt";
        startNWFileName = "pfo_00310.txt";
        evolutionModel = NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT;
        double dependenceProbability = 0.25;
        //getNeighboursCount(dirName, refNWFileName, null, true);
        //simulateNetworkEvolution(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, dirName, refNWFileName, startNWFileName, null, null, 0.6,0.3, dependenceProbability);
//        calculatedExpectedNumberOfEvents(evolutionModel,dirName, refNWFileName, startNWFileName, 0.0625, 0.5);
//        calculatedExpectedNumberOfEvents(evolutionModel,dirName, refNWFileName, startNWFileName, 0.125, 0.5);
//        calculatedExpectedNumberOfEvents(evolutionModel,dirName, refNWFileName, startNWFileName, 0.25, 0.5);
//        calculatedExpectedNumberOfEvents(evolutionModel,dirName, refNWFileName, startNWFileName, 0.5, 0.5);
//        calculatedExpectedNumberOfEvents(evolutionModel,dirName, refNWFileName, startNWFileName, 1, 0.5);
//        calculatedExpectedNumberOfEvents(evolutionModel,dirName, refNWFileName, startNWFileName, 2, 0.5);
//        calculatedExpectedNumberOfEvents(evolutionModel,dirName, refNWFileName, startNWFileName, 4, 0.5);


        refNWFileName = "ref_00030.txt";
        endNWFileName = "pfo_00030.txt";
        nw3Filename = "pae_00030.txt";
        nwFilenames = new String[]{startNWFileName, endNWFileName, nw3Filename};
        String coreFilename = "core_00030.txt";
        String prohibFilename = "prohibited_00030.txt";
        dependenceProbability = -1.0;
        //getNeighboursCount(dirName, refNWFileName, endNWFileName, true);
//        simulateNetworkEvolution(NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE, dirName, refNWFileName,
//                endNWFileName, coreFilename, prohibFilename, 1.6802, 1.9563, 1, 60000, 10, 0.0);
//        simulateNetworkEvolution(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, dirName, refNWFileName,
//                endNWFileName, coreFilename, prohibFilename, 1.38941, 0.21218, 1, 60000, 10, 1.0);
//        System.out.println("Hyrbid Model");
//        simulateNetworkEvolution(NetworkEvolver.EVOL_MODEL_HYBRID, dirName, refNWFileName,
//                endNWFileName, coreFilename, prohibFilename, 0.7444, 0.5353, 1, 60000, 10, 0.4608);
//        runGibbsSamplerOnPhylogeny(evolutionModel, dirName, refNWFileName, nwFilenames, coreFilename, prohibFilename, dependenceProbability);
//        getNeighboursCount_ToyNetwork(dirName, refNWFileName, null, true);
//        getNeighboursCount_ToyNetwork(dirName, refNWFileName, startNWFileName, true);
//        simulateNetworkEvolution_ToyNetworks(NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE, dirName, refNWFileName, startNWFileName);
//        emAlgorithm_ToyNetworks_SimulatedData(NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE, dirName, refNWFileName, startNWFileName);

//        runGibbsSampler_ToyNetworks(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, dirName, refNWFileName, startNWFileName, endNWFileName);


//        generateRateMatrix_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, 2.97, 1.02);
//        generateRateMatrix_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, 0.15, 0.3);
//        generateRateMatrix_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, 0.6,  0.3);
//        generateRateMatrix_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, 1.2,  0.3);
//        generateRateMatrix_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName, 2.4,  0.3);

    // MCMC to calculation of P(H1->H2) & P(H2->H1) as well as generate rate matrices for different rate combinations
//       testPathSampler_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName);


//        testEquilibriumProbability(evolutionModel, dirName, refNWFileName, startNWFileName);
//        testEquilibriumProbability(evolutionModel, dirName, refNWFileName, endNWFileName);
//        testEquilibriumProbability_SubNetworkSize(evolutionModel, dirName, refNWFileName, startNWFileName);
//        testEquilibriumProbability_SubNetworkSize(evolutionModel, dirName, refNWFileName, endNWFileName);
//        String[] endNWFileNames = new String[]{"toy_network_H2_1.txt", "toy_network_H2_2.txt", "toy_network_H2_3.txt", "toy_network_H2_4.txt", "toy_network_H2_5.txt"};
//        for (int en=0; en < endNWFileNames.length; en++) {
//            testEquilibriumProbability_SubNetworkSize_ToyNetworks(evolutionModel, dirName, refNWFileName, endNWFileNames[en]);
//        }


//        generateRateMatrices_LH_Surface_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName);

    //        simulateNetworkEvolution_ToyNetworks(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, dirName, refNWFileName, startNWFileName);
//        runGibbsSampler_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName);
//        runGibbsSampler_LH_Surface_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName);
//        equilibriumProbability_ToyNetworks_surface(evolutionModel, dirName, refNWFileName, startNWFileName);
//        String[] endNWFileNames = new String[]{"toy_network_H2_1.txt", "toy_network_H2_2.txt", "toy_network_H2_3.txt", "toy_network_H2_4.txt", "toy_network_H2_5.txt"};
//        profileGibbsSampler_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileNames);
//        for (int en=0; en < endNWFileNames.length; en++) {
//            generateRateMatrix_ToyNetworks(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, dirName, refNWFileName, startNWFileName, endNWFileNames[en], 0.7, 0.9);
//        }

//    runGibbsSampler_ToyNetworks(evolutionModel, dirName, refNWFileName, startNWFileName, endNWFileName);
    //runMCMC_ToyNetworks_SimulatedData(evolutionModel, dirName, refNWFileName, startNWFileName);

//        runMCMCOnPhylogeny_ToyNetworks(evolutionModel, dirName, refNWFileName, nwFilenames);

//        exponentiation_ToyNetworks(NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, dirName, refNWFileName, startNWFileName, endNWFileName);
//        equilibriumProbability_ToyNetworks(evolutionModel, dirName, refNWFileName, endNWFileName);


//        printSummaries(orgIds1, pathwayCodes, "Data for path length distribution/MCMC_toyNetwork20080624_200114_0.05_0.03_model_2.out");
//        printSummaries(orgIds1, pathwayCodes, "Data for path length distribution/MCMC_toyNetwork20080624_203114_0.05_0.03_model_2.out");
//        printSummaries(orgIds1, pathwayCodes, "Data for path length distribution/MCMC_toyNetwork20080625_051517_0.05_0.03_model_2.out");
//        printSummaries(orgIds1, pathwayCodes, "Data for path length distribution/MCMC_toyNetwork20080625_051517_0.05_0.03_model_2.out");
//        printSummaries(orgIds1, pathwayCodes, "Data for path length distribution/MCMC_toyNetwork20080625_103720_0.05_0.03_model_2.out");

//        evolutionModel = NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT;
//        studyNetworkCharacteristics(evolutionModel, dirName, refNWFileName, refNWFileName);



//        try {
//            File path = new File(".");
//            String[] list = path.list();
//            for (int i = 0; i < list.length; i++) {
//                String dataFilename = list[i];
//
//                if (dataFilename.endsWith(".out")) {
//                    ObjectInputStream in =
//                            new ObjectInputStream(
//                            new FileInputStream(new File(dirName, dataFilename)));
//
//                    ArrayList pathList = (ArrayList) in.readObject(); // pathList
//                    ArrayList networkList = (ArrayList) in.readObject(); // networkList
//                    BigDecimal[] pathProbList = (BigDecimal[]) in.readObject(); // pathProbList
//                    Double[] insRates = (Double[]) in.readObject(); // insRates
//                    Double[] delRates = (Double[]) in.readObject(); // delRates
//                    int acceptCount = in.readInt();
//
//                    int nBurning = 10000;
//                    HashMap insRateDistribution = NetworkEvolver.rateDistribution(insRates, nBurning, true, "Insertion Rate:");
//                    HashMap delRateDistribution = NetworkEvolver.rateDistribution(delRates, nBurning, true, "Deletion Rate:");
//                    NetworkEvolver.pathLengthDistribution(pathList, nBurning, true);
//                    HashMap likelihood = new NetworkMCMC().calculateLikelihood(pathList, pathProbList, insRates, delRates, nBurning, true);
//
//                    in.close();
//                } // end if
//            } // for each file
//        } catch (Exception e) {
//            e.printStackTrace();
//        }

    }

    static public void runExternalCommand(String[] args) {
        if (args != null && args.length > 0) {
            try {
                if (args[0].equalsIgnoreCase("simulate")) {
                    if (args.length < 2 || args[1].equalsIgnoreCase("help")) {
                        System.out.println("Usage: NetworkEvolution.jar simulate <job_parameters>");
                        System.out.println();
                        System.out.println("Job Parameters:");
                        System.out.println("\t1. Evolution Model - 1: Independent Edge Model, 2: Neighbour-Dependent Model");
                        System.out.println("\t2. Directory containing network and other files (. for current directory)");
                        System.out.println("\t3. Reference Network Filename");
                        System.out.println("\t4. Start Network Filename");
                        System.out.println("\t5. Core Hyperedges Filename (\"\" for none)");
                        System.out.println("\t6. Prohibited Hyperedges Filename (\"\" for none)");
                        System.out.println("\t7. Insertion Rate");
                        System.out.println("\t8. Deletion Rate");
                        System.out.println("\t9. Number of Simulation Runs");
                        System.out.println("\t10. Number of Iterations per Run");
                        System.out.println("\t11. Update Interval for Output Display");
                        System.out.println("\t12. Neighbour Dependence Probability");

                    } else {
                        int evolutionModel = Integer.parseInt(args[1]);
                        String dirName = (args[2].equals(".") ? new File(".").getCanonicalPath() : args[2]);
                        String refNWFilename = args[3];
                        String startNWFilename = args[4];
                        String coreFilename = args[5];
                        String prohibFilename = args[6];
                        double insRate = Double.parseDouble(args[7]);
                        double delRate = Double.parseDouble(args[8]);
                        int nRun = Integer.parseInt(args[9]);
                        int nIter = Integer.parseInt(args[10]);
                        int updateInterval = Integer.parseInt(args[11]);
                        double dependenceProbability = Double.parseDouble(args[12]);

                        simulateNetworkEvolution(evolutionModel, dirName, refNWFilename, startNWFilename,
                                coreFilename, prohibFilename, insRate, delRate, nRun, nIter, updateInterval, dependenceProbability);
                    }

                } else if (args[0].equalsIgnoreCase("gibbs")) {
                    if (args.length < 2 || args[1].equalsIgnoreCase("help")) {
                        System.out.println("Usage: NetworkEvolution.jar gibbs <job_parameters>");
                        System.out.println();
                        System.out.println("Job Parameters:");
                        System.out.println("\t1. Evolution Model - 1: Independent Edge Model, 2: Neighbour-Dependent Model");
                        System.out.println("\t2. Directory containing network and other files (. for current directory)");
                        System.out.println("\t3. Reference Network Filename");
                        System.out.println("\t4. Start Network Filename");
                        System.out.println("\t5. End Network Filename");
                        System.out.println("\t6. Core Hyperedges Filename (\"\" for none)");
                        System.out.println("\t7. Prohibited Hyperedges Filename (\"\" for none)");
                        System.out.println("\t8. Insertion Rate (-ve value for random starting value)");
                        System.out.println("\t9. Deletion Rate (-ve value for random starting value)");
                        System.out.println("\t10. Number of Runs");
                        System.out.println("\t11. Number of Iterations per Run");
                        System.out.println("\t12. Number of Iterations for Burn-in Period");
                        System.out.println("\t13. Update Interval for Output Display");
                        System.out.println("\t14. Neighbour Dependence Probability (-ve value for random starting value)");

                    } else {
                        int evolutionModel = Integer.parseInt(args[1]);
                        String dirName = (args[2].equals(".") ? new File(".").getCanonicalPath().replace("\\", "/") : args[2]);
                        String refNWFilename = args[3];
                        String startNWFilename = args[4];
                        String endNWFilename = args[5];
                        String coreFilename = args[6];
                        String prohibFilename = args[7];
                        double insRate = Double.parseDouble(args[8]);
                        double delRate = Double.parseDouble(args[9]);
                        int nRun = Integer.parseInt(args[10]);
                        int nIter = Integer.parseInt(args[11]);
                        int nBurning = Integer.parseInt(args[12]);
                        int updateInterval = Integer.parseInt(args[13]);
                        double dependenceProbability = Double.parseDouble(args[14]);

                        runGibbsSampler(evolutionModel, dirName, refNWFilename, startNWFilename,
                                endNWFilename, coreFilename, prohibFilename, insRate, delRate,
                                nRun, nIter, nBurning, updateInterval, dependenceProbability);
                    }

                } else if (args[0].equalsIgnoreCase("path")) {
                    if (args.length < 2 || args[1].equalsIgnoreCase("help")) {
                        System.out.println("Usage: NetworkEvolution.jar path <job_parameters>");
                        System.out.println();
                        System.out.println("Job Parameters:");
                        System.out.println("\t1. Evolution Model - 1: Independent Edge Model, 2: Neighbour-Dependent Model");
                        System.out.println("\t2. Directory containing network and other files (. for current directory)");
                        System.out.println("\t3. Reference Network Filename");
                        System.out.println("\t4. Start Network Filename");
                        System.out.println("\t5. End Network Filename");
                        System.out.println("\t6. Core Hyperedges Filename (\"\" for none)");
                        System.out.println("\t7. Prohibited Hyperedges Filename (\"\" for none)");
                        System.out.println("\t8. Insertion Rate");
                        System.out.println("\t9. Deletion Rate");
                        System.out.println("\t10. Number of Runs");
                        System.out.println("\t11. Number of Iterations per Run");
                        System.out.println("\t12. Number of Iterations for Burn-in Period");
                        System.out.println("\t13. Update Interval for Output Display");
                        System.out.println("\t14. Neighbour Dependence Probability");

                    } else {
                        int evolutionModel = Integer.parseInt(args[1]);
                        String dirName = (args[2].equals(".") ? new File(".").getCanonicalPath().replace("\\", "/") : args[2]);
                        String refNWFilename = args[3];
                        String startNWFilename = args[4];
                        String endNWFilename = args[5];
                        String coreFilename = args[6];
                        String prohibFilename = args[7];
                        double insRate = Double.parseDouble(args[8]);
                        double delRate = Double.parseDouble(args[9]);
                        int nRun = Integer.parseInt(args[10]);
                        int nIter = Integer.parseInt(args[11]);
                        int nBurning = Integer.parseInt(args[12]);
                        int updateInterval = Integer.parseInt(args[13]);
                        double dependenceProbability = Double.parseDouble(args[14]);

                        runPathSampler(evolutionModel, dirName, refNWFilename, startNWFilename,
                                endNWFilename, coreFilename, prohibFilename, insRate, delRate,
                                nRun, nIter, nBurning, updateInterval, dependenceProbability);
                    }

                } else if (args[0].equalsIgnoreCase("phylo")) {
                    if (args.length < 2 || args[1].equalsIgnoreCase("help")) {
                        System.out.println("Usage: NetworkEvolution.jar phylo <job_parameters>");
                        System.out.println();
                        System.out.println("Job Parameters:");
                        System.out.println("\t1. Evolution Model - 1: Independent Edge Model, 2: Neighbour-Dependent Model");
                        System.out.println("\t2. Directory containing network and other files (. for current directory)");
                        System.out.println("\t3. Reference Network Filename");
                        System.out.println("\t4. Network Filenames (seperated by comma)");
                        System.out.println("\t5. Phylogeny");
                        System.out.println("\t6. Core Hyperedges Filename (\"\" for none)");
                        System.out.println("\t7. Prohibited Hyperedges Filename (\"\" for none)");
                        System.out.println("\t8. Insertion Rate (-ve value for random starting value)");
                        System.out.println("\t9. Deletion Rate (-ve value for random starting value)");
                        System.out.println("\t10. Number of Simulation Runs");
                        System.out.println("\t11. Number of Iterations per Run");
                        System.out.println("\t12. Number of Iterations for Burn-in Period");
                        System.out.println("\t13. Update Interval for Output Display");
                        System.out.println("\t14. Neighbour Dependence Probability (-ve value for random starting value)");

                    } else {
                        int evolutionModel = Integer.parseInt(args[1]);
                        String dirName = (args[2].equals(".") ? new File(".").getCanonicalPath().replace("\\", "/") : args[2]);
                        String refNWFilename = args[3];
                        String[] nwFilenames = Utilities.toArray(args[4], ",");
                        String phylogeny = args[5];
                        String coreFilename = args[6];
                        String prohibFilename = args[7];
                        double insRate = Double.parseDouble(args[8]);
                        double delRate = Double.parseDouble(args[9]);
                        int nRun = Integer.parseInt(args[10]);
                        int nIter = Integer.parseInt(args[11]);
                        int nBurning = Integer.parseInt(args[12]);
                        int updateInterval = Integer.parseInt(args[13]);
                        double dependenceProbability = Double.parseDouble(args[14]);

                        runGibbsSamplerOnPhylogeny(evolutionModel, dirName, refNWFilename,
                                nwFilenames, phylogeny, coreFilename, prohibFilename,
                                insRate, delRate, nRun, nIter, nBurning, updateInterval, dependenceProbability);
                    }
                } else if (args[0].equalsIgnoreCase("help")) {
                    System.out.println("Usage: NetworkEvolution.jar <job> <job_parameters>");
                    System.out.println();
                    System.out.println("Jobs:");
                    System.out.println("\tsimulate - Simulate Network Evolution");
                    System.out.println("\tpath - MCMC for Path Sampling");
                    System.out.println("\tgibbs - Gibbs Sampler for Parameter Estimation");
                    System.out.println();
                    System.out.println("Type NetworkEvolution.jar <job> help to see parameters for the job");


                } else {
                    System.out.println("Invalid arguments.Type NetworkEvolution.jar help to see valid arguments");
                }
            } catch (Exception ex) {
                ex.printStackTrace();
                System.out.println("Invalid arguments. Type NetworkEvolution.jar help to see valid arguments");
            }
        } else {
            System.out.println("Invalid arguments. Type NetworkEvolution.jar help to see valid arguments");
        }

    }

    /* ********************************************************************** *
     *                  D I S T A N C E     M A T R I X                       *
     * ********************************************************************** */
    static private double[][] createDistanceMatrix(String[] orgIds, String ignoreMetabolitesFile,
            boolean filterNetwork, boolean enzymeBased) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        String ignoreDir;
        ArrayList ignoreMetabolites;
        try {
            ignoreDir = new File(".").getCanonicalPath();
            if (ignoreMetabolitesFile != null && !ignoreMetabolitesFile.isEmpty())
                ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreMetabolitesFile, false));
            else
                ignoreMetabolites = null;
        } catch (Exception ex) {
            ignoreDir = null;
            ignoreMetabolites = null;
        }

        ArrayList lstNetworks = new ArrayList();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways("");
        String[] pathwayIds = Utilities.extractColumn(keggPathways, 0);

        for (int i = 0; i < orgIds.length; i++) {
            String[][] arrOrg = database.getOrganisms(new String[]{orgIds[i]});
            Organism org = new Organism(arrOrg[0][0], arrOrg[0][1]);
            ArrayList orgList = new ArrayList();
            orgList.add(org);

            String[][] orgReactions = database.getOrganismReactions(new String[]{org.getID()}, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(new String[]{org.getID()}, pathwayIds);

            MetabolicNetwork orgNetwork = MetabolicNetwork.buildMetabolicNetwork(org.getID(),
                    org.getName(), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            // filter network if required
            if (filterNetwork) {
                int neighbourhoodCutoff = 1;
                orgNetwork = orgNetwork.filterNetworkByNeighbourhood(neighbourhoodCutoff, ignoreMetabolites);
            }

            // add to the list of networks
            lstNetworks.add(orgNetwork);
        } // end for all organisms

        // create distance matrix based on reaction or enzyme data (as requested)
        double[][] distMatrix;
        if (enzymeBased)
            distMatrix = MetabolicNetwork.createDistanceMatrixByEnzymes(lstNetworks);
        else
            distMatrix = MetabolicNetwork.createDistanceMatrix(lstNetworks);

        // print the matrix
        for (int i = 0; i < distMatrix.length; i++) {
            double[] ds = distMatrix[i];
            System.out.print(((MetabolicNetwork) lstNetworks.get(i)).getID() + "\t");
            for (int j = 0; j < ds.length; j++) {
                double d = ds[j];
                System.out.print(d);
                System.out.print("\t");
            }
            System.out.println();
        }

        return distMatrix;
    }

    /* ********************************************************************** *
     *                        S I M U L A T I O N                             *
     * ********************************************************************** */
    static private void simulateNetworkEvolution(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String coreFilename,
            String prohibFilename, double insRate, double delRate, double modelSwitchProbability) {

        // Simulation parameters
        int nIter = 100000;
        int updateInterval = 0;
        int nRun = 3;

        simulateNetworkEvolution(evolutionModel, dirName, refNWFilename, startNWFilename,
                coreFilename, prohibFilename, insRate, delRate, nRun, nIter, updateInterval, saveData, modelSwitchProbability);
    }

    // external function
    static public void simulateNetworkEvolution(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String coreFilename,
            String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int updateInterval, double modelSwitchProbability) {

        // call the main function with false for save data
        simulateNetworkEvolution(evolutionModel, dirName, refNWFilename, startNWFilename,
                coreFilename, prohibFilename, insRate, delRate, nRun, nIter, updateInterval, false, modelSwitchProbability);

    }

    // Main Function
    static private void simulateNetworkEvolution(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String coreFilename,
            String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int updateInterval, boolean save_datafile, double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkSequence();
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // core edges
        ArrayList coreEdges = new ArrayList();
        if (coreFilename != null && coreFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, coreFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                coreEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        // prohibited edges
        ArrayList prohibEdges = new ArrayList();
        if (prohibFilename != null && prohibFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, prohibFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                prohibEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        
        // variables to hold summary data 
        int[][] insFreq = new int[startNetwork.getDirectedReactions().size()][nRun]; // number of insertions for each edge
        int[][] delFreq = new int[startNetwork.getDirectedReactions().size()][nRun]; // number of deletions for each edge

        NetworkSimulator networkSimulator = new NetworkSimulator();
        for (int rn = 0; rn < nRun; rn++) {
            System.out.println("Run: " + Integer.toString(rn + 1));
            // call the simulation function
            ArrayList lstNetworks = networkSimulator.simulateNetworkEvolution(startNetwork, evolutionModel, nIter, coreEdges, prohibEdges, insRate, delRate, updateInterval, dependenceProbability);
            // variables to hold summary data 
            int[] insFreqRun = new int[startNetwork.getDirectedReactions().size()]; // number of insertions for each edge
            int[] delFreqRun = new int[startNetwork.getDirectedReactions().size()]; // number of deletions for each edge

            NetworkEvolver.summariseEvolution(lstNetworks, insFreqRun, delFreqRun, null, null);
            for (int i = 0; i < insFreqRun.length; i++) {
                insFreq[i][rn] = insFreqRun[i];
                delFreq[i][rn] = delFreqRun[i];
            }

            if (save_datafile) {
                try {
                    String filename = "simulation_";
                    filename += "_" + Double.toString(insRate) + "_" + Double.toString(delRate);
                    filename += "_model_" + Integer.toString(evolutionModel);
                    filename += "_run_" + Integer.toString(rn + 1) + ".jdat";
                    ObjectOutputStream out =
                            new ObjectOutputStream(
                            new FileOutputStream(filename));
                    out.writeObject(startNetwork.getReactionSequence());
                    out.writeObject(lstNetworks);
                    out.writeObject(insFreq);
                    out.writeObject(delFreq);
                    out.close(); // Also flushes output
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

        }// for rn = 1 to nRun


        // print the output
        System.out.println("Insertion Frequency:");
        System.out.print("ID" + "\t");
        for (int i = 1; i <= nRun; i++)
            System.out.print(i + "\t");
        System.out.println();
        for (int i = 0; i < insFreq.length; i++) {
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < insFreq[i].length; j++) {
                System.out.print(insFreq[i][j]);
                System.out.print("\t");
            }
            System.out.println();
        }

        System.out.println("Deletion Frequency:");
        System.out.print("ID" + "\t");
        for (int i = 1; i <= nRun; i++)
            System.out.print(i + "\t");
        System.out.println();
        for (int i = 0; i < delFreq.length; i++) {
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < delFreq[i].length; j++) {
                System.out.print(delFreq[i][j]);
                System.out.print("\t");
            }
            System.out.println();
        }

    }


    static private void simulateNetworkEvolutionNeighboursCount(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String coreFilename,
            String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int updateInterval, double dependenceProbability) {
            // calls the function that prints neighbours count at each iteration instead of networks

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkSequence();
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // core edges
        ArrayList coreEdges = new ArrayList();
        if (coreFilename != null && coreFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, coreFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                coreEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        // prohibited edges
        ArrayList prohibEdges = new ArrayList();
        if (prohibFilename != null && prohibFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, prohibFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                prohibEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }


        // variables to hold summary data
        int[][] insFreq = new int[startNetwork.getDirectedReactions().size()][nRun]; // number of insertions for each edge
        int[][] delFreq = new int[startNetwork.getDirectedReactions().size()][nRun]; // number of deletions for each edge

        NetworkSimulator networkSimulator = new NetworkSimulator();
        for (int rn = 0; rn < nRun; rn++) {
            System.out.println("Run: " + Integer.toString(rn + 1));
            // call the simulation function
            ArrayList lstNetworks = networkSimulator.simulateNetworkEvolutionNeighboursCount(startNetwork, evolutionModel, nIter, 0, coreEdges, prohibEdges, insRate, delRate, updateInterval, dependenceProbability);
            // variables to hold summary data
            int[] insFreqRun = new int[startNetwork.getDirectedReactions().size()]; // number of insertions for each edge
            int[] delFreqRun = new int[startNetwork.getDirectedReactions().size()]; // number of deletions for each edge

            NetworkEvolver.summariseEvolution(lstNetworks, insFreqRun, delFreqRun, null, null);
            for (int i = 0; i < insFreqRun.length; i++) {
                insFreq[i][rn] = insFreqRun[i];
                delFreq[i][rn] = delFreqRun[i];
            }

        }// for rn = 1 to nRun


        // print the output
        System.out.println("Insertion Frequency:");
        System.out.print("ID" + "\t");
        for (int i = 1; i <= nRun; i++)
            System.out.print(i + "\t");
        System.out.println();
        for (int i = 0; i < insFreq.length; i++) {
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < insFreq[i].length; j++) {
                System.out.print(insFreq[i][j]);
                System.out.print("\t");
            }
            System.out.println();
        }

        System.out.println("Deletion Frequency:");
        System.out.print("ID" + "\t");
        for (int i = 1; i <= nRun; i++)
            System.out.print(i + "\t");
        System.out.println();
        for (int i = 0; i < delFreq.length; i++) {
            System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
            for (int j = 0; j < delFreq[i].length; j++) {
                System.out.print(delFreq[i][j]);
                System.out.print("\t");
            }
            System.out.println();
        }

    }

    /* ********************************************************************** *
     *                    G I B B S    S A M P L E R                          *
     * ********************************************************************** */
    static private void runGibbsSampler(int evolutionModel, String[] orgIds1, String[] orgIds2,
            String[] pathwayIds, String[] orgIdsForCoreAndProhib) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null) {
            pathwayIds = Utilities.extractColumn(keggPathways, 0);
        }

        String ignoreDir;
        try {
            ignoreDir = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            ignoreDir = null;
        }
        String ignoreFile = "ignoreMetabolites.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);
        refNetwork.setupNetworkSequence();
//        refNetwork.calculateReactionConnectivity(true);

        // create the metabolic networks
        MetabolicNetwork orgNetwork1, orgNetwork2;
        // network 1
        if (orgIds1 == null || orgIds1.length == 0) {
            orgNetwork1 = refNetwork;
        } else {
            String[][] arrOrg = database.getOrganisms(orgIds1);
            ArrayList orgList = Organism.buildOrganismList(arrOrg);

            String[][] orgReactions = database.getOrganismReactions(orgIds1, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIds1, pathwayIds);

            orgNetwork1 = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
                    Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            orgNetwork1.addInactiveReactions(refNetwork.getReactions());
            orgNetwork1.removeInactiveMetabolites();
            orgNetwork1.getMetabolites().removeAll(ignoreMetabolites);


            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE) {
                orgNetwork1.setupNetworkMatrices();
            }

            orgNetwork1.setupNetworkSequence();
        }

        // network 2
        if (orgIds2 == null || orgIds2.length == 0) {
            orgNetwork2 = refNetwork.clone();
        } else {
            String[][] arrOrg = database.getOrganisms(orgIds2);
            ArrayList orgList = Organism.buildOrganismList(arrOrg);

            String[][] orgReactions = database.getOrganismReactions(orgIds2, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIds2, pathwayIds);

            orgNetwork2 = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
                    Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

            orgNetwork2.addInactiveReactions(refNetwork.getReactions());
            orgNetwork2.removeInactiveMetabolites();
            orgNetwork2.getMetabolites().removeAll(ignoreMetabolites);

            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE) {
                orgNetwork2.setupNetworkMatrices();
            }

            orgNetwork2.setupNetworkSequence();
        }

        // get core and prohibted edges
        ArrayList coreEdges, prohibEdges;
        MetabolicNetwork unionNetwork;
        if (orgIdsForCoreAndProhib == null || orgIdsForCoreAndProhib.length == 0) {
            coreEdges = new ArrayList();
            prohibEdges = new ArrayList();
            unionNetwork = refNetwork;
        } else {
            String[][] arrOrg = database.getOrganisms(orgIdsForCoreAndProhib);
            ArrayList orgList = Organism.buildOrganismList(orgIdsForCoreAndProhib);

            String[][] orgReactions = database.getOrganismReactions(orgIdsForCoreAndProhib, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIdsForCoreAndProhib, pathwayIds);

            // CORE EDGES
            MetabolicNetwork coreNetwork = MetabolicNetwork.buildMetabolicNetwork("core",
                    "Core Network", orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);
            // delete the reactions not present in all organisms
            coreNetwork.deleteReactionsNotPresentInAllOrganisms();
            // setup hyperedges
            coreNetwork.setupNetworkSequence();
            // get the core edges
            coreEdges = coreNetwork.getDirectedReactions();

            // PROHIBITED EDGES
            unionNetwork = MetabolicNetwork.buildMetabolicNetwork("prohib",
                    "Prohibited Network", orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);
            // setup hyperedges
            unionNetwork.setupNetworkSequence();
            // get the prohibited edges
            prohibEdges = refNetwork.getDirectedReactions();
            prohibEdges.removeAll(unionNetwork.getDirectedReactions());

        }

        // MCMC parameters
        int nIter = 60000, nBurning = 10000;
        int updateInterval = 10;
        int nRun = 3;
        double evolTime = 1;
        boolean enumStates = false;
        NetworkMCMC networkMCMC = new NetworkMCMC();

        BigDecimal insRate = new BigDecimal("1.0");
        BigDecimal delRate = new BigDecimal("2.0");
        BigDecimal dependenceProbability = new BigDecimal("0.5");
        System.out.println(new Date());
        System.out.println("Start: " + orgNetwork1.toString() + " End: " + orgNetwork2.toString() + " Pathway: " + Utilities.toString(pathwayIds));

        for (int rn = 0; rn < nRun; rn++) {
            System.out.println("Run: " + Integer.toString(rn + 1));

            // call the MCMC function
            boolean sampleRates = true;
            boolean sampleDP = true;

            if (sampleRates) {
                // generate random rates
                insRate = new BigDecimal(Double.toString(Math.random()));
                delRate = new BigDecimal(Double.toString(Math.random()));
            } else {
                insRate = new BigDecimal("1.0");
                delRate = new BigDecimal("2.0");
            }

            if (sampleDP)
                dependenceProbability = new BigDecimal(Double.toString(Math.random()));
            else
                dependenceProbability = new BigDecimal("1.0");
            // calll the MCMC function
            Date startTime = new Date();
            MCMCOutput mcmcOutput = networkMCMC.networkMCMC(orgNetwork1, orgNetwork2, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), evolTime, nIter, nBurning, enumStates, updateInterval, evolutionModel, sampleRates, dependenceProbability.doubleValue(), sampleDP);
            Date endTime = new Date();

            ArrayList pathList = mcmcOutput.getPaths();
            ArrayList networkList = mcmcOutput.getVisitedNetworks();
            BigDecimal[] pathProbList = mcmcOutput.getPathProbabilities();
            Double[] insRates = mcmcOutput.getInsRates();
            Double[] delRates = mcmcOutput.getDelRates();
            Double[] dependenceProbabilities = mcmcOutput.getDependenceProbabilities();
            int acceptCountPath = mcmcOutput.getAcceptCountPath();
            int acceptCountParameter = mcmcOutput.getAcceptCountParameter();
            int acceptCountDependenceProbability = mcmcOutput.getAcceptCountDependenceProbability();
            //System.out.println(new Date());
            //NetworkEvolver.summariseEvolution(orgNetwork1, networkList);
            System.out.println(new Date());
            //HashMap pathDistribution = NetworkEvolver.pathDistribution(pathList, nBurning, true);
            //HashMap networkDistribution = NetworkEvolver.networkDistribution(orgNetwork1, new ArrayList(pathDistribution.keySet()), true);
            HashMap insRateDistribution = NetworkEvolver.rateDistribution(insRates, nBurning, true, "Insertion Rate:");
            HashMap delRateDistribution = NetworkEvolver.rateDistribution(delRates, nBurning, true, "Deletion Rate:");
            HashMap dpDistribution = NetworkEvolver.rateDistribution(dependenceProbabilities, nBurning, true, "Dependence Probabilities:");
            HashMap likelihood = networkMCMC.calculateLikelihoodDP(pathList, pathProbList, insRates, delRates, dependenceProbabilities, nBurning, true);
            System.out.println("Acceptance Count - Path: " + Integer.toString(acceptCountPath));
            System.out.println("Acceptance Count - Parameter: " + Integer.toString(acceptCountParameter));
            System.out.println("Acceptance Count - Dependence Probability: " + Integer.toString(acceptCountDependenceProbability));
            System.out.println("Total Time: " + Long.toString(endTime.getTime() - startTime.getTime()) + " milliseconds");
            System.out.println();
           // networkMCMC.calculateAutocorrelation(mcmcOutput, true);

            if (saveData) {
                try {
                    String filename = "Gibbs_" + orgNetwork1.getID() + "_" + orgNetwork2.getID() + "_" + Utilities.toString(pathwayIds, "");
                    filename += "_Model_" + Integer.toString(evolutionModel);
                    filename += (coreEdges.size() == 0 && prohibEdges.size() == 0 ? "" : "_CnP");
                    filename += "_Run_" + Integer.toString(rn) + ".jdat";
                    ObjectOutputStream out =
                            new ObjectOutputStream(
                            new FileOutputStream(new File(filename)));
                    out.writeInt(nIter);
                    out.writeInt(nBurning);
                    out.writeInt(evolutionModel);
                    out.writeDouble(evolTime);
                    out.writeObject(orgNetwork1.getReactionSequence());
                    out.writeObject(orgNetwork2.getReactionSequence());
                    out.writeObject(pathList);
                    out.writeObject(networkList);
                    out.writeObject(pathProbList);
                    out.writeObject(insRates);
                    out.writeObject(delRates);
                    out.writeInt(acceptCountPath);
                    out.writeInt(acceptCountParameter);
                    out.writeObject(dependenceProbabilities);
                    out.writeInt(acceptCountDependenceProbability);
                    out.close(); // Also flushes output
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }// for rn = 1to nRun
    }

    // external function
    static public void runGibbsSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename,
            String coreFilename, String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int nBurning, int updateInterval, double dependenceProbability) {

        // call the main function with false for save data
        runGibbsSampler(evolutionModel, dirName, refNWFilename, startNWFilename, endNWFilename,
                coreFilename, prohibFilename, insRate, delRate, nRun, nIter, nBurning,
                updateInterval, false, dependenceProbability);
    }

    static private void runGibbsSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename,
            String coreFilename, String prohibFilename, double insRate, double delRate,
            double dependenceProbability) {

        // MCMC parameters
        int nIter = 110000, nBurning = 10000;
        int nRun = 1;
        int updateInterval = 10;

        runGibbsSampler(evolutionModel, dirName, refNWFilename, startNWFilename, endNWFilename,
                coreFilename, prohibFilename, insRate, delRate, nRun, nIter, nBurning,
                updateInterval, saveData, dependenceProbability);
    }

    // Main Function
    static private void runGibbsSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename,
            String coreFilename, String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int nBurning, int updateInterval, boolean save_datafile, 
            double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkSequence();

        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        // MCMC parameters
        double evolTime = 1.0;
        boolean enumStates = false;
        NetworkMCMC networkMCMC = new NetworkMCMC();

        // core edges
        ArrayList coreEdges = new ArrayList();
        if (coreFilename != null && coreFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, coreFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                coreEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        // prohibited edges
        ArrayList prohibEdges = new ArrayList();
        if (prohibFilename != null && prohibFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, prohibFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                prohibEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        int startIdx, endIdx;
        try {
            // start and end network indices
            String strStartNetwork = Utilities.toString(startNetwork.getReactionSequence(), "");
            String strEndNetwork = Utilities.toString(endNetwork.getReactionSequence(), "");
            startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
            endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);
        } catch (Exception ex) {
            startIdx = 0;
            endIdx = 0;
        }

        for (int rn = 0; rn < nRun; rn++) {
            System.out.println(new Date());
            System.out.println("Run: " + Integer.toString(rn + 1));

            boolean sampleRates = true;
            boolean sampleDP = false;
            if (evolutionModel == NetworkEvolver.EVOL_MODEL_HYBRID)
                sampleDP = true;
            else if (evolutionModel == NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT)
                dependenceProbability = 1.0;

            if (sampleRates) {
                // generate random rates
                if (insRate < 0)
                    insRate = Math.random();
                if (delRate < 0)
                    delRate = Math.random();
            }
            if (sampleDP) {
                // generate random probability value
                if (dependenceProbability < 0)
                    dependenceProbability = Math.random();                
            }            
            
            // calll the MCMC function
            Date startTime = new Date();
            MCMCOutput mcmcOutput = networkMCMC.networkMCMC(startNetwork, endNetwork, coreEdges, prohibEdges, insRate, delRate, evolTime, nIter, nBurning, enumStates, updateInterval, evolutionModel, sampleRates, dependenceProbability, sampleDP);
            Date endTime = new Date();
            ArrayList pathList = mcmcOutput.getPaths();
            ArrayList networkList = mcmcOutput.getVisitedNetworks();
            BigDecimal[] pathProbList = mcmcOutput.getPathProbabilities();
            Double[] insRates = mcmcOutput.getInsRates();
            Double[] delRates = mcmcOutput.getDelRates();
            Double[] dependenceProbabilities = mcmcOutput.getDependenceProbabilities();
            int acceptCountPath = mcmcOutput.getAcceptCountPath();
            int acceptCountParameter = mcmcOutput.getAcceptCountParameter();
            int acceptCountDependenceProbability = mcmcOutput.getAcceptCountDependenceProbability();
            //NetworkEvolver.summariseEvolution(startNetwork, networkList);
            //System.out.println(new Date());
            //HashMap networkDistribution = NetworkEvolver.networkDistribution(startNetwork, new ArrayList(pathDistribution.keySet()), true);
            NetworkEvolver.pathDistribution(pathList, nBurning, true);
            NetworkEvolver.rateDistribution(insRates, nBurning, true, "Insertion Rate:");
            NetworkEvolver.rateDistribution(delRates, nBurning, true, "Deletion Rate:");
            NetworkEvolver.pathLengthDistribution(pathList, nBurning, true);
            networkMCMC.calculateLikelihoodDP(pathList, pathProbList, insRates, delRates, dependenceProbabilities, nBurning, true);
            //networkMCMC.calculateAutocorrelation(mcmcOutput, true);
            System.out.println("Acceptance Count - Path: " + Integer.toString(acceptCountPath));
            System.out.println("Acceptance Count - Parameter: " + Integer.toString(acceptCountParameter));
            System.out.println("Acceptance Count - Dependence Probability: " + Integer.toString(acceptCountDependenceProbability));
            System.out.println("Total Time: " + Long.toString(endTime.getTime() - startTime.getTime()) + " milliseconds");
            System.out.println();

            if (save_datafile) {

                try {
                    String filename = "Gibbs_toyNetwork";// + dateFormat.format(new Date());
                    filename += "_from_" + Integer.toString(startIdx) + "_to_" + Integer.toString(endIdx);
                    filename += "_model_" + Integer.toString(evolutionModel);
                    filename += "_run_" + Integer.toString(rn + 1) + ".jdat";
                    ObjectOutputStream out =
                            new ObjectOutputStream(
                            new FileOutputStream(filename));
                    out.writeInt(nIter);
                    out.writeInt(nBurning);
                    out.writeInt(evolutionModel);
                    out.writeDouble(evolTime);
                    out.writeInt(startIdx);
                    out.writeObject(startNetwork.getReactionSequence());
                    out.writeInt(endIdx);
                    out.writeObject(endNetwork.getReactionSequence());
                    out.writeObject(pathList);
                    out.writeObject(networkList);
                    out.writeObject(pathProbList);
                    out.writeObject(insRates);
                    out.writeObject(delRates);
                    out.writeInt(acceptCountPath);
                    out.writeInt(acceptCountParameter);
                    out.writeObject(dependenceProbabilities);
                    out.writeInt(acceptCountDependenceProbability);
                    out.close(); // Also flushes output
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }// end for rn=1 to nRun
    }

    static private void runGibbsSampler_LH_Surface(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename, double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        // MCMC parameters
        int nIter = 110000, nBurning = 10000;
        int updateInterval = 1000;
        double evolTime = 1;
        boolean enumStates = false;
        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();
        NetworkMCMC networkMCMC = new NetworkMCMC();

        //String saveDirName = "/users/mithani";
        String saveDirName;
        try {
            saveDirName = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            saveDirName = "";
        }

        BigDecimal insRate = new BigDecimal("0.50");
        BigDecimal stepSizeIR = new BigDecimal("0.10");
        BigDecimal stepSizeDR = new BigDecimal("0.05");
        int nDim = 20;
        BigDecimal[][] LH_MCMC = new BigDecimal[nDim][nDim];

        try {
            String outFilename = "LH_MCMC_toy_model" + Integer.toString(evolutionModel) + ".txt";
            // create the output file
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(saveDirName, outFilename)));

            for (int ir = 0; ir < nDim; ir++) {
                BigDecimal delRate = new BigDecimal("0.30");
                for (int dr = 0; dr < nDim; dr++) {
                    System.out.println(new Date());
                    System.out.println("Insertion Rate: " + insRate.toString() + ", Deletion Rate: " + delRate.toString());

                    // call the MCMC function
                    MCMCOutput mcmcOutput = networkMCMC.networkMCMC(startNetwork, endNetwork, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), evolTime, nIter, nBurning, enumStates, updateInterval, evolutionModel, false, dependenceProbability,false);
                    ArrayList pathList = mcmcOutput.getPaths();
                    ArrayList networkList = mcmcOutput.getVisitedNetworks();
                    BigDecimal[] pathProbList = mcmcOutput.getPathProbabilities();
                    Double[] insRates = mcmcOutput.getInsRates();
                    Double[] delRates = mcmcOutput.getDelRates();
                    int acceptCountPath = mcmcOutput.getAcceptCountPath();
                    int acceptCountParameter = mcmcOutput.getAcceptCountParameter();
                    HashMap likelihood = networkMCMC.calculateLikelihood(pathList, pathProbList, insRates, delRates, nBurning, true);

                    LH_MCMC[ir][dr] = (BigDecimal) likelihood.get(new ArrayList(likelihood.keySet()).get(0));

                    if (saveData) {
                        try {
                            String filename = "Gibbs_toyNetwork";
                            filename += "_" + insRate.toString() + "_" + delRate.toString();
                            filename += "_model_" + Integer.toString(evolutionModel) + ".jdat";
                            ObjectOutputStream outData =
                                    new ObjectOutputStream(
                                    new FileOutputStream(new File(saveDirName, filename)));
                            outData.writeObject(pathList);
                            outData.writeObject(networkList);
                            outData.writeObject(pathProbList);
                            outData.writeObject(insRates);
                            outData.writeObject(delRates);
                            outData.writeInt(acceptCountPath);
                            outData.writeInt(acceptCountParameter);
                            outData.close(); // Also flushes output
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }

                    delRate = delRate.add(stepSizeDR);
                }

                // write the likelihood values for current insertion rate
                BigDecimal[] ds = LH_MCMC[ir];
                for (int j = 0; j < ds.length - 1; j++) {
                    out.write(ds[j].round(new MathContext(16)).toString() + "\t");
                }
                // write the last element
                out.write(ds[ds.length - 1].round(new MathContext(16)).toString());
                out.write("\r\n");
                out.flush();


                insRate = insRate.add(stepSizeIR);

            } // end for each insertion rate 

            out.close(); // Also flushes output
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println(Utilities.toString(LH_MCMC, "\t"));

    }

    static private void calculate_LH_Surface_Phylo(int evolutionModel, String dirName,
            String refNWFilename, String[] nwFilenames) {

        String phylogeny = "((H1:1.0,H2:1.0):1.0,H3:1.0)";
        // load the metabolic networks
        ArrayList lstNetworks = new ArrayList(nwFilenames.length);
        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        for (int i = 0; i < nwFilenames.length; i++) {
            // network 1
            MetabolicNetwork theNetwork = loadNetwork(dirName, nwFilenames[i]);
            theNetwork.addInactiveReactions(refNetwork.getReactions());
            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                theNetwork.setupNetworkMatrices();
            theNetwork.setupNetworkSequence();

            String strNetwork = Utilities.toString(theNetwork.getReactionSequence(), "");
            int nwIdx = Integer.parseInt(new StringBuffer(strNetwork).reverse().toString(), 2);

            System.out.println("Network Idx - " + theNetwork.getID() + ": " + Integer.toString(nwIdx));

            // add to the list of networks
            lstNetworks.add(theNetwork);
        }

        ArrayList lstOrganism = Organism.buildOrganismList(PhyloTree.extractNodesFromPhylogeny(phylogeny));
        PhyloTree phylotree = PhyloTree.buildPhylogeneticTree(phylogeny, lstOrganism);
        // Build Metabolic network on leaves
        assignNetworksOnTheLeaves(phylotree.getRoot(), lstNetworks);

        // MCMC parameters
        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();

        String saveDirName;
        try {
            saveDirName = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            saveDirName = "";
        }

        // get indices of core and prohibited edges
        Iterator itCore = coreEdges.iterator();
        Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iCoreEdges.length; i++)
            iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        while (itCore.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itCore.next();
            iCoreEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }
        Iterator itProhib = prohibEdges.iterator();
        Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iProhibEdges.length; i++)
            iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        while (itProhib.hasNext()) {
            DirectedReaction reaction = (DirectedReaction) itProhib.next();
            iProhibEdges[refNetwork.getDirectedReactions().indexOf(reaction)] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
        }

        int edgeCount = refNetwork.getDirectedReactions().size();
        // array to store unalterable edges
        Byte[] unalterableEdges = new Byte[edgeCount];
        // define core and prohibited edges as unalterable edges.
        // all the rest are alterable edges.
        for (int i = 0; i < edgeCount; i++) {
            if (iCoreEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            else if (iProhibEdges[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT))
                unalterableEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
            else
                unalterableEdges[i] = -1;
        } // end for all edges

        // enumerate all possible networks
        NetworkEvolver nwEvolver = new NetworkEvolver();
        ArrayList lstAllNetworks = nwEvolver.enumerateAllNetworks(edgeCount, unalterableEdges);

        BigDecimal insRate = new BigDecimal("0.40");
        BigDecimal stepSizeIR = new BigDecimal("0.10");
        BigDecimal stepSizeDR = new BigDecimal("0.04");
        int nDim = 20;
        BigDecimal[][] LH_Full = new BigDecimal[nDim][nDim];
        try {
            String outFilename = "LH_Full_Phylo.txt";
            // create the output file
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(saveDirName, outFilename)));

            for (int ir = 0; ir < nDim; ir++) {
                BigDecimal delRate = new BigDecimal("0.10");
                for (int dr = 0; dr < nDim; dr++) {
                    System.out.println(new Date());
                    System.out.println("Insertion Rate: " + insRate.toString() + ", Deletion Rate: " + delRate.toString());

                    String filename = "toy_rateMatrix_" + insRate.toString() + "_" + delRate.toString() + ".txt";
                    String[][] strRM = Utilities.readList2(saveDirName, filename, false);
                    double[][] rateMatrix = new double[strRM.length][strRM[0].length];
                    for (int i = 0; i < strRM.length; i++) {
                        for (int j = 0; j < strRM[0].length; j++) {
                            rateMatrix[i][j] = Double.parseDouble(strRM[i][j]);
                        }
                    }
                    // calculate the transition probability matrix
                    double[][] transProbs = nwEvolver.calculateTransitionProbabilityMatrix(rateMatrix, 1.0);
                    // calculate the equilibrium probabilities
                    double[] eqProbs = nwEvolver.calculateEquilibriumProbabilityVector(rateMatrix);

//                    // read the transition probability matrix
//                    String filename = "toy_TransProb_" + insRate.toString() + "_" + delRate.toString() + ".txt";
//                    String[][] strTP = Utilities.readList2(saveDirName, filename, false);
//                    double[][] transProbs = new double[strTP.length][strTP[0].length];
//                    for (int i = 0; i < strTP.length; i++) {
//                        for (int j = 0; j < strTP[0].length; j++) {
//                            transProbs[i][j] = Double.parseDouble(strTP[i][j]);
//                        }
//                    }
//                    
//                    // calculate the equilibrium probabilities
//                    filename = "toy_EqProb_" + insRate.toString() + "_" + delRate.toString() + ".txt";
//                    String[] strEq = Utilities.readList(saveDirName, filename, false);
//                    double[] eqProbs = new double[strEq.length];
//                    for (int i = 0; i < strEq.length; i++) {
//                            eqProbs[i] = Double.parseDouble(strEq[i]);
//                    }

                    // caclulate the likelihood
                    double likelihood = nwEvolver.calculateFullLikelihood(phylotree.getRoot(), refNetwork, lstAllNetworks, null, insRate.doubleValue(), delRate.doubleValue(), transProbs, eqProbs);

                    // save
                    try {
                        LH_Full[ir][dr] = new BigDecimal(likelihood);
                    } catch (Exception ex) {
                        LH_Full[ir][dr] = new BigDecimal("0.0");
                    }
                    System.out.println("Likelihood: " + LH_Full[ir][dr].round(new MathContext(16)).toString());

                    delRate = delRate.add(stepSizeDR);
                }

                // write the likelihood values for current insertion rate
                BigDecimal[] ds = LH_Full[ir];
                for (int j = 0; j < ds.length - 1; j++) {
                    out.write(ds[j].round(new MathContext(16)).toString() + "\t");
                }
                // write the last element
                out.write(ds[ds.length - 1].round(new MathContext(16)).toString());
                out.write("\r\n");
                out.flush();


                insRate = insRate.add(stepSizeIR);

            } // end for each insertion rate 

            out.close(); // Also flushes output
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println(Utilities.toString(LH_Full, "\t"));

    }

    // external function
    static private void runPathSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename,
            String coreFilename, String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int nBurning, int updateInterval, double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkSequence();
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        // core edges
        ArrayList coreEdges = new ArrayList();
        if (coreFilename != null && coreFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, coreFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                coreEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        // prohibited edges
        ArrayList prohibEdges = new ArrayList();
        if (prohibFilename != null && prohibFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, prohibFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                prohibEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        // call the main function with false for save data
        runPathSampler(evolutionModel, startNetwork, endNetwork, coreEdges, prohibEdges,
                insRate, delRate, nRun, nIter, nBurning, updateInterval, false, dependenceProbability);

    }

    static private void runPathSampler(int evolutionModel,
            MetabolicNetwork startNetwork, MetabolicNetwork endNetwork,
            ArrayList coreEdges, ArrayList prohibEdges,
            double insRate, double delRate, double dependenceProbability) {

        // MCMC parameters
        int nIter = 110000, nBurning = 10000;
        int nRun = 3;
        int updateInterval = 0;

        runPathSampler(evolutionModel, startNetwork, endNetwork, coreEdges, prohibEdges,
                insRate, delRate, nRun, nIter, nBurning, updateInterval, saveData, dependenceProbability);
    }

    // Main function
    static private void runPathSampler(int evolutionModel,
            MetabolicNetwork startNetwork, MetabolicNetwork endNetwork,
            ArrayList coreEdges, ArrayList prohibEdges, double insRate, double delRate,
            int nRun, int nIter, int nBurning, int updateInterval, boolean save_datafile, 
            double dependenceProbability) {

        int startIdx, endIdx;
        try {
            // start and end network indices
            String strStartNetwork = Utilities.toString(startNetwork.getReactionSequence(), "");
            String strEndNetwork = Utilities.toString(endNetwork.getReactionSequence(), "");
            startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
            endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);
        } catch (Exception ex) {
            startIdx = 0;
            endIdx = 0;
        }

        double evolTime = 1.0;
        boolean enumStates = false;
        NetworkMCMC networkMCMC = new NetworkMCMC();

        //System.out.println(new Date());
        BigDecimal[] likelihood = new BigDecimal[nRun];
        for (int rn = 0; rn < nRun; rn++) {
            //System.out.println("Run: " + Integer.toString(rn + 1));

            // call the MCMC function
            boolean sampleRates = false;

            MCMCOutput mcmcOutput = networkMCMC.networkMCMC(startNetwork, endNetwork, coreEdges, prohibEdges, insRate, delRate, evolTime, nIter, nBurning, enumStates, updateInterval, evolutionModel, sampleRates, dependenceProbability,false);
            ArrayList pathList = mcmcOutput.getPaths();
            ArrayList networkList = mcmcOutput.getVisitedNetworks();
            BigDecimal[] pathProbList = mcmcOutput.getPathProbabilities();
            Double[] insRates = mcmcOutput.getInsRates();
            Double[] delRates = mcmcOutput.getDelRates();
            int acceptCountPath = mcmcOutput.getAcceptCountPath();
            //System.out.println(new Date());
            NetworkEvolver.pathLengthDistribution(pathList, nBurning, true);
            HashMap lh = networkMCMC.calculateLikelihood(pathList, pathProbList, insRates, delRates, nBurning, false);

            likelihood[rn] = (BigDecimal) new ArrayList(lh.values()).get(0);

            if (save_datafile) {

                try {
                    String filename = "PathMCMC_toyNetwork";
                    filename += "_from_" + Integer.toString(startIdx) + "_to_" + Integer.toString(endIdx);
                    filename += "_rates_" + Double.toString(insRate) + "_" + Double.toString(delRate);
                    filename += "_run_" + Integer.toString(rn + 1) + ".jdat";
                    ObjectOutputStream out =
                            new ObjectOutputStream(
                            new FileOutputStream(filename));
                    out.writeInt(nIter);
                    out.writeInt(nBurning);
                    out.writeInt(evolutionModel);
                    out.writeDouble(evolTime);
                    out.writeDouble(insRate);
                    out.writeDouble(delRate);
                    out.writeInt(startIdx);
                    out.writeObject(startNetwork.getReactionSequence());
                    out.writeInt(endIdx);
                    out.writeObject(endNetwork.getReactionSequence());
                    out.writeObject(pathList);
                    out.writeObject(networkList);
                    out.writeObject(pathProbList);
                    out.writeObject(insRates);
                    out.writeObject(delRates);
                    out.writeInt(acceptCountPath);
                    out.close(); // Also flushes output
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }// end for rn=1 to nRun
        System.out.print(Integer.toString(startIdx) + "\t" + Integer.toString(endIdx));
        System.out.print("\tlambda = " + Double.toString(insRate) + ", mu = " + Double.toString(delRate));
        for (int i = 0; i < likelihood.length; i++) {
            System.out.print("\t" + likelihood[i].round(new MathContext(16)));
        }
        System.out.println();
    }

    static private void profileGibbsSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String[] endNWFilenames, double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // MCMC parameters
        int nBurning = 1000;
        int nRun = 3;
        int updateInterval = 5000;
        double evolTime = 1.0;
        int[] iterations = new int[]{5000, 10000, 20000, 40000, 80000};
        boolean enumStates = false;
        boolean sampleRates = false;
        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();
        NetworkMCMC networkMCMC = new NetworkMCMC();
        // rates
        BigDecimal insRate = new BigDecimal("0.7");
        BigDecimal delRate = new BigDecimal("0.9");
        System.out.println("Insertion Rate: " + insRate.toString() + ", Deletion Rate: " + delRate.toString());

        // variable to store profiling results
        BigDecimal[][][] profile_likelihood = new BigDecimal[endNWFilenames.length][iterations.length][nRun];
        int[][][] profile_accept_count = new int[endNWFilenames.length][iterations.length][nRun];
        long[][][] profile_cpu_time = new long[endNWFilenames.length][iterations.length][nRun];
        int[] differences = new int[endNWFilenames.length];


        for (int en = 0; en < endNWFilenames.length; en++) {

            // end network
            String endNWFilename = endNWFilenames[en];
            MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
            endNetwork.addInactiveReactions(refNetwork.getReactions());
            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                endNetwork.setupNetworkMatrices();
            endNetwork.setupNetworkSequence();
            int nDiff = Integer.parseInt(endNWFilename.substring(endNWFilename.lastIndexOf("_") + 1).replace(".txt", ""));
            System.out.println("Differences: " + Integer.toString(nDiff));
            //store the difference
            differences[en] = nDiff;

            // number of iterations
            for (int it = 0; it < iterations.length; it++) {
                int nIter = iterations[it] + nBurning;
                System.out.println("Iterations: " + Integer.toString(nIter));

                // need to average over multiple runs
                for (int rn = 0; rn < nRun; rn++) {
                    System.out.println("Run: " + Integer.toString(rn + 1));

                    // calll the MCMC function
                    Date startTime = new Date();
                    // call the MCMC function
                    MCMCOutput mcmcOutput = networkMCMC.networkMCMC(startNetwork, endNetwork,
                            coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(),
                            evolTime, nIter, nBurning, enumStates, updateInterval, evolutionModel,
                            sampleRates, dependenceProbability, false);
                    Date endTime = new Date();
                    Long cpuTime = endTime.getTime() - startTime.getTime();

                    ArrayList pathList = mcmcOutput.getPaths();
                    ArrayList networkList = mcmcOutput.getVisitedNetworks();
                    BigDecimal[] pathProbList = mcmcOutput.getPathProbabilities();
                    Double[] insRates = mcmcOutput.getInsRates();
                    Double[] delRates = mcmcOutput.getDelRates();
                    int acceptCountPath = mcmcOutput.getAcceptCountPath();
                    int acceptCountParameter = mcmcOutput.getAcceptCountParameter();
                    HashMap likelihood = networkMCMC.calculateLikelihood(pathList, pathProbList, insRates, delRates, nBurning, true);
                    System.out.println("Acceptance Count - Path: " + Integer.toString(acceptCountPath));
                    System.out.println("Total Time: " + Long.toString(cpuTime) + " milliseconds");

                    // store the profile results
                    profile_likelihood[en][it][rn] = (BigDecimal) new ArrayList(likelihood.values()).get(0);
                    profile_accept_count[en][it][rn] = acceptCountPath;
                    profile_cpu_time[en][it][rn] = cpuTime;

                    if (saveData) {
                        try {
                            String filename = "Gibbs_profile_toy_model_" + Integer.toString(evolutionModel);
                            filename += "_d_" + Integer.toString(nDiff);
                            filename += "_it_" + Integer.toString(iterations[it]);
                            filename += "_rn_" + Integer.toString(rn + 1);
                            filename += ".jdat";
                            ObjectOutputStream out =
                                    new ObjectOutputStream(
                                    new FileOutputStream(filename));
                            out.writeObject(pathList);
                            out.writeObject(networkList);
                            out.writeObject(pathProbList);
                            out.writeObject(insRates);
                            out.writeObject(delRates);
                            out.writeInt(acceptCountPath);
                            out.writeInt(acceptCountParameter);
                            out.close(); // Also flushes output
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }// end for rn=1 to nRun

            }

        }

        // print the likelihood
        System.out.println("Likelihood: ");
        for (int i = 0; i < profile_likelihood.length; i++) { // differences
            String diff = Integer.toString(differences[i]);
            for (int j = 0; j < profile_likelihood[i].length; j++) { // iterations
                String iter = "\t" + Integer.toString(iterations[j]);
                System.out.print(diff + iter);
                for (int k = 0; k < profile_likelihood[i][j].length; k++) {
                    System.out.print("\t" + profile_likelihood[i][j][k].round(new MathContext(16)));
                }
                System.out.println();
            }
            System.out.println();
        }

        // print the acceptance percentage
        System.out.println("Acceptance Percentage: ");
        for (int i = 0; i < profile_accept_count.length; i++) { // differences
            String diff = Integer.toString(differences[i]);
            for (int j = 0; j < profile_accept_count[i].length; j++) { // iterations
                String iter = "\t" + Integer.toString(iterations[j]);
                System.out.print(diff + iter);
                for (int k = 0; k < profile_accept_count[i][j].length; k++) {
                    System.out.print("\t" + Double.toString((double) profile_accept_count[i][j][k] / (double) iterations[j]));
                }
                System.out.println();
            }
            System.out.println();
        }

        // print the cpu time
        System.out.println("CPU time: ");
        for (int i = 0; i < profile_cpu_time.length; i++) { // differences
            String diff = Integer.toString(differences[i]);
            for (int j = 0; j < profile_cpu_time[i].length; j++) { // iterations
                String iter = "\t" + Integer.toString(iterations[j]);
                System.out.print(diff + iter);
                for (int k = 0; k < profile_cpu_time[i][j].length; k++) {
                    System.out.print("\t" + Long.toString(profile_cpu_time[i][j][k]));
                }
                System.out.println();
            }
            System.out.println();
        }

    }

    static private void testGibbsSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, double dependenceProbability) {


        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // core and prohibited edges
        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();
        // Rates
        BigDecimal insRate = new BigDecimal("0.06");
        BigDecimal delRate = new BigDecimal("0.03");

        int nBurningSim = 50000;
//        int nIterSim = 100000;
        int nIterSim = 50000;
        ArrayList lstSimulation = new NetworkSimulator().simulateNetworkEvolutionWithTime(startNetwork, evolutionModel, nIterSim, nBurningSim, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), dependenceProbability);
        ArrayList lstSimulatedNetworks = (ArrayList) lstSimulation.get(0);
        ArrayList lstSimulatedTimes = (ArrayList) lstSimulation.get(1);
        Object[] arrSimulatedTimes = lstSimulatedTimes.toArray();

        try {
            String resultSaveDirName = new File(".").getCanonicalPath();
            String resultfilename = "parameters_simulatedData.txt";
            // create the output file
            BufferedWriter outResult = new BufferedWriter(new FileWriter(new File(resultSaveDirName, resultfilename)));

            Iterator itTimes = lstSimulatedTimes.iterator();
            Iterator itNetworks = lstSimulatedNetworks.iterator();
            int idx = 0;
            while (itNetworks.hasNext()) {
                outResult.write(Integer.toString(idx++) + ":\t" + (Double) itTimes.next() + "\t" + Utilities.toString((Byte[]) itNetworks.next(), " ") + "\r\n");
            }
            outResult.flush();

            System.out.println("Estimating parameters");
            System.out.println("Sample\t" + "lambda\t" + "mu\t" + "time\t" + "ratio\t" + "lambda*time\t" + "mu*time");
            outResult.write("Sample\t" + "lambda\t" + "mu\t" + "time\t" + "ratio\t" + "lambda*time\t" + "mu*time" + "\r\n");

            int stepSize = 10;
            //Double[] arrRatio = new Double[nIterSim/stepSize];
            //Iterator itTimes = lstSimulatedTimes.iterator();
            //Iterator itNetworks = lstSimulatedNetworks.iterator();
            for (int i = -1; i < nIterSim - 1; i += stepSize) {
                //System.out.println("Sample: " + Integer.toString(i+1) + " to " + Integer.toString(i + stepSize));

                // get the start network and time
                //Double t1 = (Double)itTimes.next();
                //startNetwork = refNetwork.getNetwork((Byte[]) itNetworks.next());
                startNetwork = refNetwork.getNetwork((Byte[]) lstSimulatedNetworks.get(i + 1));
                startNetwork.addInactiveReactions(refNetwork.getReactions());
                if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                    startNetwork.setupNetworkMatrices();
                startNetwork.setupNetworkSequence();

                // get the end network and time
                //Double t2 = (Double)itTimes.next();
                //MetabolicNetwork endNetwork = refNetwork.getNetwork((Byte[]) itNetworks.next());
                MetabolicNetwork endNetwork = refNetwork.getNetwork((Byte[]) lstSimulatedNetworks.get(i + stepSize));
                endNetwork.addInactiveReactions(refNetwork.getReactions());
                if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                    endNetwork.setupNetworkMatrices();
                endNetwork.setupNetworkSequence();

                //double evolTime = t2 - t1;
                double evolTime = (Double) arrSimulatedTimes[i + stepSize] - (Double) arrSimulatedTimes[i + 1];
                double evolTimeMCMC = 1.0; //evolTime;


                // MCMC parameters
                int nIter = 11000, nBurning = 1000;
                int nRun = 1;
                int updateInterval = 0;
                boolean enumStates = false;
                NetworkMCMC networkMCMC = new NetworkMCMC();

                //System.out.println(new Date());
                //System.out.println("Insertion Rate: " + insRate.toString() + ", Deletion Rate: " + delRate.toString());

                for (int rn = 0; rn < nRun; rn++) {
                    //System.out.println("Run: " + Integer.toString(rn + 1));

                    // generate random rates
                    insRate = new BigDecimal(Double.toString(Math.random()));
                    delRate = new BigDecimal(Double.toString(Math.random()));

                    // call the MCMC function
                    MCMCOutput mcmcOutput = networkMCMC.networkMCMC(startNetwork, endNetwork, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), evolTimeMCMC, nIter, nBurning, enumStates, updateInterval, evolutionModel, true, dependenceProbability, false);
                    ArrayList pathList = mcmcOutput.getPaths();
                    ArrayList networkList = mcmcOutput.getVisitedNetworks();
                    BigDecimal[] pathProbList = mcmcOutput.getPathProbabilities();
                    Double[] insRates = mcmcOutput.getInsRates();
                    Double[] delRates = mcmcOutput.getDelRates();
                    int acceptCountPath = mcmcOutput.getAcceptCountPath();
                    int acceptCountParameter = mcmcOutput.getAcceptCountParameter();
                    //System.out.println(new Date());
                    //NetworkEvolver.summariseEvolution(startNetwork, networkList);
                    //System.out.println(new Date());
                    //HashMap pathDistribution = NetworkEvolver.pathDistribution(pathList, nBurning, true);
                    //HashMap networkDistribution = NetworkEvolver.networkDistribution(startNetwork, new ArrayList(pathDistribution.keySet()), true);
                    //HashMap insRateDistribution = NetworkEvolver.rateDistribution(insRates, nBurning, false, "Insertion Rate:");
                    //HashMap delRateDistribution = NetworkEvolver.rateDistribution(delRates, nBurning, false, "Deletion Rate:");
                    //NetworkEvolver.pathLengthDistribution(pathList, nBurning, true);
                    HashMap likelihood = networkMCMC.calculateLikelihood(pathList, pathProbList, insRates, delRates, nBurning, false);

                    // find max likelihood 
                    Iterator it = likelihood.keySet().iterator();
                    BigDecimal maxLH = new BigDecimal("0.0");
                    Double[] rates = new Double[2];
                    while (it.hasNext()) {
                        String strRate = (String) it.next();
                        BigDecimal lh = (BigDecimal) likelihood.get(strRate);
                        if (lh.compareTo(maxLH) > 0) {
                            String[] strRates = strRate.split("_");
                            rates[0] = new Double(strRates[0]);
                            rates[1] = new Double(strRates[1]);
                        }
                    }
                    Double ratio;
                    if (rates[0] > rates[1])
                        ratio = rates[0] / rates[1];
                    else
                        ratio = rates[1] / rates[0];
                    System.out.println(Integer.toString(i + 1) + " to " + Integer.toString(i + stepSize) + "\t" + rates[0].toString() + "\t" + rates[1].toString() + "\t" + Double.toString(evolTime) + "\t" + Double.toString(Utilities.round(ratio, 6)) + "\t" + Double.toString(rates[0] * evolTime) + "\t" + Double.toString(rates[1] * evolTime));
                    // save the ratio
                    //arrRatio[(i+1)/stepSize] = ratio;
                    outResult.write(Integer.toString(i + 1) + " to " + Integer.toString(i + stepSize) + "\t" + rates[0].toString() + "\t" + rates[1].toString() + "\t" + Double.toString(evolTime) + "\t" + Double.toString(Utilities.round(ratio, 6)) + "\t" + Double.toString(rates[0] * evolTime) + "\t" + Double.toString(rates[1] * evolTime));
                    outResult.write("\r\n");
                    outResult.flush();

                }// end for rn=1 to nRun
            } // for i = stepSize to nIterSim 

            outResult.close(); // Also flushes output
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    static private void testPathSamplerForPathLengths(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename, double dependenceProbability) {

        double[] insRates = new double[]{0.05, 0.5, 5.0};
        double[] delRates = new double[]{0.03, 0.3, 3.0};

        for (int i = 0; i < insRates.length; i++)
            testPathSamplerForPathLengths(evolutionModel, dirName, refNWFilename, startNWFilename, endNWFilename, insRates[i], delRates[i], dependenceProbability);
    }

    static private void testPathSamplerForPathLengths(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename,
            double insRate, double delRate, double dependenceProbability) {
        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();

        // calculate the likelihood conditioned on start network (H1->H2)
        runPathSampler(evolutionModel, startNetwork, endNetwork, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
    }

    static private void testPathSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename, double dependenceProbability) {

        double[] insRates = new double[]{0.1, 0.3, 0.5, 0.7, 0.9};
        // mu = 2*lambda
        for (int i = 0; i < insRates.length; i++)
            testPathSampler(evolutionModel, dirName, refNWFilename, startNWFilename, endNWFilename, insRates[i], insRates[i] * 2.0, dependenceProbability);
        // mu = lambda
        for (int i = 0; i < insRates.length; i++)
            testPathSampler(evolutionModel, dirName, refNWFilename, startNWFilename, endNWFilename, insRates[i], insRates[i], dependenceProbability);
        // mu = lambda/2
        for (int i = 0; i < insRates.length; i++)
            testPathSampler(evolutionModel, dirName, refNWFilename, startNWFilename, endNWFilename, insRates[i], insRates[i] / 2.0, dependenceProbability);
    }

    static private void testPathSampler(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename,
            double insRate, double delRate, double dependenceProbability) {
        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();


        // calculate the likelihood conditioned on start network (H1->H2 & H2->H1)
        runPathSampler(evolutionModel, startNetwork, endNetwork, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
        runPathSampler(evolutionModel, endNetwork, startNetwork, coreEdges, prohibEdges, insRate, delRate, dependenceProbability);
        // generate rate matrix for stationary distribtion calculation
        generateRateMatrix_ToyNetworks(evolutionModel, refNetwork, insRate, delRate, dependenceProbability);
    }

    /* ********************************************************************** *
     *            R A T E    M A T R I X    G E N E R A T I O N               *
     * ********************************************************************** */
    static private void generateRateMatrix_ToyNetworks(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename,
            double insRate, double delRate, double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        String strStartNetwork = Utilities.toString(startNetwork.getReactionSequence(), "");
        String strEndNetwork = Utilities.toString(endNetwork.getReactionSequence(), "");
        int startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
        int endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);

        generateRateMatrix_ToyNetworks(evolutionModel, refNetwork, insRate, delRate, dependenceProbability);

        System.out.println("Start network index: " + Integer.toString(startIdx));
        System.out.println("End network index: " + Integer.toString(endIdx));

    }

    static private void generateRateMatrix_ToyNetworks(int evolutionModel,
            MetabolicNetwork refNetwork, double insRate, double delRate, double dependenceProbability) {

        int edgeCount = refNetwork.getDirectedReactions().size();
        ArrayList lstNetworks = new NetworkEvolver().enumerateAllNetworks(edgeCount);

        double[][] rateMatrix = new NetworkEvolver().generateRateMatrix(lstNetworks, refNetwork, evolutionModel, insRate, delRate, dependenceProbability);

        try {
            String filename = "RateMatrix_toyNetwork_" + Double.toString(insRate) + "_" + Double.toString(delRate) + ".txt";
            // create the output file
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(filename)));
            for (int i = 0; i < rateMatrix.length; i++) {
                double[] ds = rateMatrix[i];
                for (int j = 0; j < ds.length - 1; j++) {
                    out.write(new BigDecimal(ds[j]).round(new MathContext(6)).toString() + "\t");
                }
                // write the last element
                out.write(new BigDecimal(ds[ds.length - 1]).round(new MathContext(6)).toString());
                out.write("\r\n");

            }
            out.close(); // Also flushes output
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    static private void generateRateMatrices_CornerCutting_ToyNetworks(int evolutionModel,
            String dirName, String refNWFilename, String startNWFilename, String endNWFilename,
            double insRate, double delRate, double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();


        for (int extraEvents = 0; extraEvents < 5; extraEvents++) {
            ArrayList lstNetworks = new NetworkEvolver().enumerateNetworks(startNetwork, endNetwork, new ArrayList(),
                    new ArrayList(), extraEvents);
            //startIdx = 0;
            //endIdx = 1;

            System.out.println("Networks enumerated. " + new Date().toString());

            double[][] rateMatrix = null;
            System.out.println("Insertion Rate: " + Double.toString(insRate) + ", Deletion Rate: " + Double.toString(delRate));
            rateMatrix = new NetworkEvolver().generateRateMatrix(lstNetworks, refNetwork, evolutionModel, insRate, delRate, dependenceProbability);

            try {
                String filename = "toy_rateMatrix_" + Double.toString(insRate) + "_" + Double.toString(delRate) + "_" + Integer.toString(lstNetworks.size()) + ".txt";
                // create the output file
                BufferedWriter out = new BufferedWriter(new FileWriter(new File(filename)));
                for (int i = 0; i < rateMatrix.length; i++) {
                    double[] ds = rateMatrix[i];
                    for (int j = 0; j < ds.length - 1; j++) {
                        out.write(new BigDecimal(ds[j]).round(new MathContext(6)).toString() + "\t");
                    }
                    // write the last element
                    out.write(new BigDecimal(ds[ds.length - 1]).round(new MathContext(6)).toString());
                    out.write("\r\n");

                }
                out.close(); // Also flushes output
            } catch (Exception e) {
                e.printStackTrace();
            }

        }
    }

    static private void generateRateMatrices_LH_Surface_ToyNetworks(int evolutionModel,
            String dirName, String refNWFilename, String startNWFilename, String endNWFilename, 
            double dependenceProbability) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        String strStartNetwork = Utilities.toString(startNetwork.getReactionSequence(), "");
        String strEndNetwork = Utilities.toString(endNetwork.getReactionSequence(), "");
        int startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
        int endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);

        System.out.println("Start Network Idx: " + Integer.toString(startIdx));
        System.out.println("End Network Idx: " + Integer.toString(endIdx));


        int edgeCount = refNetwork.getDirectedReactions().size();
        ArrayList lstNetworks = new NetworkEvolver().enumerateAllNetworks(edgeCount);

        System.out.println("Networks enumerated. " + new Date().toString());


        final double dummy_ir = 1.0;
        final double dummy_dr = 2.0;
        double[][] rmIE = null;
        String filename = "IE_toy_rateMatrix_" + Double.toString(dummy_ir) + "_" + Double.toString(dummy_dr) + ".jdat";
        try {
            if (new File(filename).exists()) {
                System.out.println("Indepdent Edge rate matrix file exists ... reading it");
                ObjectInputStream in =
                        new ObjectInputStream(
                        new FileInputStream(new File(filename)));
                rmIE = (double[][]) in.readObject();
                in.close();

            } else {
                rmIE = new NetworkEvolver().generateRateMatrix(lstNetworks, refNetwork, NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE, dummy_ir, dummy_dr, 0.0);

                ObjectOutputStream out =
                        new ObjectOutputStream(
                        new FileOutputStream(new File(filename)));
                out.writeObject(rmIE);
                out.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        double[][] rmModel = null;
        filename = "ND_toy_rateMatrix_" + Double.toString(1.0) + "_" + Double.toString(1.0) + ".jdat";
        try {
            if (new File(filename).exists()) {
                System.out.println("Neighbour weight file exists ... reading it");
                ObjectInputStream in =
                        new ObjectInputStream(
                        new FileInputStream(new File(filename)));
                rmModel = (double[][]) in.readObject();
                in.close();
            } else {
                rmModel = new NetworkEvolver().generateRateMatrix(lstNetworks, refNetwork, evolutionModel, 1.0, 1.0, dependenceProbability);

                ObjectOutputStream out =
                        new ObjectOutputStream(
                        new FileOutputStream(new File(filename.replace(".txt", ".jdat"))));
                out.writeObject(rmModel);
                out.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        BigDecimal insRate = new BigDecimal("0.4");
        BigDecimal stepSizeIR = new BigDecimal("0.1");
        BigDecimal stepSizeDR = new BigDecimal("0.04");
        int nDim = 20;
        for (int ir = 0; ir < nDim; ir++) {
            BigDecimal delRate = new BigDecimal("0.10");
            for (int dr = 0; dr < nDim; dr++) {
                System.out.println("Insertion Rate: " + insRate.toString() + ", Deletion Rate: " + delRate.toString());
                double[][] rateMatrix = new double[rmIE.length][rmIE[0].length];
                for (int i = 0; i < rateMatrix.length; i++) {
                    for (int j = 0; j < rateMatrix[0].length; j++) {
                        double rate = 0.0;
                        if (i == j)
                            rate = 0.0;
                        else if (rmIE[i][j] == dummy_ir)
                            rate = insRate.doubleValue();
                        else if (rmIE[i][j] == dummy_dr)
                            rate = delRate.doubleValue();
                        rateMatrix[i][j] = (1- dependenceProbability) * rate + dependenceProbability * rate * rmModel[i][j];
                    }
                    // set the diagonal element
                    rateMatrix[i][i] = -Utilities.sum(rateMatrix[i]);
                }
                try {
                    filename = "toy_rateMatrix_" + insRate.toString() + "_" + delRate.toString() + "_" + Double.toString(dependenceProbability) + ".txt";
                    // create the output file
                    String resultDirName = new File(".").getCanonicalPath();
                    BufferedWriter out = new BufferedWriter(new FileWriter(new File(resultDirName, filename)));
                    for (int i = 0; i < rateMatrix.length; i++) {
                        double[] ds = rateMatrix[i];
                        for (int j = 0; j < ds.length - 1; j++) {
                            out.write(new BigDecimal(ds[j]).round(new MathContext(6)).toString() + "\t");
                        }
                        // write the last element
                        out.write(new BigDecimal(ds[ds.length - 1]).round(new MathContext(6)).toString());
                        out.write("\r\n");
                    }
                    out.close(); // Also flushes output
                } catch (Exception e) {
                    e.printStackTrace();
                }
                delRate = delRate.add(stepSizeDR);
            } // end for (deletion rates)
            insRate = insRate.add(stepSizeIR);
        } // end for (insertion rates)

    }

    /* ********************************************************************** *
     *            E Q U I L I B R I U M     P R O B A B I L I T Y             *
     * ********************************************************************** */
    static private void testEquilibriumProbability(int evolutionModel, String dirName,
            String refNWFilename, String nwFilename, double dependenceProbability) {

        double[] insRates = new double[]{0.1, 0.3, 0.5, 0.7, 0.9};
        // mu = 2*lambda
        for (int i = 0; i < insRates.length; i++)
            equilibriumProbability(evolutionModel, dirName, refNWFilename, nwFilename, insRates[i], insRates[i] * 2.0, dependenceProbability);
        // mu = lambda
        for (int i = 0; i < insRates.length; i++)
            equilibriumProbability(evolutionModel, dirName, refNWFilename, nwFilename, insRates[i], insRates[i], dependenceProbability);
        // mu = lambda/2
        for (int i = 0; i < insRates.length; i++)
            equilibriumProbability(evolutionModel, dirName, refNWFilename, nwFilename, insRates[i], insRates[i] / 2.0, dependenceProbability);
    }

    static private void testEquilibriumProbability_SubNetworkSize(int evolutionModel,
            String dirName, String refNWFilename, String nwFilename, double dependenceProbability) {

        double insRate = 0.1;
        int[] subNetworkSizes = new int[]{2, 4, 6, 8};
        // mu = 2*lambda
        for (int i = 0; i < subNetworkSizes.length; i++)
            equilibriumProbability(evolutionModel, dirName, refNWFilename, nwFilename, insRate, insRate * 2.0, subNetworkSizes[i], dependenceProbability);
        // mu = lambda
        for (int i = 0; i < subNetworkSizes.length; i++)
            equilibriumProbability(evolutionModel, dirName, refNWFilename, nwFilename, insRate, insRate, subNetworkSizes[i], dependenceProbability);
        // mu = lambda/2
        for (int i = 0; i < subNetworkSizes.length; i++)
            equilibriumProbability(evolutionModel, dirName, refNWFilename, nwFilename, insRate, insRate / 2.0, subNetworkSizes[i], dependenceProbability);
    }

    static private void equilibriumProbability(int evolutionModel, String dirName,
            String refNWFilename, String nwFilename, double insRate, double delRate, double dependenceProbability) {

        equilibriumProbability(evolutionModel, dirName, refNWFilename, nwFilename, insRate, delRate, -1, dependenceProbability);
    }

    static private void equilibriumProbability(int evolutionModel, String dirName,
            String refNWFilename, String nwFilename, double insRate, double delRate,
            int subNetworkSize, double dependenceProbability) {

        // load the metabolic networks
        // reference network
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        // given network
        MetabolicNetwork theNetwork = loadNetwork(dirName, nwFilename);
        theNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            theNetwork.setupNetworkMatrices();
        theNetwork.setupNetworkSequence();

        if (subNetworkSize < 0)
            subNetworkSize = refNetwork.getReactionSequence().length;

//        String strTheNetwork = Utilities.toString(theNetwork.getReactionSequence(), "");
//        int nwIdx = Integer.parseInt(new StringBuffer(strTheNetwork).reverse().toString(), 2);

        // Equillibrium probability
        double eqProb = new NetworkEvolver().approximateEquilibriumProbability(theNetwork, refNetwork,
                evolutionModel, insRate, delRate, subNetworkSize, dependenceProbability);

//        System.out.print(Integer.toString(nwIdx) + "\t" + Utilities.toString(theNetwork.getReactionSequence()));
        System.out.print("\tlambda = " + Double.toString(insRate) + ", mu = " + Double.toString(delRate));
        System.out.println("\t" + Double.toString(eqProb));

    }


    /* ********************************************************************** *
     *             T R A N S I T I O N     P R O B A B I L I T Y              *
     * ********************************************************************** */
    static private void transitionProbability(int evolutionModel, String dirName, String refNWFilename,
            String startNWFilename, String endNWFilename, String coreFilename, String prohibFilename,
            double insRate, double delRate, double evolTime, int subNetworkSize, double dependenceProbability) {

        // load the metabolic networks
        // reference network
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        int startIdx, endIdx;
        try {
            // start and end network indices
            String strStartNetwork = Utilities.toString(startNetwork.getReactionSequence(), "");
            String strEndNetwork = Utilities.toString(endNetwork.getReactionSequence(), "");
            startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
            endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);
        } catch (Exception ex) {
            startIdx = 0;
            endIdx = 0;
        }

        // core edges indices
        Byte[] iCoreEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iCoreEdges.length; i++)
            iCoreEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        if (coreFilename != null && coreFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, coreFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                iCoreEdges[idx] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            }
        }

        // prohibited edges indices
        Byte[] iProhibEdges = new Byte[refNetwork.getDirectedReactions().size()];
        for (int i = 0; i < iProhibEdges.length; i++)
            iProhibEdges[i] = MetabolicNetwork.SEQ_ENTRY_ABSENT;
        if (prohibFilename != null && prohibFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, prohibFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                iProhibEdges[idx] = MetabolicNetwork.SEQ_ENTRY_PRESENT;
            }
        }

        if (subNetworkSize < 0)
            subNetworkSize = refNetwork.getReactionSequence().length;


        // Transition probability
        Date start = new Date();
        double transProb = new NetworkEvolver().approximateTransitionProbability(startNetwork,
                endNetwork, refNetwork, iCoreEdges, iProhibEdges, evolutionModel, insRate, delRate,
                evolTime, subNetworkSize, dependenceProbability);
        Date end = new Date();

        System.out.print("Start: " + Integer.toString(startIdx) + "\t" + Utilities.toString(startNetwork.getReactionSequence()));
        System.out.print("\tEnd: " + Integer.toString(endIdx) + "\t" + Utilities.toString(endNetwork.getReactionSequence()));
        System.out.print("\tlambda = " + Double.toString(insRate) + ", mu = " + Double.toString(delRate));
        //System.out.print("\t" + Long.toString(end.getTime() - start.getTime()));
        System.out.println("\t" + Double.toString(transProb));

    }

//    static public void transitionProbability(int evolutionModel, MetabolicNetwork refNetwork,
//            MetabolicNetwork startNetwork, MetabolicNetwork endNetwork, Byte[] iCoreEdges, 
//            Byte[] iProhibEdges, double insRate, double delRate, double evolTime, int subNetworkSize) {
//        
//        transProb = new NetworkEvolver().approximateTransitionProbability(startNetwork,
//                endNetwork, refNetwork, iCoreEdges, iProhibEdges, evolutionModel, insRate, delRate,
//                evolTime, subNetworkSize);
//    }
    /* ********************************************************************** *
     *                          P H Y L O G E N Y                             *
     * ********************************************************************** */
    static public void runGibbsSamplerOnPhylogeny(int evolutionModel, String dirName,
            String refNWFilename, String[] nwFilenames, String phylogeny,
            String coreFilename, String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int nBurning, int updateInterval, double dependenceProbability) {

        // call the main function with false for save data
        runGibbsSamplerOnPhylogeny(evolutionModel, dirName, refNWFilename, nwFilenames,
                phylogeny, coreFilename, prohibFilename, insRate, delRate, nRun, nIter, nBurning,
                updateInterval, false, dependenceProbability);
    }

    static private void runGibbsSamplerOnPhylogeny(int evolutionModel, String dirName,
            String refNWFilename, String[] nwFilenames, double dependenceProbability) {

        String phylogeny = "((H1:1.0,H2:1.0):1.0,H3:1.0)";

        // MCMC parameters
        int nIter = 11000, nBurning = 1000;
        int updateInterval = 10;
        int nRun = 3;
        // random starting values
        double insRate = -1.0;
        double delRate = -1.0;

        saveData = false;
        runGibbsSamplerOnPhylogeny(evolutionModel, dirName, refNWFilename, nwFilenames, phylogeny,
                null, null, insRate, delRate, nRun, nIter, nBurning, updateInterval, saveData, 
                dependenceProbability);
    }

    static private void runGibbsSamplerOnPhylogeny(int evolutionModel, String dirName,
            String refNWFilename, String[] nwFilenames, String coreFilename, String prohibFilename,
            double dependenceProbability) {

        String phylogeny = "((pfo:1.0,pst:1.0):1.0,pae:1.0)";

        // MCMC parameters
        int nIter = 5000, nBurning = 1000;
        int updateInterval = 10;
        int nRun = 1;
        // random starting values
        double insRate = -1.0;
        double delRate = -1.0;

        saveData = false;
        runGibbsSamplerOnPhylogeny(evolutionModel, dirName, refNWFilename, nwFilenames, phylogeny,
                coreFilename, prohibFilename, insRate, delRate, nRun, nIter, nBurning, updateInterval, saveData,
                dependenceProbability);
    }

    // main function
    static private void runGibbsSamplerOnPhylogeny(int evolutionModel, String dirName,
            String refNWFilename, String[] nwFilenames, String phylogeny,
            String coreFilename, String prohibFilename, double insRate, double delRate,
            int nRun, int nIter, int nBurning, int updateInterval, boolean save_datafile, 
            double dependenceProbability) {

        ArrayList lstNetworks = new ArrayList(nwFilenames.length);
        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        for (int i = 0; i < nwFilenames.length; i++) {
            // network 1
            MetabolicNetwork theNetwork = loadNetwork(dirName, nwFilenames[i]);
            theNetwork.addInactiveReactions(refNetwork.getReactions());
            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                theNetwork.setupNetworkMatrices();
            theNetwork.setupNetworkSequence();

            /*            String strNetwork = Utilities.toString(theNetwork.getReactionSequence(), "");
            int nwIdx = Integer.parseInt(new StringBuffer(strNetwork).reverse().toString(), 2);
            System.out.println("Network Idx - " + theNetwork.getID() + ": " + Integer.toString(nwIdx));
             */
            // add to the list of networks
            lstNetworks.add(theNetwork);
        }

        // core edges
        ArrayList coreEdges = new ArrayList();
        if (coreFilename != null && coreFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, coreFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                coreEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        // prohibited edges
        ArrayList prohibEdges = new ArrayList();
        if (prohibFilename != null && prohibFilename.length() != 0) {
            String[] edges = Utilities.readList(dirName, prohibFilename, false);
            for (int i = 0; i < edges.length; i++) {
                int rxnIdx = refNetwork.getReactions().indexOf(new Reaction(edges[i]));
                int idx = refNetwork.getDirectedReactions().indexOf(new DirectedReaction((Reaction) refNetwork.getReactions().get(rxnIdx), 'F'));
                prohibEdges.add(refNetwork.getDirectedReactions().get(idx));
            }
        }

        ArrayList lstOrganism = Organism.buildOrganismList(PhyloTree.extractNodesFromPhylogeny(phylogeny));
        PhyloTree phylotree = PhyloTree.buildPhylogeneticTree(phylogeny, lstOrganism);
        // Build Metabolic network on leaves
        assignNetworksOnTheLeaves(phylotree.getRoot(), lstNetworks);

        NetworkMCMC networkMCMC = new NetworkMCMC();

//        System.out.println(new Date());
//        System.out.println("Insertion Rate: " + insRate.toString() + ", Deletion Rate: " + delRate.toString());

        boolean sampleRates = true;
        boolean sampleDP = false;
        if (evolutionModel == NetworkEvolver.EVOL_MODEL_HYBRID)
            sampleDP = true;
        else if (evolutionModel == NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT)
            dependenceProbability = 1.0;
        
        int decimalPlaces = 5;
        int subNetworkSize = 4;
        for (int rn = 0; rn < nRun; rn++) {
            System.out.println(new Date());
            System.out.println("Run: " + Integer.toString(rn + 1));
            if (sampleRates) {
                // generate random rates
                if (insRate < 0)
                    insRate = Utilities.round(Math.random(), decimalPlaces);
                if (delRate < 0)
                    delRate = Utilities.round(Math.random(), decimalPlaces);
            }
            if (sampleDP) {
                // generate random probability value
                if (dependenceProbability < 0)
                    dependenceProbability = Utilities.round(Math.random(), decimalPlaces);
            }

            networkMCMC.setDistinctTrees(new HashMap((int) (nIter * 0.5)));
            // call the MCMC function
            Date startTime = new Date();
            networkMCMC.networkMCMC(phylotree, refNetwork, coreEdges, prohibEdges, insRate, delRate, nIter, updateInterval, evolutionModel, nBurning, sampleRates, subNetworkSize, dependenceProbability, sampleDP);
            Date endTime = new Date();
            // extract insertion and deletion rates from the root
            MCMCOutput mcmcOutput = (MCMCOutput) phylotree.getRoot().getData();
            Double[] insRates = mcmcOutput.getInsRates();
            Double[] delRates = mcmcOutput.getDelRates();
            int acceptCountParameter = mcmcOutput.getAcceptCountParameter();
            int acceptCountDP = mcmcOutput.getAcceptCountDependenceProbability();
            NetworkEvolver.rateDistribution(insRates, nBurning, true, "Insertion Rate:");
            NetworkEvolver.rateDistribution(delRates, nBurning, true, "Deletion Rate:");
            //networkMCMC.calculateLikelihood(phylotree, refNetwork, coreEdges, prohibEdges, evolutionModel, nBurning, true, dependenceProbability);
            //networkMCMC.calculateAutocorrelationPhylo(mcmcOutput, true);
            System.out.println("Acceptance Count - Parameter: " + Integer.toString(acceptCountParameter));
            System.out.println("Acceptance Count - Dependence Probabilty: " + Integer.toString(acceptCountDP));
            System.out.println("Total Time: " + Long.toString(endTime.getTime() - startTime.getTime()) + " milliseconds");
            System.out.println();

            if (save_datafile) {
            // save data?
            }
        }// end for rn=1 to nRun
//        if (saveData) {
//            try {
//                out.close(); // Also flushes output
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        }

    }

    static private void runGibbsSamplerOnPhylogeny_LH_Surface(int evolutionModel, String dirName,
            String refNWFilename, String[] nwFilenames, double dependenceProbability) {

        String phylogeny = "((H1:1.0,H2:1.0):1.0,H3:1.0)";
        // load the metabolic networks
        ArrayList lstNetworks = new ArrayList(nwFilenames.length);
        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        for (int i = 0; i < nwFilenames.length; i++) {
            // network 1
            MetabolicNetwork theNetwork = loadNetwork(dirName, nwFilenames[i]);
            theNetwork.addInactiveReactions(refNetwork.getReactions());
            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                theNetwork.setupNetworkMatrices();
            theNetwork.setupNetworkSequence();

            String strNetwork = Utilities.toString(theNetwork.getReactionSequence(), "");
            int nwIdx = Integer.parseInt(new StringBuffer(strNetwork).reverse().toString(), 2);

            System.out.println("Network Idx - " + theNetwork.getID() + ": " + Integer.toString(nwIdx));

            // add to the list of networks
            lstNetworks.add(theNetwork);
        }

        ArrayList lstOrganism = Organism.buildOrganismList(PhyloTree.extractNodesFromPhylogeny(phylogeny));
        PhyloTree phylotree = PhyloTree.buildPhylogeneticTree(phylogeny, lstOrganism);
        // Build Metabolic network on leaves
        assignNetworksOnTheLeaves(phylotree.getRoot(), lstNetworks);

        // MCMC parameters
        int nIter = 60000, nBurning = 10000;
        int updateInterval = 1000;
        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();
        boolean sampleRates = false;
        int subNetworkSize = 6;
        NetworkMCMC networkMCMC = new NetworkMCMC();

        String saveDirName;
        try {
            saveDirName = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            saveDirName = "";
        }
        BigDecimal insRate = new BigDecimal("0.40");
        BigDecimal stepSizeIR = new BigDecimal("0.10");
        BigDecimal stepSizeDR = new BigDecimal("0.04");
        int nDim = 20;
        BigDecimal[][] LH_MCMC = new BigDecimal[nDim][nDim];

        try {
            String outFilename = "LH_MCMC_Phylo.txt";
            // create the output file
            BufferedWriter out = new BufferedWriter(new FileWriter(new File(saveDirName, outFilename)));

            for (int ir = 0; ir < nDim; ir++) {
                BigDecimal delRate = new BigDecimal("0.10");
                for (int dr = 0; dr < nDim; dr++) {
                    System.out.println(new Date());
                    System.out.println("Insertion Rate: " + insRate.toString() + ", Deletion Rate: " + delRate.toString());

                    // call the MCMC function
                    networkMCMC.setDistinctTrees(new HashMap((int) (nIter * 0.5)));
                    networkMCMC.networkMCMC(phylotree, refNetwork, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), nIter, updateInterval, evolutionModel, nBurning, sampleRates, subNetworkSize, dependenceProbability, false);
                    double likelihood = networkMCMC.calculateLikelihood(phylotree, refNetwork, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), evolutionModel, nBurning, dependenceProbability);

                    //LH_MCMC[ir][dr] = (BigDecimal) likelihood.get(new ArrayList(likelihood.keySet()).get(0));
                    LH_MCMC[ir][dr] = new BigDecimal(likelihood);

                    System.out.println("Likelihood: " + LH_MCMC[ir][dr].round(new MathContext(16)).toString());

                    if (saveData) {
                    }

                    delRate = delRate.add(stepSizeDR);
                }

                // write the likelihood values for current insertion rate
                BigDecimal[] ds = LH_MCMC[ir];
                for (int j = 0; j < ds.length - 1; j++) {
                    out.write(ds[j].round(new MathContext(16)).toString() + "\t");
                }
                // write the last element
                out.write(ds[ds.length - 1].round(new MathContext(16)).toString());
                out.write("\r\n");
                out.flush();


                insRate = insRate.add(stepSizeIR);

            } // end for each insertion rate 

            out.close(); // Also flushes output
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println(Utilities.toString(LH_MCMC, "\t"));

    }

    static private void assignNetworksOnTheLeaves(PhyloNode phylonode, ArrayList lstNetworks) {

        if (!phylonode.isLeaf()) {
            assignNetworksOnTheLeaves(phylonode.getLeftSon(), lstNetworks);
            assignNetworksOnTheLeaves(phylonode.getRightSon(), lstNetworks);
        } else {
            Organism org = (Organism) phylonode.getNetworkObject();
            int idx = lstNetworks.indexOf(org);
            if (idx >= 0)
                phylonode.setMetabolicNetwork((MetabolicNetwork) lstNetworks.get(idx));
        }
    }

//    static private void testGibbsSamplerOnPhylogeny(int evolutionModel, String dirName,
//            String refNWFilename, String[] nwFilenames) {
//
//        ArrayList lstNetworks = new ArrayList(nwFilenames.length);
//        // load the metabolic networks
//        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
//        refNetwork.setupNetworkMatrices();
//        refNetwork.setupNetworkSequence();
//        for (int i = 0; i < nwFilenames.length; i++) {
//            // network 1
//            MetabolicNetwork theNetwork = loadNetwork(dirName, nwFilenames[i]);
//            theNetwork.addInactiveReactions(refNetwork.getReactions());
//            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
//                theNetwork.setupNetworkMatrices();
//            theNetwork.setupNetworkSequence();
//
//            String strNetwork = Utilities.toString(theNetwork.getReactionSequence(), "");
//            int nwIdx = Integer.parseInt(new StringBuffer(strNetwork).reverse().toString(), 2);
//
//            System.out.println("Network Idx - " + theNetwork.getID() + ": " + Integer.toString(nwIdx));
//
//            // add to the list of networks
//            lstNetworks.add(theNetwork);
//        }
//
//        // core edges
//        ArrayList coreEdges = new ArrayList();
//        // prohibited edges
//        ArrayList prohibEdges = new ArrayList();
//
//        String phylogeny = "((H1:1.0,H2:1.0):1.0,H3:1.0)";
//        ArrayList lstOrganism = Organism.buildOrganismList(PhyloTree.extractNodesFromPhylogeny(phylogeny));
//        PhyloTree phylotree = PhyloTree.buildPhylogeneticTree(phylogeny, lstOrganism);
//        // Build Metabolic network on leaves
//        assignNetworksOnTheLeaves(phylotree.getRoot(), lstNetworks);
//
//        NetworkMCMC networkMCMC = new NetworkMCMC();
//
//        // Rates
//        double insRate = 0.06;
//        double delRate = 0.03;
//
//        // Simulation parameters
//        int nBurningSim = 2500;
//        int nIterSim = 5000;
//        int nUpdateSim = 10;
//        int subNetworkSize = 6;
//        boolean sampleRates = false;
//        // simulate the network evolution
//        networkMCMC.networkMCMC(phylotree, refNetwork, coreEdges, prohibEdges, insRate, delRate, nIterSim, nUpdateSim, evolutionModel, nBurningSim, sampleRates, subNetworkSize);
//
//        // Initialise the intermediary variables
//        networkMCMC.setDistinctNetworkPairs(new HashMap(nIterSim));
//        networkMCMC.setDistinctTrees(new HashMap(nIterSim));
//        
//        System.out.println();
//        System.out.println("Re-estimating parameters ...");
//        // MCMC Parameters
//        int nIter = 600;
//        int nBurning = 100;
//        int stepSize = 10;
//        for (int iter = nIterSim - nBurningSim + 1; iter < nIterSim; iter += stepSize) {
//            // set the networks at internal nodes
//            networkMCMC.setNetworksAtInternalNodes(phylotree.getRoot(), refNetwork, iter);
//            // estimate the rates
//            Double[] mlEstimate = networkMCMC.estimateRates(phylotree, refNetwork, coreEdges, prohibEdges, nIter, evolutionModel, nBurning, subNetworkSize);
//
//            insRate = Utilities.round(mlEstimate[0], 6);
//            delRate = Utilities.round(mlEstimate[1], 6);
//            System.out.print("Iteration: " + Integer.toString(iter));
//            System.out.println("\t" + Double.toString(insRate) + "\t" + Double.toString(delRate));
//        }
//
//    }
    static private void testGibbsSamplerOnPhylogeny(int evolutionModel, String dirName,
            String refNWFilename, String[] nwFilenames, double dependenceProbability) {

        ArrayList lstNetworks = new ArrayList(nwFilenames.length);
        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        for (int i = 0; i < nwFilenames.length; i++) {
            // network 1
            MetabolicNetwork theNetwork = loadNetwork(dirName, nwFilenames[i]);
            theNetwork.addInactiveReactions(refNetwork.getReactions());
            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                theNetwork.setupNetworkMatrices();
            theNetwork.setupNetworkSequence();

            String strNetwork = Utilities.toString(theNetwork.getReactionSequence(), "");
            int nwIdx = Integer.parseInt(new StringBuffer(strNetwork).reverse().toString(), 2);

            System.out.println("Network Idx - " + theNetwork.getID() + ": " + Integer.toString(nwIdx));

            // add to the list of networks
            lstNetworks.add(theNetwork);
        }

        // core edges
        ArrayList coreEdges = new ArrayList();
        // prohibited edges
        ArrayList prohibEdges = new ArrayList();

        String phylogeny = "((H1:1.0,H2:1.0):1.0,H3:1.0)";
        ArrayList lstOrganism = Organism.buildOrganismList(PhyloTree.extractNodesFromPhylogeny(phylogeny));
        PhyloTree phylotree = PhyloTree.buildPhylogeneticTree(phylogeny, lstOrganism);
        // Build Metabolic network on leaves
        assignNetworksOnTheLeaves(phylotree.getRoot(), lstNetworks);

        NetworkMCMC networkMCMC = new NetworkMCMC();

        // Rates
        double insRate = 0.03;
        double delRate = 0.03;

        // Simulation parameters
        int nBurningSim = 25000;
        int nIterSim = 75000;
        int nUpdateSim = 1000;
        int subNetworkSize = 6;
        boolean sampleRates = false;

        networkMCMC.setDistinctTrees(new HashMap(1000));
        // simulate the network evolution
        networkMCMC.networkMCMC(phylotree, refNetwork, coreEdges, prohibEdges, insRate, delRate, nIterSim, nUpdateSim, evolutionModel, nBurningSim, sampleRates, subNetworkSize, dependenceProbability, false);

        // Initialise the intermediary variables
        networkMCMC.setDistinctNetworkPairs(new HashMap(nIterSim));
        networkMCMC.setDistinctTrees(new HashMap(nIterSim));

        int nUpdate = 10;
        System.out.println();
        System.out.println("Re-estimating parameters ...");
        networkMCMC.estimateRates(phylotree, refNetwork, coreEdges, prohibEdges, nIterSim, nUpdate, evolutionModel, nBurningSim, subNetworkSize, dependenceProbability);

    }

    /* ********************************************************************** *
     *     EXPECTED EVENTS VS MINIMUM DISTANCE - NEIGBOUR DEPENDENT MODEL     *
     * ********************************************************************** */
    static private void calculatedExpectedNumberOfEvents(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, double insRate, double delRate, double dependenceProbability) {

        System.out.println("Insertion Rate: " + Double.toString(insRate) + ", Deletion Rate: " + Double.toString(delRate));

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();

        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        int nIter = 60000, nBurning = 10000;
        int updateInterval = 0;
        NetworkSimulator simulator = new NetworkSimulator();
        ArrayList lstNetworks = simulator.simulateNetworkEvolution(startNetwork, evolutionModel, nIter, nBurning,
                new ArrayList(), new ArrayList(), insRate, delRate, updateInterval, dependenceProbability);
        Object[] arrNetworks = lstNetworks.toArray();

        System.out.println("Compilation Started");
        int MAX_DIFF = (int) (refNetwork.getDirectedReactions().size() / 2);
        int MAX_EVENTS = MAX_DIFF * 10;
        int[][] diffVsEvents = new int[MAX_DIFF + 1][MAX_EVENTS + 1];
        diffVsEvents[1][1] = arrNetworks.length - 2; // subtract two since we ignore the first element and number of pairs is 'count-1'
        for (int i = 1; i < arrNetworks.length; i++) { // ignore the first element
            Byte[] network1 = (Byte[]) arrNetworks[i];

            for (int e = 2; e <= MAX_EVENTS; e++) { // no point in doing it for e=1 since in this case d=1
                if (i + e >= arrNetworks.length)
                    break;

                Byte[] network2 = (Byte[]) arrNetworks[i + e];

                int d = MetabolicNetwork.findDifferenceCount(network1, network2);
                if (d > MAX_DIFF)
                    continue;

                diffVsEvents[d][e]++;
            }
        }

        // calculate the weighted average
        double[][] expectedEvents = new double[diffVsEvents.length - 1][2];
        for (int d = 1; d < diffVsEvents.length; d++) {
            int[] ds = diffVsEvents[d];
            double sumprod = 0.0;
            double sum = 0.0;
            for (int j = 1; j < ds.length; j++) {
                sumprod += j * ds[j];
                sum += ds[j];
            }
            expectedEvents[d - 1][0] = d;
            expectedEvents[d - 1][1] = sumprod / sum;
        }

        System.out.println(Utilities.toString(expectedEvents));
    }

    /* ********************************************************************** *
     *                      E    M    A L G O R I T H M   ??                  *
     * ********************************************************************** */
    static private void emAlgorithm_ToyNetworks_SimulatedData(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename) {


        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // core and prohibited edges
        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();
        // Rates
        Double insRate = new Double("0.9");
        Double delRate = new Double("0.3");
        double dependenceProbability = 0.5;
        // Simulate the data
        int nIterSim = 1000;
        int nBurningSim = 50000;
        //ArrayList lstSimulation = new NetworkSimulator().simulateNetworkEvolutionWithTime(startNetwork, evolutionModel, nIterSim, nBurningSim, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue());
        ArrayList lstSimulation = new NetworkSimulator().simulateNetworkEvolutionWithTime(startNetwork, evolutionModel, nIterSim, nBurningSim, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), dependenceProbability);
        ArrayList lstSimulatedNetworks = (ArrayList) lstSimulation.get(0);
        ArrayList lstSimulatedTimes = (ArrayList) lstSimulation.get(1);

        try {
            String resultSaveDirName = new File(".").getCanonicalPath();
            DateFormat dateFormat = new SimpleDateFormat("yyyyMMdd_HHmmss");
            String resultfilename = "em_simulatedData_" + dateFormat.format(new Date()) + ".txt";
            // create the output file
            BufferedWriter outResult = new BufferedWriter(new FileWriter(new File(resultSaveDirName, resultfilename)));

            Iterator itTimes = lstSimulatedTimes.iterator();
            Iterator itNetworks = lstSimulatedNetworks.iterator();
            int idx = 0;
            while (itNetworks.hasNext()) {
                outResult.write(Integer.toString(idx++) + ":\t" + (Double) itTimes.next() + "\t" + Utilities.toString((Byte[]) itNetworks.next(), " ") + "\r\n");
            }
            outResult.flush();

            System.out.println("Estimating parameters");
            System.out.println("Iteration\t" + "sampleSize\t" + "lambda\t" + "mu");
            outResult.write("Iteration\t" + "sampleSize\t" + "lambda\t" + "mu");

            //System.out.println("Sample: " + Integer.toString(i+1) + " to " + Integer.toString(i + stepSize));

            MetabolicNetwork endNetwork = refNetwork.getNetwork((Byte[]) lstSimulatedNetworks.get(lstSimulatedNetworks.size() - 1));
            endNetwork.addInactiveReactions(refNetwork.getReactions());
            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                endNetwork.setupNetworkMatrices();
            endNetwork.setupNetworkSequence();

            double evolTime = (Double) lstSimulatedTimes.get(lstSimulatedTimes.size() - 1);
            double evolTimeMCMC = 1.0;

            // MCMC parameters
            int nBurning = 100;
            int updateInterval = 0;
            boolean enumStates = false;
            NetworkMCMC networkMCMC = new NetworkMCMC();

            // generate random rates
            Double newInsRate = new Double(Math.random());
            Double newDelRate = new Double(Math.random());

            // EM Parameters
            double tolerance = 0.0001;
            int blockSize = 5;
            int sampleSize = 0;
            int sampleSizeIncrement = 50;
            int iter = 0;
            insRate = 0.0;
            delRate = 0.0;
            while (Math.abs(insRate - newInsRate) > tolerance || Math.abs(delRate - newDelRate) > tolerance) {

                insRate = newInsRate;
                delRate = newDelRate;

                if (iter++ % blockSize == 0)
                    sampleSize += sampleSizeIncrement;
                /* ***************************************************** *
                 *                   E   -   S T E P                     *
                 * ***************************************************** */

                // call the MCMC function
                int nIter = sampleSize + nBurning;
                MCMCOutput mcmcOutput = networkMCMC.networkMCMC(startNetwork, endNetwork, coreEdges, prohibEdges, insRate.doubleValue(), delRate.doubleValue(), evolTimeMCMC, nIter, nBurning, enumStates, updateInterval, evolutionModel, false, dependenceProbability, false);
                ArrayList pathList = mcmcOutput.getPaths();
                //int acceptCount = mcmcOutput.getAcceptCount();

                int[] condMeans = calculateConditionalMeans(startNetwork, pathList, coreEdges, prohibEdges);

                /* ***************************************************** *
                 *                   E   -   S T E P                     *
                 * ***************************************************** */
                newInsRate = (double) condMeans[0] / ((double) condMeans[2]);
                newDelRate = (double) condMeans[1] / ((double) condMeans[3]);

                System.out.println(Integer.toString(iter) + "\t" + Integer.toString(sampleSize) + "\t" + newInsRate.toString() + "\t" + newDelRate.toString());
                outResult.write(Integer.toString(iter) + "\t" + Integer.toString(sampleSize) + "\t" + newInsRate.toString() + "\t" + newDelRate.toString());
                outResult.write("\r\n");
                outResult.flush();
            }// end while (convergence)
            outResult.close(); // Also flushes output
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    static private int[] calculateConditionalMeans(MetabolicNetwork startNetwork,
            ArrayList pathList, ArrayList coreEdges, ArrayList prohibEdges) {

        int nIns = 0;
        int nDel = 0;
        int nInsTotal = 0;
        int nDelTotal = 0;

        Iterator itPath = pathList.iterator();
        while (itPath.hasNext()) {
            int[] path = (int[]) itPath.next();

            // initialise
            int nInsTotalPath = 0;
            int nDelTotalPath = 0;

            MetabolicNetwork theNetwork = startNetwork.clone();
            for (int i = 0; i < path.length; i++) {
                // get the edge corresponding to ith event
                DirectedReaction reaction = (DirectedReaction) theNetwork.getDirectedReactions().get(path[i] - 1);

                // increment the appropriate insertion/deletion count
                if (theNetwork.getReactionSequence()[path[i] - 1] == MetabolicNetwork.SEQ_ENTRY_PRESENT)
                    nDel++;
                else
                    nIns++;

                int totalPresent = theNetwork.getActiveDirectedReactions().size();
                int totalAbsent = theNetwork.getDirectedReactions().size() - totalPresent;

                nInsTotalPath += totalAbsent - prohibEdges.size(); // insertable edges = no of edges absent - no. of prohib edges
                nDelTotalPath += totalPresent - coreEdges.size(); // deletable edges = no of edges present - no. of core edges

                // Update the network
                theNetwork = theNetwork.updateNetwork(reaction);
            } // end for each event in the path

            nInsTotal += (int) Math.round((double) nInsTotalPath / (double) path.length);
            nDelTotal += (int) Math.round((double) nDelTotalPath / (double) path.length);
        } // while itPath.hasNext()

        int[] condMeans = new int[4];
        condMeans[0] = (int) Math.round((double) nIns / (double) pathList.size());
        condMeans[1] = (int) Math.round((double) nDel / (double) pathList.size());
        condMeans[2] = (int) Math.round((double) nInsTotal / (double) pathList.size());
        condMeans[3] = (int) Math.round((double) nDelTotal / (double) pathList.size());

        return condMeans;
    }

    /* ********************************************************************** *
     *                      M I S C E L L A N E O U S                         *
     * ********************************************************************** */
    static private MetabolicNetwork loadNetwork(String dirName, String fileName) {
        String[][] data = Utilities.readList2(dirName, fileName, false);

        String networkId = data[0][0];

        // extract the nodes
        String[][] nodes = new String[data.length - 1][3];
        // first (row, col) entry is the network name
        for (int i = 1; i < data.length; i++) {
            String node = data[i][0];
            nodes[i - 1][0] = node;
            nodes[i - 1][1] = node;
            nodes[i - 1][2] = node;
        }

        // process the network column-wise
        String[] edgeRow = new String[data[0].length - 1]; // first column contains the node names
        ArrayList lstEdges = new ArrayList(data.length - 1);
        for (int j = 1; j <= edgeRow.length; j++) {
            String rxnId = data[0][j];
            boolean reactionPresent = false;
            String substrate = "", product = "";
            for (int i = 1; i < data.length; i++) {
                Integer value = new Integer(data[i][j]);
                if (value.intValue() != 0) {
                    reactionPresent = true;
                    if (value.intValue() < 0) {// substrate
                        substrate += " + " + data[i][0];
                    } else {// product
                        product += " + " + data[i][0];
                    }
                } // if value != 0
            } // for each row

            if (reactionPresent) {
                // there is extra "+" in substrates and products
                substrate = substrate.trim().substring(1).trim();
                product = product.trim().substring(1).trim();
                String[] edge = new String[]{rxnId, rxnId, substrate + " => " + product};
                lstEdges.add(edge);
            }
        } // for each column

        // extract the edges
        String[][] edges = new String[lstEdges.size()][4];
        Iterator it = lstEdges.iterator();
        int idx = -1;
        while (it.hasNext()) {
            String[] edge = (String[]) it.next();
            idx++;

            edges[idx][0] = edge[0];
            edges[idx][1] = edge[1];
            edges[idx][2] = edge[2];
            edges[idx][3] = "N"; // irreversible
        } // while it.hasNext()

        // create dummy pathway
        String[][] pathways = new String[][]{{"P1", "P1"}};
        // create dummy rxn pathway
        String[][] rxnPathways = new String[edges.length][2];
        for (int i = 0; i < rxnPathways.length; i++) {
            rxnPathways[i][0] = edges[i][0];
            rxnPathways[i][1] = "P1";
        }

        return MetabolicNetwork.buildMetabolicNetwork(networkId, networkId, new ArrayList(), nodes, edges, rxnPathways, null, pathways, null, null, null, null);
    }

    static private int[] getNeighboursCount(String dirName, String refNWFilename,
            String nwFilename, boolean printOutput) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        refNetwork.removeInactiveMetabolites();
        refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        MetabolicNetwork theNetwork = null;
        if (nwFilename != null && !nwFilename.isEmpty()) {
            // start network
            theNetwork = loadNetwork(dirName, nwFilename);
            theNetwork.addInactiveReactions(refNetwork.getReactions());
            theNetwork.removeInactiveMetabolites();
            theNetwork.setupNetworkMatrices();
            theNetwork.setupNetworkSequence();
        }

        int[] neighboursCount = new int[refNetwork.getDirectedReactions().size()];
        if (theNetwork == null) {
            // neighbours count for reference network
            neighboursCount = refNetwork.getNeighboursCount();
        } else {
            // the neighbours count for the given network
            neighboursCount = theNetwork.getNeighboursCount(refNetwork);
        }

        if (printOutput) {
            System.out.println("Neighbours Count:");
            for (int i = 0; i < neighboursCount.length; i++) {
                System.out.print(((DirectedReaction) refNetwork.getDirectedReactions().get(i)).getID() + "\t");
                System.out.println(neighboursCount[i]);
            }
        }
        return neighboursCount;
    }

    static private int[] getNeighboursCount(String[] orgIds, String[] pathwayIds, boolean printOutput) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String ignoreDir;
        try {
            ignoreDir = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            ignoreDir = null;
        }
        String ignoreFile = "ignoreMetabolites.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        // remove current metabolites
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);

        // create the metabolic network
        MetabolicNetwork orgNetwork;
        // network
        if (orgIds == null || orgIds.length == 0)
            orgNetwork = refNetwork;
        else {
            String[][] arrOrg = database.getOrganisms(orgIds);
            ArrayList orgList = Organism.buildOrganismList(arrOrg);

            String[][] orgReactions = database.getOrganismReactions(orgIds, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIds, pathwayIds);

            orgNetwork = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
                    Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);
            orgNetwork.removeInactiveMetabolites();
            // print reaction connectivity
            System.out.println(orgNetwork.getName() + ":");
            // remove current metabolties
            orgNetwork.getMetabolites().removeAll(ignoreMetabolites);
        }
        
        refNetwork.setupNetworkMatrices();
        orgNetwork.addInactiveReactions(refNetwork.getReactions());
        orgNetwork.setupNetworkMatrices();
        int[] neighboursCount = orgNetwork.getNeighboursCount(refNetwork);
        
        if (printOutput) {
            System.out.println("Neighbours Count:");
            Iterator it = refNetwork.getDirectedReactions().iterator();
            while (it.hasNext()) {
                DirectedReaction theReaction = (DirectedReaction) it.next();
                int idx = orgNetwork.getDirectedReactions().indexOf(theReaction);
                int nc = 0;
                if (idx != -1)
                    nc = neighboursCount[idx];
                System.out.println(theReaction.getID() + "\t" + Integer.toString(nc));
            }
        }
        
        return neighboursCount;        
    }

    static private void printMetabolicNetworks(String[] orgIds, String[] pathwayIds,
            String[] orgIdsForCoreAndProhib, boolean writeToFile) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null) {
            pathwayIds = Utilities.extractColumn(keggPathways, 0);
        }

        String ignoreDir;
        try {
            ignoreDir = new File(".").getCanonicalPath();
            String ignoreFile = "ignoremetabolites.txt";
            ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

            String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
            MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                    "Reference", new ArrayList(), database.getMetabolites(""),
                    database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                    null, null, null, null);
            refNetwork.removeInactiveMetabolites();
            refNetwork.getMetabolites().removeAll(ignoreMetabolites);
            refNetwork.setupNetworkMatrices();

            System.out.println("Reference Network:");
            BufferedWriter outRef;
            if (writeToFile) {
                String filename = "ref_" + Utilities.toString(pathwayIds, "_") + ".txt";
                outRef = new BufferedWriter(new FileWriter(new File(filename)));
            } else
                outRef = null;
            HashSet ignoreReactions = refNetwork.printNetwork(outRef);
            outRef.close();

            // create the metabolic networks
            MetabolicNetwork orgNetwork;
            // network 1
            if (orgIds != null && orgIds.length != 0) {
                for (int i = 0; i < orgIds.length; i++) {

                    String[][] arrOrg = database.getOrganisms(new String[]{orgIds[i]});
                    ArrayList orgList = Organism.buildOrganismList(arrOrg);

                    String[][] orgReactions = database.getOrganismReactions(new String[]{orgIds[i]}, pathwayIds);
                    String[][] rxnEnzymes = database.getReactionEnzymes(new String[]{orgIds[i]}, pathwayIds);

                    orgNetwork = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
                            Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
                            database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                            null, null, null, rxnEnzymes);

                    orgNetwork.addInactiveReactions(refNetwork.getReactions());
                    orgNetwork.removeInactiveMetabolites();
                    orgNetwork.getMetabolites().removeAll(ignoreMetabolites);

                    orgNetwork.setupNetworkMatrices();

                    System.out.println(orgNetwork.getName() + ":");

                    BufferedWriter out;
                    if (writeToFile) {
                        String filename = orgNetwork.getID() + "_" + Utilities.toString(pathwayIds, "_") + ".txt";
                        out = new BufferedWriter(new FileWriter(new File(filename)));
                    } else
                        out = null;

                    orgNetwork.printNetwork(out, ignoreReactions);
                    out.close();
                } // end for each organism
            }


            // get core and prohibted edges
            ArrayList coreEdges, prohibEdges;
            if (orgIdsForCoreAndProhib != null & orgIdsForCoreAndProhib.length != 0) {
//            String[][] arrOrg = database.getOrganisms(orgIdsForCoreAndProhib);
                ArrayList orgList = Organism.buildOrganismList(orgIdsForCoreAndProhib);

                String[][] orgReactions = database.getOrganismReactions(orgIdsForCoreAndProhib, pathwayIds);
                String[][] rxnEnzymes = database.getReactionEnzymes(orgIdsForCoreAndProhib, pathwayIds);

                // CORE EDGES
                MetabolicNetwork coreNetwork = MetabolicNetwork.buildMetabolicNetwork("core",
                        "Core Network", orgList, database.getMetabolites(""),
                        database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                        null, null, null, rxnEnzymes);
                // delete the reactions not present in all organisms
                coreNetwork.deleteReactionsNotPresentInAllOrganisms();
                // setup hyperedges
                coreNetwork.setupNetworkSequence();
                // get the core edges
                coreEdges = coreNetwork.getDirectedReactions();

                System.out.println("Core Edges:");
                BufferedWriter outCore;
                if (writeToFile) {
                    String filename = "core_" + Utilities.toString(pathwayIds, "_") + ".txt";
                    outCore = new BufferedWriter(new FileWriter(new File(filename)));
                } else
                    outCore = null;

                Iterator itCore = coreEdges.iterator();
                while (itCore.hasNext()) {
                    DirectedReaction theReaction = (DirectedReaction) itCore.next();

                    if (ignoreReactions.contains(theReaction))
                        continue;
                    if (outCore == null)
                        System.out.println(theReaction.getID());
                    else
                        outCore.write(theReaction.getID() + "\r\n");
                }
                outCore.close();

                // PROHIBITED EDGES
                MetabolicNetwork unionNetwork = MetabolicNetwork.buildMetabolicNetwork("prohib",
                        "Prohibited Network", orgList, database.getMetabolites(""),
                        database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                        null, null, null, rxnEnzymes);
                // setup hyperedges
                unionNetwork.setupNetworkSequence();
                // get the prohibited edges
                prohibEdges = refNetwork.getDirectedReactions();
                prohibEdges.removeAll(unionNetwork.getDirectedReactions());

                System.out.println("Prohibited Edges:");
                BufferedWriter outProhib;
                if (writeToFile) {
                    String filename = "prohibited_" + Utilities.toString(pathwayIds, "_") + ".txt";
                    outProhib = new BufferedWriter(new FileWriter(new File(filename)));
                } else
                    outProhib = null;
                Iterator itProhib = prohibEdges.iterator();
                while (itProhib.hasNext()) {
                    DirectedReaction theReaction = (DirectedReaction) itProhib.next();

                    if (ignoreReactions.contains(theReaction))
                        continue;

                    if (outProhib == null)
                        System.out.println(theReaction.getID());
                    else
                        outProhib.write(theReaction.getID() + "\r\n");
                }
                outProhib.close();
            }
        } catch (Exception ex) {
            ignoreDir = null;
        }

    }

    static private void printReferenceReactionsAndStatus(String[] pathwayIds,
            String[] orgIdsForCoreAndProhib, boolean writeToFile) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null) {
            pathwayIds = Utilities.extractColumn(keggPathways, 0);
        }

        String ignoreDir;
        try {
            ignoreDir = new File(".").getCanonicalPath();
            String ignoreFile = "ignoremetabolites.txt";
            ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

            String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
            MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                    "Reference", new ArrayList(), database.getMetabolites(""),
                    database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                    null, null, null, null);
            refNetwork.removeInactiveMetabolites();
            refNetwork.getMetabolites().removeAll(ignoreMetabolites);
            refNetwork.setupNetworkMatrices();

            System.out.println("Reference Network:");
            HashSet ignoreReactions = refNetwork.printNetwork(null);

            // get core and prohibted edges
            ArrayList coreEdges, prohibEdges;
            if (orgIdsForCoreAndProhib != null & orgIdsForCoreAndProhib.length != 0) {
//            String[][] arrOrg = database.getOrganisms(orgIdsForCoreAndProhib);
                ArrayList orgList = Organism.buildOrganismList(orgIdsForCoreAndProhib);

                String[][] orgReactions = database.getOrganismReactions(orgIdsForCoreAndProhib, pathwayIds);
                String[][] rxnEnzymes = database.getReactionEnzymes(orgIdsForCoreAndProhib, pathwayIds);

                // CORE EDGES
                MetabolicNetwork coreNetwork = MetabolicNetwork.buildMetabolicNetwork("core",
                        "Core Network", orgList, database.getMetabolites(""),
                        database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                        null, null, null, rxnEnzymes);
                // delete the reactions not present in all organisms
                coreNetwork.deleteReactionsNotPresentInAllOrganisms();
                // setup hyperedges
                coreNetwork.setupNetworkSequence();
                // get the core edges
                coreEdges = coreNetwork.getDirectedReactions();

                // PROHIBITED EDGES
                MetabolicNetwork unionNetwork = MetabolicNetwork.buildMetabolicNetwork("prohib",
                        "Prohibited Network", orgList, database.getMetabolites(""),
                        database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                        null, null, null, rxnEnzymes);
                // setup hyperedges
                unionNetwork.setupNetworkSequence();
                // get the prohibited edges
                prohibEdges = refNetwork.clone().getDirectedReactions();
                prohibEdges.removeAll(unionNetwork.getDirectedReactions());

            } else {
                coreEdges = new ArrayList();
                prohibEdges = new ArrayList();
            }

            BufferedWriter outRef;
            if (writeToFile) {
                String filename = "rxn" + Utilities.toString(pathwayIds, "_") + ".txt";
                outRef = new BufferedWriter(new FileWriter(new File(filename)));
            } else
                outRef = null;

            Iterator it = refNetwork.getDirectedReactions().iterator();
            while (it.hasNext()) {
                DirectedReaction theReaction = (DirectedReaction) it.next();

                if (ignoreReactions.contains(theReaction))
                    continue;

                String str = theReaction.getID();
                str += "\t" + (coreEdges.contains(theReaction) ? "1" : "0");
                str += "\t" + (prohibEdges.contains(theReaction) ? "1" : "0") + "\r\n";
                if (outRef != null)
                    outRef.write(str);
                else
                    System.out.print(str);
            }
            outRef.close();

        } catch (Exception ex) {
            ignoreDir = null;
        }

    }

    /* ********************************************************************** *
     *                              U N U S E D                               *
     * ********************************************************************** */
    static private void printPathProbabilities(String[] pathwayIds, String dataFilename) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.setupNetworkSequence();

        ArrayList pathList, networkList;
        BigDecimal[] pathProbList;
        try {
            ObjectInputStream in =
                    new ObjectInputStream(
                    new FileInputStream(dataFilename));
            pathList = (ArrayList) in.readObject();
            networkList = (ArrayList) in.readObject();
            pathProbList = (BigDecimal[]) in.readObject();
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        NetworkEvolver.summariseEvolution(refNetwork, networkList);

        int idx = 0;
        Iterator it = pathList.iterator();
        while (it.hasNext()) {
            // get the next path
            int[] path = (int[]) it.next();
            // convert into edge indices
            int[] edgeIndices = new int[path.length];
            for (int i = 0; i < edgeIndices.length; i++)
                edgeIndices[i] = path[i] - 1;

            System.out.print(Utilities.toString(refNetwork.getDirectedReactionIds(edgeIndices), ", ") + "\t");
            System.out.println(pathProbList[idx++].round(MathContext.DECIMAL128).toString());

        }
    }

    static private void printSummaries(String[] orgIds1, String[] pathwayIds, String dataFilename) {

        /*        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();
        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
        pathwayIds = Utilities.extractColumn(keggPathways, 0);
        String[][] rxnPathways = database.getReactionPathways(pathwayIds);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
        "Reference", new ArrayList(), database.getMetabolites(""),
        database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
        null, null, null, null);
        refNetwork.setupNetworkSequence();
        String ignoreDir = "E:/MyData/Oxford/DPhil/Code/KEGG DB/DataFiles";
        String ignoreFile = "ignoreCompoundsTCA.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));
        // create the metabolic networks
        MetabolicNetwork orgNetwork1;
        // network 1
        if (orgIds1 == null || orgIds1.length == 0)
        orgNetwork1 = refNetwork;
        else {
        String[][] arrOrg = database.getOrganisms(orgIds1);
        ArrayList orgList = Organism.buildOrganismList(arrOrg);
        String[][] orgReactions = database.getOrganismReactions(orgIds1, pathwayIds);
        String[][] rxnEnzymes = database.getReactionEnzymes(orgIds1, pathwayIds);
        orgNetwork1 = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
        Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
        database.getReactions(pathwayIds), null, orgReactions, keggPathways,
        null, null, null, rxnEnzymes);
        orgNetwork1.addInactiveReactions(refNetwork.getReactions());
        orgNetwork1.removeInactiveMetabolites();
        orgNetwork1.getMetabolites().removeAll(ignoreMetabolites);
        //orgNetwork1.setupNetworkMatrices();
        orgNetwork1.setupNetworkSequence();
        }
         */
        ArrayList pathList, networkList;
        double[] insRates1, delRates1;
        BigDecimal[] pathProbList;
        int nBurning = 10000;
        try {
            ObjectInputStream in =
                    new ObjectInputStream(
                    new FileInputStream(dataFilename));
            pathList = (ArrayList) in.readObject();
//            networkList = (ArrayList) in.readObject();
            in.readObject(); // we dont need the network list
            pathProbList = (BigDecimal[]) in.readObject();
            insRates1 = (double[]) in.readObject();
            delRates1 = (double[]) in.readObject();
            Double[] insRates = new Double[insRates1.length];
            Double[] delRates = new Double[delRates1.length];
            for (int i = 0; i < insRates.length; i++)
                insRates[i] = new Double(insRates1[i]);
            for (int i = 0; i < delRates.length; i++)
                delRates[i] = new Double(delRates1[i]);

//            NetworkEvolver.summariseEvolution(orgNetwork1, networkList);
            System.out.println(new Date());
//            NetworkEvolver.pathLengthDistribution(pathList, nBurning, true);
//            HashMap pathDistribution = NetworkEvolver.pathDistribution(pathList, nBurning, true);
//            HashMap networkDistribution = NetworkEvolver.networkDistribution(orgNetwork1, new ArrayList(pathDistribution.keySet()), true);
//            HashMap insRateDistribution = NetworkEvolver.rateDistribution(insRates, nBurning, true, "Insertion Rate:");
//            HashMap delRateDistribution = NetworkEvolver.rateDistribution(delRates, nBurning, true, "Deletion Rate:");
            HashMap likelihood = new NetworkMCMC().calculateLikelihood(pathList, pathProbList, insRates, delRates, nBurning, true);
            int block = 100000, start = 0, end = block;
            while (end <= pathList.size()) {
//                System.out.println("Size:" + Integer.toString(end));
//                HashMap likelihood = NetworkEvolver.calculateLikelihood(new ArrayList(pathList.subList(start, end)), pathProbList, insRates, delRates, nBurning, true);
                end += block;
            }


            in.close();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

    }

    static private void studyKOCascading(int evolutionModel, String[] orgIds, String[] pathwayIds) {

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String ignoreDir;
        try {
            ignoreDir = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            ignoreDir = null;
        }
        String ignoreFile = "ignoremetabolites.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);

        // create the metabolic network
        MetabolicNetwork orgNetwork;
        if (orgIds == null || orgIds.length == 0)
            orgNetwork = refNetwork;
        else {
            String[][] arrOrg = database.getOrganisms(orgIds);
            ArrayList orgList = Organism.buildOrganismList(arrOrg);

            String[][] orgReactions = database.getOrganismReactions(orgIds, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIds, pathwayIds);

            orgNetwork = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
                    Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);

//            orgNetwork.addInactiveReactions(refNetwork.getReactions());
            orgNetwork.removeInactiveMetabolites();
            orgNetwork.getMetabolites().removeAll(ignoreMetabolites);

            if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
                orgNetwork.setupNetworkMatrices();

            orgNetwork.setupNetworkSequence();
        }


        // MCMC parameters
        double insRate = 0.05, delRate = 0.03;
        double dependenceProbability = 0.5;
        double evolTime = 1.0;
        ArrayList coreEdges = new ArrayList();
        ArrayList prohibEdges = new ArrayList();

        for (int i = 0; i < 100; i++) {
            new NetworkSimulator().simulateKOCascading(orgNetwork, coreEdges, prohibEdges, evolutionModel, insRate, delRate, ignoreMetabolites, dependenceProbability);
        }
    }

    static private void networkProperties(String[] orgIds, String[] pathwayIds) {

        HashMap rxnConnectivity, rxnClustCoeff;

        Database database = new Database(DB_URL, DB_NAME, DB_USER, DB_PASSWORD, DB_DRIVER);
        database.readApplicationParameters();

        // read the kegg pathways from the db
        String[][] keggPathways = database.getKEGGPathways(pathwayIds);
        if (pathwayIds == null)
            pathwayIds = Utilities.extractColumn(keggPathways, 0);

        String ignoreDir;
        try {
            ignoreDir = new File(".").getCanonicalPath();
        } catch (Exception ex) {
            ignoreDir = null;
        }
        String ignoreFile = "ignoreMetabolites.txt";
        ArrayList ignoreMetabolites = Metabolite.buildList(Utilities.readList2(ignoreDir, ignoreFile, false));

        String[][] rxnPathways = database.getPathwayReactions(pathwayIds, false, false);
        MetabolicNetwork refNetwork = MetabolicNetwork.buildMetabolicNetwork("ref",
                "Reference", new ArrayList(), database.getMetabolites(""),
                database.getReactions(pathwayIds), rxnPathways, null, keggPathways,
                null, null, null, null);
        refNetwork.removeInactiveMetabolites();
        // print reaction connectivity
        System.out.println("Reference Network:");
        rxnConnectivity = refNetwork.calculateReactionConnectivity(false);
        printConnectivitySummary(rxnConnectivity);
        rxnClustCoeff = refNetwork.calculateReactionClusteringCoefficient(false);
        printClusteringCoefficientSummary(rxnClustCoeff);
        // remove current metabolites
        refNetwork.getMetabolites().removeAll(ignoreMetabolites);
        // print reaction connectivity
        System.out.println("Reference Network - Without Current Metabolites:");
        rxnConnectivity = refNetwork.calculateReactionConnectivity(false);
        printConnectivitySummary(rxnConnectivity);
        rxnClustCoeff = refNetwork.calculateReactionClusteringCoefficient(false);
        printClusteringCoefficientSummary(rxnClustCoeff);


        // create the metabolic network
        MetabolicNetwork orgNetwork;
        // network
        if (orgIds == null || orgIds.length == 0)
            orgNetwork = refNetwork;
        else {
            String[][] arrOrg = database.getOrganisms(orgIds);
            ArrayList orgList = Organism.buildOrganismList(arrOrg);

            String[][] orgReactions = database.getOrganismReactions(orgIds, pathwayIds);
            String[][] rxnEnzymes = database.getReactionEnzymes(orgIds, pathwayIds);

            orgNetwork = MetabolicNetwork.buildMetabolicNetwork(Utilities.toString(NetworkObject.getIDs(orgList), ","),
                    Utilities.toString(NetworkObject.getNames(orgList), ","), orgList, database.getMetabolites(""),
                    database.getReactions(pathwayIds), null, orgReactions, keggPathways,
                    null, null, null, rxnEnzymes);
            orgNetwork.removeInactiveMetabolites();
            // print reaction connectivity
            System.out.println(orgNetwork.getName() + ":");
            rxnConnectivity = orgNetwork.calculateReactionConnectivity(false);
            printConnectivitySummary(rxnConnectivity);
            rxnClustCoeff = orgNetwork.calculateReactionClusteringCoefficient(false);
            printClusteringCoefficientSummary(rxnClustCoeff);
            // remove current metabolties
            orgNetwork.getMetabolites().removeAll(ignoreMetabolites);
            // print reaction connectivity
            System.out.println(orgNetwork.getName() + " - Without Current Metabolites:");
            rxnConnectivity = orgNetwork.calculateReactionConnectivity(false);
            printConnectivitySummary(rxnConnectivity);
            rxnClustCoeff = orgNetwork.calculateReactionClusteringCoefficient(false);
            printClusteringCoefficientSummary(rxnClustCoeff);

            refNetwork.setupNetworkMatrices();
            orgNetwork.addInactiveReactions(refNetwork.getReactions());
            orgNetwork.setupNetworkMatrices();
            int[] neighboursCount = orgNetwork.getNeighboursCount(refNetwork);
            System.out.println("Neighbours Count:");
            Iterator it = refNetwork.getDirectedReactions().iterator();
            while (it.hasNext()) {
                DirectedReaction theReaction = (DirectedReaction) it.next();
                int idx = orgNetwork.getDirectedReactions().indexOf(theReaction);
                int nc = 0;
                if (idx != -1)
                    nc = neighboursCount[idx];
                System.out.println(theReaction.getID() + "\t" + Integer.toString(nc));
            }

        }
    }

    static private void printConnectivitySummary(HashMap rxnConnectivity) {
        Iterator it = rxnConnectivity.values().iterator();
        Integer[] maxDegree = new Integer[]{0, 0};
        Double[] avgDegree = new Double[]{0.0, 0.0};
        while (it.hasNext()) {
            Integer[] degree = (Integer[]) it.next();
            // set the maximum value
            maxDegree[0] = Math.max(degree[0], maxDegree[0]);
            maxDegree[1] = Math.max(degree[1], maxDegree[1]);

            // add to the total (for average)
            avgDegree[0] = avgDegree[0] + degree[0];
            avgDegree[1] = avgDegree[1] + degree[1];
        }
        avgDegree[0] /= rxnConnectivity.size();
        avgDegree[1] /= rxnConnectivity.size();

        System.out.println("Average Connectivity: " + avgDegree[0] + "\t" + avgDegree[1]);
        System.out.println("Maximum Connectivity: " + maxDegree[0] + "\t" + maxDegree[1]);
    }

    static private void printClusteringCoefficientSummary(HashMap rxnClustCoeff) {
        Iterator it = rxnClustCoeff.values().iterator();
        Double maxClustCoeff = 0.0;
        Double avgClustCoeff = 0.0;
        while (it.hasNext()) {
            Double clustCoeff = (Double) it.next();
            // set the maximum value
            maxClustCoeff = Math.max(clustCoeff, maxClustCoeff);

            // add to the total (for average)
            avgClustCoeff = avgClustCoeff + clustCoeff;
        }
        avgClustCoeff /= rxnClustCoeff.size();

        System.out.println("Average Clustering Coefficieint: " + avgClustCoeff);
        System.out.println("Maximum Clustering Coefficieint: " + maxClustCoeff);
    }

    static private void studyNetworkCharacteristics(int evolutionModel, String dirName,
            String refNWFilename, String nwFilename) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        // start network
        MetabolicNetwork theNetwork = loadNetwork(dirName, nwFilename);
        theNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            theNetwork.setupNetworkMatrices();
        theNetwork.setupNetworkSequence();

        System.out.println("Starting Network:");
        System.out.println(Utilities.toString(theNetwork.getNeighboursCount(), "\t"));
        int nIter = 50000;
        int updateInterval = 5000;
        double insRate = 0.9;
        double delRate = 0.1;
        double dependenceProbability = 0.5;
        ArrayList lstNetworks = new NetworkSimulator().simulateNetworkEvolution(theNetwork, evolutionModel, nIter, new ArrayList(), new ArrayList(), insRate, delRate, updateInterval, dependenceProbability);
        int sampleSize = 5000;
        int start = lstNetworks.size() - sampleSize - 1;
        Iterator it = lstNetworks.iterator();
        int idx = 0;
        int[] neighbourCount = new int[theNetwork.getReactionSequence().length];
        int[] edgeCount = new int[theNetwork.getReactionSequence().length];
        while (it.hasNext()) {
            Byte[] nwSeq = (Byte[]) it.next();
            if (idx++ < start)
                continue;

            theNetwork = theNetwork.getNetwork(nwSeq);
            int[] neighbours = theNetwork.getNeighboursCount();
            for (int i = 0; i < nwSeq.length; i++) {
                if (nwSeq[i].equals(MetabolicNetwork.SEQ_ENTRY_PRESENT)) {
                    edgeCount[i]++;
                    neighbourCount[i] += neighbours[i];
                }
            }
        }

        for (int i = 0; i < neighbourCount.length; i++) {
            if (edgeCount[i] > 0)
                neighbourCount[i] = Math.round((float) neighbourCount[i] / (float) edgeCount[i]);
            else
                neighbourCount[i] = -1;
        }

        System.out.println(Utilities.toString(neighbourCount, "\t"));
    }

    static private void exponentiation_ToyNetworks(int evolutionModel, String dirName,
            String refNWFilename, String startNWFilename, String endNWFilename) {

        // load the metabolic networks
        MetabolicNetwork refNetwork = loadNetwork(dirName, refNWFilename);
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            refNetwork.setupNetworkMatrices();
        refNetwork.setupNetworkSequence();
        // start network
        MetabolicNetwork startNetwork = loadNetwork(dirName, startNWFilename);
        startNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            startNetwork.setupNetworkMatrices();
        startNetwork.setupNetworkSequence();

        // end network
        MetabolicNetwork endNetwork = loadNetwork(dirName, endNWFilename);
        endNetwork.addInactiveReactions(refNetwork.getReactions());
        if (evolutionModel != NetworkEvolver.EVOL_MODEL_INDEPENDENT_EDGE)
            endNetwork.setupNetworkMatrices();
        endNetwork.setupNetworkSequence();

        BigDecimal insRate = new BigDecimal("0.05"), delRate = new BigDecimal("0.03");
        double evolTime = 1.0;
        double dependenceProbability = 0.5;

        String strStartNetwork = Utilities.toString(startNetwork.getReactionSequence(), "");
        String strEndNetwork = Utilities.toString(endNetwork.getReactionSequence(), "");
        int startIdx = Integer.parseInt(new StringBuffer(strStartNetwork).reverse().toString(), 2);
        int endIdx = Integer.parseInt(new StringBuffer(strEndNetwork).reverse().toString(), 2);

        int edgeCount = refNetwork.getDirectedReactions().size();
        ArrayList lstNetworks = new NetworkEvolver().enumerateAllNetworks(edgeCount);
        //ArrayList lstNetworks = new NetworkEvolver().enumerateNetworks(theNetwork, endNetwork, new ArrayList(), new ArrayList(), 0);
        //startIdx = 0;
        //endIdx = 1;

        System.out.println("Networks enumerated. " + new Date().toString());
//        BigDecimal[][] rateMatrix = new NetworkEvolver().generateRateMatrix(lstNetworks, 
//                refNetwork, NetworkEvolver.EVOL_MODEL_NEIGHBOUR_DEPENDENT, insRate, delRate);
//        System.out.println("Rate matrix calculated. " + new Date().toString());
//        System.out.println(Utilities.toString(rateMatrix, " "));
//        BigDecimal[][] transProb = Utilities.exponentiate(rateMatrix, evolTime, 10, 301);
//        System.out.println("Rate matrix exponentiated. " + new Date().toString());
//        System.out.println(transProb[startIdx][endIdx].toString());

        double[][] rateMatrix = null;
        String rmFilename = "";
        if (rmFilename.isEmpty()) {
            rateMatrix = new NetworkEvolver().generateRateMatrix(lstNetworks,
                    refNetwork, evolutionModel, insRate.doubleValue(), delRate.doubleValue(), dependenceProbability);

            if (1 == 1) {
                try {
                    String filename = "toy_rateMatrix_" + Integer.toString(lstNetworks.size()) + ".txt";
                    // create the output file
                    BufferedWriter out = new BufferedWriter(new FileWriter(new File(filename)));
                    for (int i = 0; i < rateMatrix.length; i++) {
                        double[] ds = rateMatrix[i];
                        for (int j = 0; j < ds.length - 1; j++) {
                            out.write(new BigDecimal(ds[j]).round(new MathContext(6)).toString() + "\t");
                        }
                        // write the last element
                        out.write(new BigDecimal(ds[ds.length - 1]).round(new MathContext(6)).toString());
                        out.write("\r\n");

                    }
                    out.close(); // Also flushes output
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        } else {
            rateMatrix = new double[lstNetworks.size()][lstNetworks.size()];
            try {
                BufferedReader in = new BufferedReader(new FileReader(new File(rmFilename)));
                String s;
                int i = -1;
                while ((s = in.readLine()) != null) {
                    String[] row = s.split("\t");
                    i++;

                    for (int j = 0; j < row.length; j++) {
                        rateMatrix[i][j] = Double.parseDouble(row[j]);
                    }
                }
                in.close(); // Also flushes output
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

//        double[][] transProb = Utilities.exponentiate(rateMatrix, evolTime, 20);
//        System.out.println(transProb[startIdx][endIdx]);
    }
}

