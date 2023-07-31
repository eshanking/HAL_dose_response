package DoseResponse;

import java.io.File;
import java.io.IOException;

import HAL.lib.CommandLine;

@CommandLine.Command(name = "2-D Agent-Based Model of Bacterial Evolution with Drug Diffusion",
        mixinStandardHelpOptions = true,
        showDefaultValues = true,
        description = "A 2D on-lattice agent-based model that simulates bacterial cells evolving with dose-dependent growth rates and drug diffusion.")


public class DoseResponseModel implements Runnable {
    DoseResponseGrid model = new DoseResponseGrid();
    @CommandLine.Option(names = { "-s", "--seed"}, description="Random number seed.") 
    int seed = model.seed;
    // ------------------------- Experimental Setup -------------------------
    @CommandLine.Option(names = { "-xDim", "--xDim"}, description="x-dimension of domain (in lattice sites)") 
    int xDim = model.xDim;
    @CommandLine.Option(names = { "-yDim", "--yDim"}, description="y-dimension of domain (in lattice sites)") 
    int yDim = model.yDim;
    @CommandLine.Option(names = { "-n", "--nReplicates"}, description="Number of replicates.") 
    int nReplicates = model.nReplicates;
    @CommandLine.Option(names = { "-w0", "--initialWidth"}, description="Initial population width") 
    int initWidth = model.initWidth;
    @CommandLine.Option(names = { "-g0", "--initialGeometry"}, description="Initial population structure geometry") 
    String initGeometry = model.initGeometry;
    @CommandLine.Option(names = { "-d0", "--initialDensity"}, description="Initial population density") 
    double initDensity = model.initDensity;
    @CommandLine.Option(names = { "-pM", "--pMutant"}, description="Initial mutant proportion") 
    double initMutantProp = model.initMutantProp;

    @CommandLine.Option(names = { "-dt", "--dt"}, description="Time step in days.") 
    double dt = model.dt;
    @CommandLine.Option(names = { "-nTSteps", "--nTSteps"}, description="Number of time steps to run for.")
    int nTSteps = model.nTSteps;

    // ------------------------- Cell Properties -------------------------
    @CommandLine.Option(names = { "-dieProb", "--dieProb"}, description="Probability of cell death")
    double dieProb = model.dieProb;
    @CommandLine.Option(names = { "-mutProb", "--mutProb"}, description="Probability of cell mutation")
    double mutProb = model.mutProb;
    // ------------------------- Diffusion Grid Properties -------------------------
    @CommandLine.Option(names = { "-diffRate", "--diffRate"}, description="Diffusion rate of drug")
    double diffRate = model.diffRate;
    @CommandLine.Option(names = { "-srcConc", "--srcConc"}, description="Drug concentration at blood vessel source")
    double[] srcConc = model.srcConc;
    @CommandLine.Option(names = {"-consumpRate", "--consumpRate"}, description="Field consumption rate")
    double consumpRate = model.consumpRate;
    // ------------------------- Output - Text -------------------------
    @CommandLine.Option(names = { "--outDir"}, description="Directory which to save output files to.") 
    String outDir = "./tmp/";
    @CommandLine.Option(names = { "--imageOutDir"}, description="Directory which to save images to.") 
    String imageOutDir = model.imageOutDir;
    // ------------------------- Output - Visualisation -------------------------
    @CommandLine.Option(names = { "--imageFrequency"}, description="Frequency at which an image of the tumour is saved. Negative number turns it off.") 
    int imageFrequency = model.imageFrequency;
    @CommandLine.Option(names = { "-visualiseB", "--visualiseB"}, description="Whether or not to show visualization.")
    Boolean visualiseB = model.visualiseB;
    @CommandLine.Option(names = { "--saveModelState"}, description="Whether or not to save the model object at the end of the simulation.") 
    Boolean saveModelState = model.saveModelState;
    @CommandLine.Option(names = { "--saveFinalDiffImg"}, description="Save final image of diffusion grid.") 
    Boolean saveFinalDiffImg = model.saveFinalDiffImg;
    @CommandLine.Option(names = { "--saveFinalDiffGrid"}, description="Save final image of diffusion grid.") 
    Boolean saveFinalDiffGrid = model.saveFinalDiffGrid;
    @CommandLine.Option(names = { "--saveFinalPopGrid"}, description="Save final population grid") 
    Boolean saveFinalPopGrid = model.saveFinalPopGrid;

    // ------------------------------------------------------------------------------------------------------------
    public void run(){
        DoseResponseGrid model;
        String outFName;
        String diffFileName;
        String popGridFileName;

        new File(outDir).mkdirs(); // create output directory if it doesn't exist

        int[] replicateIdList = new int[] {seed};
        int replicateId;
        if (nReplicates>1) {
            replicateIdList = new int[nReplicates];
            for (int i=0; i<nReplicates; i++) {replicateIdList[i]=i;}
        }

        for (int replicateIdx = 0; replicateIdx < nReplicates; replicateIdx++) {
            replicateId = replicateIdList[replicateIdx];
            
            // System.out.println(xDim);
            model = new DoseResponseGrid(xDim, yDim);

            // Set the random number seed. Behaviour depends no whether this is a single run or part of a series of nReplicate runs. By default will assign every replicate the value ```seed=replicateId```
            if (seed != -1) {
                if (nReplicates == 1) {
                    model.SetSeed(seed);
                } else {
                    model.SetSeed(replicateId);
                }
            } else {
                model.SetSeed(replicateId);
            }

            outFName = outDir + "RepId_" + replicateId + ".csv";
            diffFileName = outDir + "RepId_" + replicateId + "_diffGrid.csv";
            popGridFileName = outDir + "RepId_" + replicateId + "_popGrid.csv";

            model.SetSaveParams(diffFileName, outFName, popGridFileName, saveFinalDiffGrid, saveFinalPopGrid);

            if (imageFrequency > 0) {
                model.visualiseB = true;
                String currImgOutDir = imageOutDir + "RepId_" + replicateId;
                new File(currImgOutDir).mkdirs();
                model.ConfigureImaging(visualiseB, imageOutDir, imageFrequency,
                                    saveFinalDiffImg);
                                    // scaleFactor, pause, colorList);
            }

            // Set parameters
            
            model.InitialiseCellLog(outFName);
            model.SetDiffParams(srcConc, diffRate, consumpRate);
            model.SetCellParams(mutProb, dieProb);
            // System.out.println("initWidth: " + initWidth);
            model.SetInitPopParams(initGeometry, initWidth,initDensity,initMutantProp);
            model.ConfigureExperiment(nReplicates, nTSteps, dt);

            try {
                model.Run();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            model.Close();
        }
    }
    public static void main(String[] args) {
        int exitCode = new CommandLine(new DoseResponseModel()).execute(args); 
        System.exit(exitCode);
    }
        
} 

