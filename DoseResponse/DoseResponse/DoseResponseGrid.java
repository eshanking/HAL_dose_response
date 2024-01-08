package DoseResponse;

import HAL.Rand;
import HAL.Util;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import static HAL.Util.*;

import java.util.ArrayList;
// import java.util.Arrays;
import java.util.List;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import HAL.Tools.FileIO;
import HAL.Gui.UIGrid;
import HAL.Interfaces.SerializableModel;

class Source extends AgentSQ2Dunstackable<DoseResponseGrid>{
    // Blood vessel source
    double conc;
    double xPos;
    double yPos;

    public Source(double conc, int x, int y) {
        super();
        Init(conc, x, y);
    }

    void Init(double conc,int x, int y){
        this.conc=conc;
        this.xPos=x;
        this.yPos=y;
    }

    void SetSource(PDEGrid2D diff){
        diff.Set(xPos,yPos,conc);
    }
}

// class Cell extends AgentSQ2Dunstackable<DoseResponseGrid> {
//     int genotype;
//     double g_drugless;
//     double ic50;
//     double hillCoef;

//     public int GetGenotype(){
//         return this.genotype;
//     }

//     public double HillEqn(double conc){
//         // convert conc to log10
//         return this.g_drugless / (1 + Math.exp((this.ic50 - Math.log10(conc)) / this.hillCoef));
//     }

//     public double GetFitness(double conc){
//         return HillEqn(conc);
//     }
// }
    

public class DoseResponseGrid extends AgentGrid2D<DoseResponseCell> implements SerializableModel{

    public int nReplicates = 1;
    int nTSteps = 1000;
    double dt = 0.1; // time step in hours
    // diffusion grid
    PDEGrid2D diff;
    Source src;
    boolean constantGrid = false; // constant concentration of drug in the grid
    double constantGridConc = 0; // constant concentration of drug in the grid
    // double[] srcConc = new double[]{1000, 1000};
    double srcConc = 1000;
    int[] srcX = new int[]{50, 50};
    int[] srcY = new int[]{25, 75};
    int vesselSep = 50;
    Source[] srcList = new Source[srcX.length];
    double diffRate = 0.1;
    double diffBoundary = 0;
    double consumpRate = 0.5;
    // static grid parameters
    boolean staticGrid = false;
    double charLength = 1; // microns, distance at half max
    double drugConcScale = 1; // makes diffusion model more stable for lower concentrations
    int pulseInterval = 7*24;// in hours
    int pulseDuration = 24; // in hours
    double drugStopTime = nTSteps*dt + 1; // in hours
    boolean drugOn = false;
    boolean pulseDosing = false;
    // int currOnTime = 0; 

    // population initialization
    int xDim = 100;
    int yDim = 100; 
    int[] divHood = Util.VonNeumannHood(false);
    // public double initRadius = 10;
    public double initMutantProp = 0.1;
    String initGeometry = "circle"; // shape of initial population
    int initWidth = 100;
    double initDensity = 0.01;

    double dieProb = 0.1;
    double mutProb = 0.0001;
    
    // Synthetic parameters
    double[] GrowthRateList = new double[]{1.28949852, 1.14399848, 1.22802236, 0.93619847};
    double[] ic50List = new double[]{-0.49205992, 1.76224515,  1.39341393,  2.84653598};
    double[] hillCoefList = new double[]{-0.6824968,-0.6824968,-0.6824968,-0.6824968};
    double[] gminList = new double[]{0,0,0,0}; // placeholder
    
    // double[] GrowthRateList = new double[]{0.0393,0.0376,0.0360,0.0355};
    // double[] gminList = new double[]{0.0179,0.0286,0.0337,0.0296};
    // double[] hillCoefList = new double[]{2.18,2.86,17.45,1.85};
    // double[] ic50List = new double[]{0.0333,0.0371,0.020,0.0530};

    int n_genotype = 4;
    int nAllele = 2;
    int[] genotypeCounts = new int[n_genotype];
    int initGenotype = 0; // default wild-type genotype
    boolean useMaxConc = false; // if true, sets a max drug concentration above which the growth rate is zero
    double maxConc = 1; // maximum drug concentration

    // p.drugless_rates = [1.28949852, 1.14399848, 1.22802236, 0.93619847]
    // p.ic50 = [-0.49205992, 1.76224515,  1.39341393,  2.84653598]

    FileIO cellCountLogFile = null;
    String cellCountLogFileName = "./data/cellCountLog.csv";
    double logCellCountFrequency = 1;
    int tIdx = 0;
    boolean saveModelState = false;

    Rand rng = new Rand();
    int seed = 0;

    // Output - visualization
    UIGrid visCells;
    UIGrid visDiff;
    int scaleFactor = 2;
    boolean visualiseB = true;
    boolean saveFinalDiffImg = false;
    boolean saveFinalDiffGrid = false;
    boolean saveFinalPopGrid = false;
    String diffFileName = "~/tmp/diffGrid.csv";
    String popGridFileName = "~/tmp/popGrid.csv";
    public String imageOutDir = "./data/images/";
    int imageFrequency = -1;
    int pause = 1;
    int[] colorList = new int[]{RGB(0,1,0), RGB(0,0,1),
                                RGB(1,1,0), RGB(1,0,0)};

    boolean threeParamHill = false;

    // ------------------------------------------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------------------------------------------

    public DoseResponseGrid(int x, int y) {
        super(x, y, DoseResponseCell.class);
        this.xDim = x;
        this.yDim = y;
    }

    public DoseResponseGrid(){
        super(100,100,DoseResponseCell.class);
    }

    public DoseResponseGrid(int x, int y , double[] paramArr, double dt){
        super(x,y,DoseResponseCell.class);
        // SetParameters(paramArr);
        this.dt = dt;
    }

    @Override
    public void SetupConstructors() {
        _PassAgentConstructor(DoseResponseCell.class);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Model setup and configuration
    // ------------------------------------------------------------------------------------------------------------

    public void SetSaveParams(String diffFileName, String cellCountLogFileName, 
                              String popGridFileName, boolean saveFinalDiffGrid,
                              boolean saveFinalPopGrid){
        this.diffFileName = diffFileName;
        this.cellCountLogFileName = cellCountLogFileName;
        this.saveFinalDiffGrid = saveFinalDiffGrid;
        this.saveFinalPopGrid = saveFinalPopGrid;
        this.popGridFileName = popGridFileName;
    }

    public void SetDiffParams(double srcConc, int[] srcX, int[] srcY, 
                              double diffRate, double diffBoundary, double consumpRate){
        this.srcConc = srcConc;
        this.srcX = srcX;
        this.srcY = srcY;
        this.diffRate = diffRate;
        this.diffBoundary = diffBoundary;
        this.srcList = new Source[srcX.length];
        this.consumpRate = consumpRate;
        
    }

    public void SetDiffParams(double srcConc,double diffRate, double consumpRate, 
                              int vesselSep, boolean staticGrid, double charLength,
                              boolean constantGrid, double constantGridConc){
        this.srcConc = srcConc;
        this.diffRate = diffRate;
        this.srcList = new Source[srcX.length];
        this.consumpRate = consumpRate;
        this.vesselSep = vesselSep;
        this.srcY = new int[]{(int) Math.round(yDim/2.0 - vesselSep/2.0),
                            (int) Math.round(yDim/2.0 + vesselSep/2.0)};
        this.staticGrid = staticGrid;
        this.charLength = charLength;
        this.constantGrid = constantGrid;
        this.constantGridConc = constantGridConc;
        // System.out.println(srcY[0] + " " + srcY[1]);
    }


    public void SetDiffParams(double srcConc,double diffRate, double consumpRate, 
                              int vesselSep, boolean staticGrid, double charLength,
                              boolean constantGrid, double constantGridConc,
                              double drugConcScale){
        this.srcConc = srcConc;
        this.diffRate = diffRate;
        this.srcList = new Source[srcX.length];
        this.consumpRate = consumpRate;
        this.vesselSep = vesselSep;
        this.srcY = new int[]{(int) Math.round(yDim/2.0 - vesselSep/2.0),
                            (int) Math.round(yDim/2.0 + vesselSep/2.0)};
        this.staticGrid = staticGrid;
        this.charLength = charLength;
        this.constantGrid = constantGrid;
        this.constantGridConc = constantGridConc;
        this.drugConcScale = drugConcScale;
        // System.out.println(srcY[0] + " " + srcY[1]);
    }

    public void SetDiffParams(double srcConc,double diffRate, double consumpRate, 
                              int vesselSep, boolean staticGrid, double charLength,
                              boolean constantGrid, double constantGridConc,
                              double drugConcScale, boolean pulseDosing,
                              int pulseInterval, int pulseDuration, double drugStopTime){
        this.srcConc = srcConc;
        this.diffRate = diffRate;
        this.srcList = new Source[srcX.length];
        this.consumpRate = consumpRate;
        this.vesselSep = vesselSep;
        this.srcY = new int[]{(int) Math.round(yDim/2.0 - vesselSep/2.0),
                            (int) Math.round(yDim/2.0 + vesselSep/2.0)};
        this.staticGrid = staticGrid;
        this.charLength = charLength;
        this.constantGrid = constantGrid;
        this.constantGridConc = constantGridConc;
        this.drugConcScale = drugConcScale;
        this.pulseDosing = pulseDosing;
        this.pulseInterval = pulseInterval;
        this.pulseDuration = pulseDuration;
        this.drugStopTime = drugStopTime;
        // System.out.println(srcY[0] + " " + srcY[1]);
    }

    public void SetCellParams(double[] GrowthRateList, double[] ic50List, 
                              int n_genotype, int nAllele, double mutProb, double dieProb){
        this.GrowthRateList = GrowthRateList;
        this.ic50List = ic50List;
        // this.hillCoefList = hillCoefList;
        this.n_genotype = n_genotype;
        this.nAllele = nAllele;
        this.mutProb = mutProb;
        this.dieProb = dieProb;
    }

    public void SetCellParams(double mutProb, double dieProb){
        this.mutProb = mutProb;
        this.dieProb = dieProb;
    }

    public void SetCellParams(double mutProb, double dieProb, boolean threeParamHill,
                                boolean useMaxConc, double maxConc){
        this.mutProb = mutProb;
        this.dieProb = dieProb;
        this.threeParamHill = threeParamHill;
        this.useMaxConc = useMaxConc;
        this.maxConc = maxConc;
    }
    /*
     * Set the initial population parameters. initGeometry is either "circle" or "square".
     */
    public void SetInitPopParams(String initGeometry, int initWidth, 
                                 double initDensity, double initMutantProp,
                                 int initGenotype){
        this.initGeometry = initGeometry;
        // this.initRadius = initRadius;
        this.initWidth = initWidth;
        this.initDensity = initDensity;
        this.initMutantProp = initMutantProp;
        this.initGenotype = initGenotype;
    }

    public void SetSeed(int seed) {
        this.rng = new Rand(seed);
        // this.rn_ICs = new Rand(seed);
    }
    // public void ConfigureVisualisation(boolean visualiseB, int pause) {
    //     this.visualiseB = visualiseB;
    //     this.pause = pause;
    // }

    public void ConfigureImaging(boolean visualiseB, String imageOutDir, int imageFrequency,
                                 int scaleFactor, int pause, int[] colorList) {
        /*
        * Configure location of and frequency at which tumour is imaged.
        */
        this.visualiseB = visualiseB;
        this.imageOutDir = imageOutDir;
        this.imageFrequency = imageFrequency;
        this.scaleFactor = scaleFactor;
        this.pause = pause;
        this.colorList = colorList;
    }

    public void ConfigureImaging(boolean visualiseB, String imageOutDir, int imageFrequency,
                                 boolean saveFinalDiffImage) {
        /*
        * Configure location of and frequency at which tumour is imaged.
        */
        this.visualiseB = visualiseB;
        this.imageOutDir = imageOutDir;
        this.imageFrequency = imageFrequency;
        this.saveFinalDiffImg = saveFinalDiffImage;

    }

    public void ConfigureExperiment(int nReplicates, int nTSteps, double dt){
        this.nReplicates = nReplicates;
        this.nTSteps = nTSteps;
        this.dt = dt;
        this.drugStopTime = nTSteps*dt + 1;
    }

    // public void SetVesselDistance(int vesselDistance){
    //     // compute the y position of the two blood vessels
    //     this.srcY = new int[]{(int) Math.round(yDim/2.0 - vesselDistance/2.0),
    //                              (int) Math.round(yDim/2.0 + vesselDistance/2.0)};
    // }

    // ------------------------------------------------------------------------------------------------------------
    // Model initialization
    // ------------------------------------------------------------------------------------------------------------

    public double Distance(double x1,double y1,double x2,double y2){
        return Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    }

    public void InitPopulation(){
        if ("circle".equals(this.initGeometry)) {
            InitPopCircle(this.initWidth, this.initMutantProp);
        }
        else if ("square".equals(this.initGeometry)) {
            InitPopSquare(this.initWidth, this.initMutantProp, this.initDensity);
        }
    }

    public void InitPopCircle(double radius, double mutantProb) {
        //get a list of indices that fill a circle at the center of the grid
        int[] tumorNeighborhood = CircleHood(true, radius);
        int hoodSize = MapHood(tumorNeighborhood, xDim / 2, yDim / 2);
        for (int i = 0; i < hoodSize; i++) {
            DoseResponseCell c = NewAgentSQ(tumorNeighborhood[i]);
            if (rng.Double() < mutantProb) {
                int genotype = rng.Int(n_genotype-1) + 1; // any genotype other than 0
                // c.genotype = genotype;
                // c.g_drugless = GrowthRateList[genotype];
                // c.ic50 = ic50List[genotype];
                // c.hillCoef = hillCoefList[genotype];
                // c.gmin = gminList[genotype];
                c.setParams(genotype);
            } else {
                int genotype = initGenotype;
                // c.genotype = genotype;
                // c.g_drugless = GrowthRateList[genotype];
                // c.ic50 = ic50List[genotype];
                // c.hillCoef = hillCoefList[genotype];
                // c.gmin = gminList[genotype];
                c.setParams(genotype);
            }
            genotypeCounts[c.genotype] += 1;
        }
    }

    public void InitPopSquare(int x, double mutantProb, double density){
        
        int hoodSize = x*x;
        for (int i = 0; i < hoodSize; i++){
            if (rng.Double() < density){
                DoseResponseCell c = NewAgentSQ(i);
                if (rng.Double() < mutantProb) {
                    int genotype = rng.Int(n_genotype-1) + 1; // any genotype other than 0
                    // c.genotype = genotype;
                    // c.g_drugless = GrowthRateList[genotype];
                    // c.ic50 = ic50List[genotype];
                    // // c.hillCoef = hillCoef;
                    // c.hillCoef = hillCoefList[genotype];
                    // c.gmin = gminList[genotype];
                    c.setParams(genotype);
                } else {
                    int genotype = initGenotype;
                    // c.genotype = genotype;
                    // c.g_drugless = GrowthRateList[genotype];
                    // c.ic50 = ic50List[genotype];
                    // c.hillCoef = hillCoefList[genotype];
                    // c.gmin = gminList[genotype];
                    c.setParams(genotype);
                }
                genotypeCounts[c.genotype] += 1;
            }
            
        }
    }

    public void InitPDEGrid() {
        diff = new PDEGrid2D(xDim, yDim);
        srcList = new Source[srcX.length];
        // initialize sources
        for (int i = 0; i < srcX.length; i++) {
            Source src = new Source(srcConc, srcX[i], srcY[i]);
            srcList[i] = src;
        }

        if (staticGrid) {
            // set the initial concentration
            // calculate the static profile for each source
            // initialize zero field

            double conc;
            double dist;
            Source src;

            if (constantGrid){
                for (int i = 0; i < diff.length; i++){
                    diff.Set(i, constantGridConc);
                }
                diff.Update();
            }
            else{
                for (int i = 0; i < diff.length; i++){
                    conc = 0;
                    for (int s = 0; s < srcList.length; s++){
                        src = srcList[s];
                        dist = Distance(i / xDim, i % yDim, src.xPos, src.yPos);
                        conc += src.conc * Math.exp(-dist * Math.log(2) / charLength);
                    }
                    diff.Set(i, conc);
                }
                diff.Update();
            }

            // for (int i = 0;i < diff.length; i++){
            //     conc = 0;
            //     for (int s = 0; s < srcList.length; s++){
            //         src = srcList[s];
            //         dist = Distance(i / xDim, i % yDim, src.xPos, src.yPos);
            //         conc += src.conc * Math.exp(-dist * Math.log(2) / charLength);
            //     }
            //     diff.Set(i, conc);
            // }
            // diff.Update();
        }
        
    }

    // ------------------------------------------------------------------------------------------------------------
    // Model step
    // ------------------------------------------------------------------------------------------------------------

    public void StepPDEGrid() {



        if (!staticGrid){
            for (Source src : srcList) {
                src.SetSource(diff); // sets the concentration at the source's x,y position
            }
            diff.Diffusion(diffRate,diffBoundary);
            // field consumption
            // for (int i = 0; i < diff.length; i++) {
            //     diff.Mul(i, -consumpRate);
            // }
            diff.MulAll(-consumpRate);
            diff.Update();
        }

        else if (pulseDosing){

            if (drugStopTime < tIdx*dt){ // if the treatment has ended
                if (drugOn){ // if the drug is currently on, turn it off
                    for (int i = 0; i < diff.length; i++){
                        diff.Set(i, 0);
                    }
                    diff.Update();
                }
            }
            
            else{ 
                boolean chgField = false;
                if (tIdx*dt % pulseInterval < pulseDuration){
                    // drugOn = true;
                    if (!drugOn){
                        chgField = true;
                        drugOn = true;
                    }
                }
                else{
                    if (drugOn){
                        chgField = true;
                        drugOn = false;
                    }
                    // drugOn = false;
                }
                if (chgField){ // change the state of the PDE grid
                    // if drugOn is false, then set the concentration to zero
                    if (!drugOn){
                        for (int i = 0; i < diff.length; i++){
                            diff.Set(i, 0);
                        }
                        diff.Update();
                    }
                    else{ // drugOn is true, set the concentration to the static profile
                        double conc;
                        double dist;
                        for (int i = 0; i < diff.length; i++){
                            conc = 0;
                            for (int s = 0; s < srcList.length; s++){
                                src = srcList[s];
                                dist = Distance(i / xDim, i % yDim, src.xPos, src.yPos);
                                conc += src.conc * Math.exp(-dist * Math.log(2) / charLength);
                            }
                            diff.Set(i, conc);
                        }
                        diff.Update();
                    }
                }
            }
        }

    }

    public List<Integer> genNeighbors(int genotype, int nAllele) {
        List<Integer> neighbors = new ArrayList<>();
        int mask = (1 << nAllele) - 1;

        for (int m = 0; m < nAllele; m++) {
            int neighbor = genotype ^ (1 << m);
            neighbors.add(neighbor & mask);
        }

        return neighbors;
    }

    public void StepModel(){
        StepCells();
        StepPDEGrid();
    }

    public void StepCells() {
        int currPos;

        for (DoseResponseCell cell : this) {
            currPos = cell.Isq();
            // cell.StepCell(dieProb);
            if (rng.Double() < dieProb * dt){
                cell.Dispose();
                genotypeCounts[cell.genotype] -= 1;
            }
            else{
                // get the drug concentration
                double currConc = diff.Get(currPos);
                double divProb = cell.GetFitness(currConc) * dt;

                if (rng.Double() < divProb) {
                    int options = MapEmptyHood(divHood, cell.Xsq(), cell.Ysq());
                    if (options > 0) {

                        int iDaughter = divHood[rng.Int(options)];
                        DoseResponseCell daughter = NewAgentSQ(iDaughter);
                        
                        if (rng.Double() < mutProb * dt){ // should be multiplied by dt?
                            // mutate
                            List<Integer> neighbors = genNeighbors(cell.genotype, nAllele);
                            int newGenotype = neighbors.get(rng.Int(neighbors.size()));
                            // daughter.genotype = newGenotype;
                            // daughter.g_drugless = GrowthRateList[newGenotype];
                            // daughter.ic50 = ic50List[newGenotype];
                            // daughter.hillCoef = hillCoefList[newGenotype];
                            // daughter.gmin = gminList[newGenotype];
                            daughter.setParams(newGenotype);
                            
                        }
                        else{
                            // inherit properties from parent
                            // daughter.genotype = cell.genotype;
                            // daughter.g_drugless = cell.g_drugless;
                            // daughter.ic50 = cell.ic50;
                            // daughter.hillCoef = cell.hillCoef;
                            // daughter.gmin = cell.gmin;
                            daughter.setParams(cell.genotype);
                            
                        }
                        genotypeCounts[daughter.genotype] += 1;
                    }
                }
            }
        }
    }
    
    public static void main(String[] args) {

        int x=100;
        int y=100;
        int timesteps=100;

        GridWindow visCells=new GridWindow(x,y,5,true,null,false);
        GridWindow visDiff=new GridWindow(x,y,5,true,null,true);
    
        DoseResponseGrid model = new DoseResponseGrid(x,y);
        model.consumpRate = 0.0001;
        model.initGeometry = "square";
        // model.initGeometry = "circle";
        // model.initWidth = 10;
        model.initDensity = 0.01;
        model.threeParamHill = false;
        model.srcConc = 10000;
        model.mutProb = 0.001;
        model.dieProb = 0.02;
        model.drugConcScale = 0.01;
        model.dt = 1;
        model.staticGrid = true;
        model.pulseDosing = true;
        model.pulseInterval = 24;
        model.pulseDuration = 10;
        
        // model.constantGrid = true;
        // model.constantGridConc =  0;
        model.InitPopulation();
        model.InitPDEGrid();
    
        for (int i = 0; i < timesteps; i++) {
            // model step
            model.StepModel();
            model.Draw(visCells,visDiff);

            System.out.println("Timestep: " + i*model.dt + " ");
            System.out.println(model.drugOn);
            // tIdx*dt % pulseInterval < pulseDuration

            System.out.println('\n');
            model.tIdx++;

        }
        // model.SaveDiffGrid(model.diff, "/Users/eshanking/repos/HAL_fds/DoseResponse/data/diffGrid.csv");
    }

    public void Run() throws IOException {
        // System.out.println(threeParamHill);
        UIGrid currVis = new UIGrid(xDim, yDim, scaleFactor, visualiseB); // For head-less run
        this.visCells = currVis;

        Boolean completedSimulationB = false;
        Boolean logged = false;
        // System.out.println(this.initGeometry);
        InitPopulation();
        // initialize the diffusion grid
        InitPDEGrid();
        if (cellCountLogFile==null && cellCountLogFileName!=null) {InitialiseCellLog(this.cellCountLogFileName);}

        SaveCurrentCellCount(0);
        SaveTumourImage(tIdx);
        
        completedSimulationB = false;
        while (!completedSimulationB) {
            // currVis.TickPause(pause);
            Draw(currVis);
            StepModel();
            logged = SaveCurrentCellCount(tIdx);
            SaveTumourImage(tIdx);
            tIdx++;
            completedSimulationB = (tIdx>nTSteps)?true:false;
        }

        if (saveFinalDiffImg){
            SaveDiffImage(tIdx);
        }
        if (saveFinalDiffGrid){
            SaveDiffGrid(diffFileName);
        }

        if (saveFinalPopGrid){
            SavePopGrid(popGridFileName);
        }
        this.Close(logged);
    }
    // ------------------------------------------------------------------------------------------------------------
    // Manage and save output
    // ------------------------------------------------------------------------------------------------------------

    public void InitialiseCellLog(String cellCountLogFileName) {
        cellCountLogFile = new FileIO(cellCountLogFileName, "w");
        WriteLogFileHeader();
        this.cellCountLogFileName = cellCountLogFileName;
        this.logCellCountFrequency = 1;
    }

    private void WriteLogFileHeader() {
        // cellCountLogFile.Write("TIdx,Time,NCells_S,NCells_R,NCells,DrugConcentration,rS,rR,mS,mR,dS,dR,dD_div_S,dD_div_R,dt");

        // cellCountLogFile.Write("TIdx,Time,NCells_S,NCells_R,NCells");
        cellCountLogFile.Write("TIdx,Time,Pop,n_gen0,n_gen1,n_gen2,n_gen3,mutProb,dieProb,diffRate,consumpRate,srcConc,src1_X,src2_X,src1_Y,src2_Y,drugOn");
        
        cellCountLogFile.Write("\n");
    }
    public void Close() {
        if (cellCountLogFile!=null) {cellCountLogFile.Close();}
    }

    public Boolean SaveCurrentCellCount(int currTimeIdx) {
        Boolean successfulLog = false;
        if ((currTimeIdx % (int) (logCellCountFrequency)) == 0 && logCellCountFrequency > 0) {
            cellCountLogFile.WriteDelimit(GetModelState(),",");
            // if (extraSimulationInfoNames!=null) {
            //     cellCountLogFile.Write(",");
            //     cellCountLogFile.WriteDelimit(extraSimulationInfo, ",");
            // }
            cellCountLogFile.Write("\n");
            successfulLog = true;
        }
        return successfulLog;
    }

    public double[] GetModelState() {

        double drugOnStatus = 0;
        if (drugOn){
            drugOnStatus = 1;
        }
   
        return new double[] {tIdx, tIdx*dt, Pop(), genotypeCounts[0], genotypeCounts[1], 
                             genotypeCounts[2], genotypeCounts[3], mutProb, dieProb, diffRate,
                             consumpRate, srcConc, srcX[0], srcX[1], srcY[0], srcY[1], drugOnStatus};
        
    }

    public void SaveTumourImage(int currTimeIdx) {
        if (imageFrequency > 0 && (currTimeIdx % (int) (imageFrequency/dt)) == 0) {
            this.visCells.ToPNG(imageOutDir +"/img_t_"+currTimeIdx*dt+".png");
        }
    }

    public void SaveDiffImage(int tIdx) {

        if (this.visDiff == null){
            this.visDiff = new UIGrid(xDim, yDim, scaleFactor, visualiseB);
        }
        Draw(this.visCells,this.visDiff);

        this.visDiff.ToPNG(imageOutDir +"/diffGrid_"+tIdx*dt+".png");

    }

    public void SaveDiffGrid(String diffFileName) throws IOException {
        try{
            File diffGridFile = new File(diffFileName);
            diffGridFile.createNewFile();
            FileWriter diffGridFileWriter = new FileWriter(diffFileName);

            if (tIdx*dt > drugStopTime){
                // reinitialize the field
                InitPDEGrid();
            }
            double[] field = diff.GetField();
            
            diffGridFileWriter.write("x,y,field\n");
            // int lineNum = 0;
            for (int i = 0; i < field.length; i++) {
                diffGridFileWriter.write(i % xDim + "," + i / xDim + "," + field[i] + "\n");
                // lineNum++;
            }
            diffGridFileWriter.close();
            // System.out.println("Wrote "+lineNum+" lines to "+diffFileName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*
     * Save the population grid to a file
     */
    public void SavePopGrid(String popFileName) throws IOException {
        try{
            File popGridFile = new File(popFileName);
            popGridFile.createNewFile();
            FileWriter popGridFileWriter = new FileWriter(popFileName);
            
            popGridFileWriter.write("x,y,genotype\n");
            int currPos;
            for (DoseResponseCell cell : this){
                currPos = cell.Isq();
                popGridFileWriter.write(currPos % xDim + "," + currPos / xDim + "," + cell.GetGenotype() + "\n");
            }
            popGridFileWriter.close();
            // System.out.println("Wrote "+lineNum+" lines to "+diffFileName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // public void SaveDiffGrid(String diffFileName) {
    //     // get the diffusion grid values
    //     double[] field = diff.GetField();
    //     // new fileio object
    //     FileIO diffGridFile = new FileIO(diffFileName, "w");
    //     // write the header
    //     diffGridFile.Write("x,y,field\n");
    //     // write the values
    //     int lineNum = 0;
    //     for (int i = 0; i < field.length; i++) {
    //         System.out.println(i % xDim + "," + i / xDim + "," + field[i]);
    //         diffGridFile.Write(i % xDim + "," + i / xDim + "," + field[i] + "\n");
    //         lineNum++;
    //     }
    //     System.out.println("Wrote "+lineNum+" lines to "+diffFileName);

    // }
    public void SaveDiffGrid(PDEGrid2D diffGrid, String diffFileName) {
        double[] field = diffGrid.GetField();
        FileIO diffGridFile = new FileIO(diffFileName, "w");
        diffGridFile.Write("x,y,field\n");
        
        for (int i = 0; i < field.length; i++) {
            diffGridFile.Write(i % xDim + "," + i / xDim + "," + field[i] + "\n");
        }
    }

    public void Close(Boolean logged) {
        if (!logged) {
            tIdx--;
            SaveCurrentCellCount(0);}
        if (cellCountLogFile!=null) {cellCountLogFile.Close();}
    }
    public void SaveModelState(String stateFileName) {
        // Can't have active pointers when saving the model. So, close everything here.
        if (cellCountLogFile!=null) {
            cellCountLogFile.Close();
            cellCountLogFile = null;
        }
        if (visCells!=null) {visCells = null;}
        SaveState(this,stateFileName);
    }
    /*
     * Draw cells and diffusion grid
     */
    public void Draw(UIGrid visCells, UIGrid visDiff){
        for (DoseResponseCell cell : this) {
            visCells.SetPix(cell.Isq(),colorList[cell.genotype]);//draw sources and sinks
        }
        for (int i = 0; i < length; i++) {//length of the Grid
            //visDiff.SetPix(i,SetAlpha(HeatMapRGB(diff.Get(i)*4),diff.Get(i)*4));
            visDiff.SetPix(i,HeatMapRGB(diff.Get(i)*4));
        }
    }
    /*
     * Draw cells only
     */
    public void Draw(UIGrid visCells){
        for (DoseResponseCell cell : this) {
            visCells.SetPix(cell.Isq(),colorList[cell.genotype]);//draw sources and sinks
        }
    }
}
